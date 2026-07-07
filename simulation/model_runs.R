library(tidyverse)
set_theme(theme_classic() + theme(strip.background = element_rect(colour = NA, fill = 'grey90')))
library(furrr)

run_model <- function(k, seed, slim_script = 'time_to_loss.slim') {
  out <- system2('slim',
    c(str_glue('-d "k={k}"'), '-l 0', '-x', str_glue('-s {seed}'), slim_script),
    stdout = TRUE
  )
  out <- str_split(out, ' ')[[1]] |> map_int(as.integer) |> sort()
  data.frame(two = out[1], one = out[2], none = out[3])
}

get_population_proportions <- function(result, generation) {
  result |>
    mutate(polymorphic_in = case_when(
      generation < two ~ 'three',
      generation < one ~ 'two',
      generation < none ~ 'one',
      TRUE ~ 'none'
    )) |>
    count(polymorphic_in) |>
    mutate(
      frac = n / sum(n),
      generation = generation
    )
}

plan(multisession, workers = 24)
seeds <- sample(1e7, 1000)

cat('Running 1000 runs of k = 20...\n')
result_k20 <- future_map_dfr(
  seeds,
  \(s) run_model(k = 20, seed = s),
  .progress = TRUE, .options = furrr_options(chunk_size = 1)
)
props_k20 <- future_map_dfr(seq(1, 1e6, 100), get_population_proportions, result = result_k20)

cat('Running 250 runs of k = 10...\n')
result_k10 <- future_map_dfr(
  seeds[1:250],
  \(s) run_model(k = 10, seed = s),
  .progress = TRUE, .options = furrr_options(chunk_size = 1)
)
props_k10 <- future_map_dfr(seq(1, 1e6, 100), get_population_proportions, result = result_k10)

cat('Running 100 runs of k = 5...\n')
result_k5 <- future_map_dfr(
  seeds[1:100],
  \(s) run_model(k = 5, seed = s),
  .progress = TRUE, .options = furrr_options(chunk_size = 1)
)
props_k5 <- future_map_dfr(seq(1, 1e6, 100), get_population_proportions, result = result_k5)
plan(sequential)

result <- bind_rows(k20 = result_k20, k10 = result_k10, k5 = result_k5, .id = 'k')
props <- bind_rows(k20 = props_k20, k10 = props_k10, k5 = props_k5, .id = 'k')

data.table::fwrite(result, 'simulation_results.csv.gz')
data.table::fwrite(props, 'simulation_results_as_proportions.csv.gz')

ggplot(result) +
  stat_ecdf(aes(two, color = 'Polymorphic in two populations')) +
  stat_ecdf(aes(one, color = 'Polymorphic in one populations')) +
  stat_ecdf(aes(none, color = 'Polymorphic in no populations')) +
  scale_x_continuous(
    limits = c(0, 1e6), labels =  \(x) x / 1e3, name = 'Generations (x1000)', expand = c(0, 0),
    sec.axis = sec_axis(\(x) x * .6, labels =  \(x) x / 1e3, name = 'Time in years (x1000)')
  ) +
  scale_y_continuous(expand = c(0, 0), name = 'Proportion of simulation runs') +
  facet_grid(cols = vars(k))

ggplot(props, aes(generation, frac, fill = factor(polymorphic_in, c('three', 'two', 'one', 'none')))) +
  geom_col(width = 100) +
  scale_x_continuous(
    limits = c(0, 1e6), labels =  \(x) x / 1e3, name = 'Generations (x1000)', expand = c(0, 0),
    sec.axis = sec_axis(\(x) x * .6, labels =  \(x) x / 1e3, name = 'Time in years (x1000)')
  ) +
  scale_fill_viridis_d(
    name = 'Polymorphic in ...',
    limits = c('three', 'two', 'one', 'none'),
    labels = c('three populations', 'two populations', 'one population', 'none'),
    direction = -1
  ) +
  scale_y_continuous(expand = c(0, 0), name = 'Proportion of simulation runs', breaks = c(0.25, 0.5, 0.75)) +
  facet_grid(rows = vars(
    k |>
      factor(c('k20', 'k10', 'k5')) |>
      fct_recode('1000 runs with k = 20' = 'k20', '250 runs with k = 10' = 'k10', '100 runs with k = 5' = 'k5')
  ))

# k=20 only
props |>
  filter(k == 'k20') |>
  ggplot(aes(generation, frac, fill = factor(polymorphic_in, c('three', 'two', 'one', 'none')))) +
  geom_col(width = 100) +
  scale_x_continuous(
    limits = c(0, 1e6), labels =  \(x) x / 1e3, name = 'Generations (x1000)', expand = c(0, 0),
    sec.axis = sec_axis(\(x) x * .6, labels =  \(x) x / 1e3, name = 'Time in years (x1000)')
  ) +
  scale_fill_viridis_d(
    name = 'Polymorphic in ...',
    limits = c('three', 'two', 'one', 'none'),
    labels = c('three populations', 'two populations', 'one population', 'none'),
    direction = -1
  ) +
  scale_y_continuous(expand = c(0, 0), name = 'Proportion of simulation runs', breaks = c(0.25, 0.5, 0.75))

ggsave('k20_simulation_proportions.pdf', w = 6, h = 6)
