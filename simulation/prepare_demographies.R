library(tidyverse)

demo <- list.files('demographies/', '.txt*', full.names = TRUE) |>
  map_dfr(data.table::fread) |>
  set_names(c('individual', 'population', 'time', 'Ne')) |>
  mutate(
    Ne = Ne * 10 ^ 4,
    generation = time / 0.6
  ) |> as_tibble()

ggplot(demo, aes(time, Ne, color = population, group = individual)) +
  geom_step() +
  #scale_x_log10() +
  facet_wrap(vars(population)) +
  coord_cartesian(xlim = c(1e2, NA), ylim = c(0, 20e4))

max(demo$time)

max_generation <- 1e6
step_size <- 1e4       # update pop size every step_size generations
min_generation <- 1e4  # flatten anything more recent than this (years)

pairs <- tibble::tribble(
  ~population, ~pair,
  "ARH","AR", "ARL","AR",
  "QUH","QU", "QUL","QU",
  "YAH","YA", "YAL","YA"
)

medlog <- function(x) exp(median(log(x), na.rm = TRUE))

windows <- tibble(window = seq(1, max_generation, step_size))

# 1) per-population summary (across individuals)
pop_curve <- demo |>
  filter(generation > min_generation) |>
  mutate(window = floor(generation / step_size) * step_size) |>
  group_by(individual, population) |>
  complete(window = seq(0, max_generation, by = step_size)) |>
  fill(Ne, .direction = "downup") |>
  ungroup() |>
  summarise(Ne = medlog(Ne), .by = c(population, window)) |>
  left_join(pairs, join_by(population)) |>
  summarise(Ne = sum(Ne), .by = c(pair, window)) |>
  group_by(pair) |>
  complete(window = seq(0, max_generation, by = step_size)) |>
  fill(Ne, .direction = "downup") |>
  mutate(
    generation = max_generation - window,
    generation = ifelse(generation == 0, 1, generation),
  ) |>
  arrange(pair, generation) |>
  select(population = pair, generation, Ne) |>
  ungroup() |>
  mutate(population = as.numeric(as.factor(population)))

ggplot(pop_curve, aes(generation, Ne, color = as.factor(population))) + geom_step() +
  ylim(0, NA)

data.table::fwrite(pop_curve, "simple_demography.csv")
