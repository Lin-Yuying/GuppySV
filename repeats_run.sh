ref=$1
prefix=$2
output=$3
genomesize=$4
################## 1. RepeatModeler #######################
### 1.1 Build Database
/Users/evolutioneco/Downloads/Software/RepeatModeler-2.0.2a/BuildDatabase -name ${prefix} ${ref} 

### 1.2 Run RepeatModeler, include LTR
/Users/evolutioneco/Downloads/Software/RepeatModeler-2.0.2a/RepeatModeler -database ${prefix} -pa 10 -LTRStruct > ${prefix}_run.out &

################## 2. RepeatMasker #######################
### 2.1 Running repeat masker, using rmblast for searching 
~/Downloads/Software/RepeatMasker/RepeatMasker -lib consensi.fa.classified -dir ${output} -s -a -nolow -html -gff ${ref}

### 2.2 align2divsum 
~/Downloads/Software/RepeatMasker/util/calcDivergenceFromAlign.pl -s ${prefix}.divsum ${prefix}.align

### 2.3 Repeat landscape
~/Downloads/Software/RepeatMasker/util/createRepeatLandscape.pl -div ${prefix}.divsum -g ${genomesize} > ${prefix}.html
