# GuppySV

Scripts for Lin, Y., van der Bijl, W., Mank, J. E. Sexual selection maintains ancient structural polymorphisms across natural populations. 2026, Under Review

detailed pipeline to be done

1. Genotyping and SNP filtering
  ```
  sh 01.genotyping.sh
  ```
2. Local PCA using lostruct
   ```
   # plot all pops 
   Rscript 02.1.plot_lostrcut.R
   # or do one population one chromsome for a single run
   Rscript 02.2.single_lostruct.R [pop] [chr]
   ```
4. PCA and admixture
   ```
   sh 03.run_PCA_Admixture.sh 
   ```
