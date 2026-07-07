# GuppySV

Scripts for Lin, Y., van der Bijl, W., Mank, J. E. Sexual selection maintains ancient structural polymorphisms across natural populations. 2026, Under Review

detailed pipeline to be done

1. Genotyping and SNP filtering using BCFtools and VCFtools
   ```
   sh 01.genotyping.sh
   ```
  
2. Running localPCA using lostruct
   ```
   # plot all pops 
   Rscript 02.1.plot_lostrcut.R
   # or do one population one chromosome for a single run
   Rscript 02.2.single_lostruct.R [pop] [chr]
   ```
   
3. PCA and admixture using PLINK2 
   ```
   sh 03.run_pca_admixture.sh 
   ```
   
4. Inferring demogrphic history using PSMC
   ```
   sh 04.run_psmc.sh
   ```
   
5. Calculating Fst and nucleotide diversity using VCFtools
   ```
   sh 05.calc_fst_pi.sh.
   ```

7. Calculating SV age
   ```
   sh 06.calc_sv_age.sh
   ```
