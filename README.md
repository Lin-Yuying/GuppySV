# GuppySV

Scripts for Lin, Y., van der Bijl, W., Mank, J. E. Sexual selection maintains ancient structural polymorphisms across natural populations. 2026, Under Review

detailed pipeline to be done

1. Genotyping and SNP filtering with [female reference genome](http://uswest.ensembl.org/Poecilia_reticulata/Info/Index) using [BWA MEM](https://github.com/lh3/bwa), [BCFtools](https://github.com/samtools/bcftools) and [VCFtools](https://vcftools.github.io/index.html)
   ```
   sh 01.genotyping.sh
   ```
  
2. Running localPCA using [lostruct](https://github.com/petrelharp/local_pca)
   ```
   # plot all pops 
   Rscript 02.1.plot_lostrcut.R
   # or do one population one chromosome for a single run
   Rscript 02.2.single_lostruct.R [pop] [chr]
   ```
   
3. PCA and admixture using [PLINK2.0](https://www.cog-genomics.org/plink/2.0/)
   ```
   sh 03.run_pca_admixture.sh 
   ```
   
4. Inferring demogrphic history using [PSMC](https://github.com/lh3/psmc)
   ```
   sh 04.run_psmc.sh
   ```
   
5. Calculating Fst and nucleotide diversity using [VCFtools](https://vcftools.github.io/index.html)
   ```
   sh 05.calc_fst_pi.sh 
   ```

7. Calculating SV age using [pixy](https://github.com/ksamuk/pixy)
   ```
   sh cal_inv_age.sh
   ```
