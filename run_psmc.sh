#bcftools=~/bin/bcftools-1.16/bcftools                                                                                                                                                                      
SAMPLE=$1                                                                                     
#AVGDEP=$2                                                                                                              
#minDep=$(expr $AVGDEP / 3) # 1/3 of avg depth                                                                          
#maxDep=$(expr $AVGDEP \* 2) # twice of avg depth    

#consensus sequences                                                                                                    
#${bcftools} mpileup -C50 -f ${REF} ./1.alignment_male/sortdedup_${SAMPLE}.bam | ${bcftools} call -c -| vcfutils.pl vcf2fq -d ${minDep} -D ${maxDep} | gzip > ${SAMPLE}.diploid.fq.gz                                                           
                                                                                                                        
utils/fq2psmcfa -q20 ${SAMPLE}.diploid.fq.gz > ${SAMPLE}.diploid.psmcfa                                                 
/home/linyuying/bin/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ${SAMPLE}.diploid.psmc ${SAMPLE}.diploid.psmcfa          
utils/psmc2history.pl ${SAMPLE}.diploid.psmc | utils/history2ms.pl > ms-cmd.sh                                          
utils/psmc_plot.pl diploid ${SAMPLE}.diploid.psmc
