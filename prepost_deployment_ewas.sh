#User script for running pre-post deployment EWAS

qsub -t1-8 -e errandout/ -o errandout/ PREPOST_runjobs.pbs

echo 1 CpG Beta SE z p > header.txt
cat header.txt data_analysis/r_assoc/all_studymax_condv0_PTSDbroad_agecelldrri_HC3_aug8_2017_*.model | awk '{print $2,$3,$4,$5,$6}' >  /home/genetics/Desktop/new_methylation/data_analysis/RESULTS/all_studymax_condv0_PTSDbroad_agecelldrri_HC3_aug8_2017.mwas
awk '{print $5}' /home/genetics/Desktop/new_methylation/data_analysis/RESULTS/all_studymax_condv0_PTSDbroad_agecelldrri_HC3_aug8_2017.mwas > /home/genetics/Desktop/new_methylation/data_analysis/RESULTS/all_studymax_condv0_PTSDbroad_agecelldrri_HC3_aug8_2017.pval

Rscript qq_plot.R data_analysis/RESULTS/all_studymax_condv0_PTSDbroad_agecelldrri_HC3_aug8_2017.pval data_analysis/RESULTS/all_studymax_condv0_PTSDbroad_agecelldrri_HC3_aug8_2017.qq 1
