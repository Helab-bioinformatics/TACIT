#makeTagDirectory HiCTagDirectory/ GSE82185_late_2cell_rep1234_allValidPairs.txt -format HiCsummary
analyzeHiC HiCTagDirectory -res 5000 -pvalue 0.05 -nomatrix -cpu 25 -washu -interactions late2cell_interactions.txt


makeTagDirectory 8cell_rep1 GSM2203837_8cell_rep1_allValidPairs.txt -format HiCsummary
analyzeHiC 8cell_rep1 -res 5000 -pvalue 0.05 -nomatrix -cpu 25 -washu -interactions 8cell_rep1_interactions.txt 

makeTagDirectory early2cell GSE82185_early_2cell_rep123_allValidPairs.txt -format HiCsummary
analyzeHiC early2cell -res 5000 -pvalue 0.05 -nomatrix -cpu 25 -washu  -interactions early2cell_interactions.txt 


makeTagDirectory 8cell_rep2 GSM2203838_8cell_rep2_allValidPairs.txt -format HiCsummary
analyzeHiC 8cell_rep2 -res 5000 -pvalue 0.05 -nomatrix -cpu 25 -washu -interactions 8cell_rep2_interactions.txt 
