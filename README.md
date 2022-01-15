# Patient stratification for simultaneous learning of molecular subtype and survival outcome using genetic algorithm-based gene set selection

Motivation: Patient stratification is a clinically important task because it allows to establish and develop efficient treatment strategies for particular groups of patients. Molecular subtypes have been successfully defined using transcriptomic profiles and they are used effectively in clinical practice, e.g., PAM50 subtypes of breast cancer. Survival prediction using transcriptome data contributed to predicting prognosis of disease and also identifying genes related to prognosis. It is desirable to stratify patients considering these two aspects simultaneously. However, there is no existing methods for patient stratification that consider molecular subtypes and survival outcomes at the same time.

Results: We devised a methodology that stratify patients considering molecular subtypes and survival outcomes simultaneously. A set of genes was selected from transcriptome data by using genetic algorithm, and their expression quantities were utilized to represent each patient as a risk score. The patients are ordered and stratified according to the score. This methodology was applied to breast cancer patients (TCGA-BRCA), and validated in an independent cohort (SCAN-B). In this experiment, our method was successful in stratifying patients in accordance with both molecular subtype and survival outcome. The robustness of the results was shown through repeated experiments, and genes related to prognosis of patients were successfully identified. Additionally, it was observed that the risk score can be used to evaluate the molecular aggressiveness of individual patients.

![image](https://github.com/BonilKoo/patient_stratification/blob/main/overview.png)

## Data preparation
Three types of data are required.
1. Gene expression profile
2. Subtype information
3. Survival data

## Usage
```
python patient_stratification.py [options] \\
                                 --ourdir <output directory> \\
                                 --expression <gene expression profile> \\
                                 --subtype <subtype information file> \\
                                 --survival <survival data file> \\
                                 --subtype_order <string>
```

### Package requirements
+ numpy
+ pandas
+ scipy
+ scikit-learn
+ lifelines (https://github.com/CamDavidsonPilon/lifelines)
+ parmap (https://github.com/zeehio/parmap)

### Mandatory arguments:
--outdir <path>
: output directory to save results. If the path doen not exist, it is automatically created.
  
--expression <file>
: gene expression data of samples. Each row is a gene and each columns is a sample. (Seperated by tab)
  | Gene  | Sample 1 | Sample2 | Sample3 | ... |
  | --- | --- | --- | --- | --- |
  | **Gene 1**  | **5.45**  | **3.80**  | **4.47**  |   |
  | **Gene 2**  | **19.35**  | **23.21**  | **21.20**  |   |
  | **Gene 3**  | **42.57**  | **57.78**  | **58.28**  |   |
  | **...**  |   |   |   |   |
  
--subtype <file>
: subtype information for samples. The first column is sample, and the second column is subtype of the sample. (No header and separted by tab)
  | Sample 1  | Subtype 1  |
  | --- | --- |
  | **Sample 2**  | **Subtype 3**  |
  | **Sample 3**  | **Subtype 2**  |
  | **...**  |   |
  
--survival <file>
: survival data for samples. The first columns is sample, the second column is event, and the third column is duration. Event is 1 if it occured, 0 if it didn't. (No header and separted by tab)
  | Sample 1  | 0  | 4047  |
  | --- | --- | --- |
  | **Sample 2**  | **1**  | **4005**  |
  | **Sample 3**  | **0**  | **1474**  |
  | **...**  |   |   |
  
--subtype_order <string>
: subtype order to be sorted and stratified. (Separated by comma)
  E.g., for breast cancer, "LumA,LumB,Her2,Basal".

### Optional arguments:
--survival_coef <float> (default=0.5)
: A coefficient to modulate the balance between the effects of the two scores on the fitness value.
  
--n_subpop <int> (default=10)
: Number of subpopulations.
  
--n_chrom <int> (default=150)
: Number of chromosomes in each subpopulation.
  
--n_init_gene <int> (default=10)
: Number of genes initially selected.
  
--elite_proportion <float> (default=0.1)
: Proportion of chromosomes to be selected as elite in each subpopulation in the selection step of genetic algorithm.
  
--num_participant <int> (default=16)
: Number of chromosomes participating in K-way tournament selection.
  
--crossover_prob <float> (default=0.3)
: Crossover probability in uniform crossover.
  
--mutation_prob <float> (default=0.0005)
: Mutation probability.
  
--patience <int> (default=10)
: Number of iterations with no improvement after which run will be terminated.
  
--seed <int> (default=42)
: Random seed.
