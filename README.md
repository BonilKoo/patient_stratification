# Risk Stratification for Breast Cancer Patient by Simultaneous Learning of Molecular Subtype and Survival Outcome Using Genetic Algorithm-Based Gene Set Selection

Patient stratification is a clinically important task because it allows us to establish and develop efficient treatment strategies for particular groups of patients. Molecular subtypes have been successfully defined using transcriptomic profiles and they are used effectively in clinical practice, e.g., PAM50 subtypes of breast cancer. Survival prediction contributed to understanding diseases and also identifying genes related to prognosis. It is desirable to stratify patients considering these two aspects simultaneously. However, there are no methods for patient stratification that consider molecular subtypes and survival outcomes at once. Here, we propose a methodology to deal with the problem. Genetic algorithm is used to select a gene set from transcriptome data, and their expression quantities are utilized to assign a risk score to each patient. The patients are ordered and stratified according to the score. A gene set was selected by our method on a breast cancer cohort (TCGA-BRCA), and we examined its clinical utility using an independent cohort (SCAN-B). In this experiment, our method was successful in stratifying patients with respect to both molecular subtype and survival outcome. We demonstrated that the orders of patients were consistent across repeated experiments, and prognostic genes were successfully nominated. Additionally, it was observed that the risk score can be used to evaluate the molecular aggressiveness of individual patients.

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

--score_threshold (default=0.75)
: When the proportion of chromosomes that exceeds a certain stratification score (*score_threshold*) exceeds a certain proportion in the population (*prop_population*), the survival score is started to be considered in the fitness function.

--prop_population (default=0.95)
: When the proportion of chromosomes that exceeds a certain stratification score (*score_threshold*) exceeds a certain proportion in the population (*prop_population*), the survival score is started to be considered in the fitness function.
  
--seed <int> (default=42)
: Random seed.

## Output Files

1. geneset.txt

| gene  | ascending  |
| --- | --- |
| Gene a  | False  |
| Gene b  | True  |
| Gene c  | False  |
| ...  |   |

The first column (gene): selected genes.

The second column (ascending): If ascending is True, the gene was selected as having relation with worse prognosis when its expression is high. If ascending is False, the gene was selected as having association with worse prognosis when its expression is low.

2. sample_risk.txt

| sample | subtype  | risk  |
| --- | --- | --- |
| Sample A  | Subtype 1  | Low  |
| ...  |   |   |
| Sample B  | Subtype 2  | Intermediate  |
| ...  |   |   |
| Sample C  | Subtype 3  | High  |
| ...  |   |   |

The first column (sample): each sample.

The second column (subtype): subtype of the sample.

The third column (risk): predicted risk within each subtype.

# Citation
Koo, Bonil, et al. "Risk Stratification for Breast Cancer Patient by Simultaneous Learning of Molecular Subtype and Survival Outcome Using Genetic Algorithm-Based Gene Set Selection." Cancers 14.17 (2022): 4120. (https://doi.org/10.3390/cancers14174120)
```

```
