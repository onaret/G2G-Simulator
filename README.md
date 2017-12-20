# Introduction
This tool has been developed to simulate*G2G* analysis.  

G2G or genome-to-genome analysis is a joint analysis of host and pathogen genomes that study side by side correlation between host and pathogen systematic variation.  

# Load tool

```R
source("G2G_simulator.R")
```

# Simulate GWAS Case-Control Study

#### Define population
Here we want two sub-populations P1 and P2.  
* From P1, 200 individuals are in case group and 400 in control group.  
* From P2, 400 individuals are in case group and 200 in control group.  

```R
my_population = generate_population_for_GWAS(list(
	`P1` = c(`case` = 200, `control` = 400), 
	`P2` = c(`case` = 400, `control`  = 200)))
```

#### Define genotyping data & run analysis
Here we want 100,040 SNPs, neutral and causal, stratified or not
* Number of neutral SNP is 100,000  
* On this 5% will be stratified  
* 20 non stratified causal SNP will be added with R coefficient between 1 and 2  
* 20 stratified causal SNP will be added with R coefficient between 1 and 2   
* Fixation coefficient for making stratification strength is 0.2  

```R
GWAS_result = GWAS_scenario(populations = my_population, 
 neutral = 100000, 
 neutral_S_rate = 0.05, 
 causal_NS = seq(1,2, by = 0.05), 
 causal_S = seq(1,2, by = 0.05), 
 fst_strat = 0.2)
```

#### Plot results
Plot the results with 3 different conditions :
* Without correction
* With human groups
* With 5 first PCs

On Manhattan plots

```R
plot_GWAS_manhattan(GWAS_result)
```

![plot_GWAS_manhattan()](doc/GWAS_m.png "Manhattan plot of GWAS without correction")

On QQ plots

```R
plot_GWAS_QQ(GWAS_result)
```

![plot_GWAS_QQ()](doc/GWAS_qq.png "QQ-plot of GWAS with computed principal component")

# Simulate G2G setup

#### Define host's and pathogen's sub-populations distribution
Here host sub-population P1 and P2 have the same amount of samples (2500 each) but in P1 1500 samples have strain A and 1000 have strain B and conversely in P2 1000 samples have strain A and 1500 have strain .  

```R
study_design = to_study_design(list(
	`P1` = c(`A` = 1500, `B` = 1000), 
	`P2` = c(`A` = 1000, `B`  = 1500)))
```

#### Setup the the G2G pattern

**rep** is the number of repetition you want to execute to draw the pvalue distribution.  

**s_stratified** is a vector defining the stratification direction you want to set for **S**NP between **human sub-populations**. The order of sub-population defined in the vector give the stratification strength direction. Similar properties apply for the 3 parameters set forth bellow.   
> eg : here there will be more alternate allele in population P1 than population P2.

**s_biased** is the stratification you want to set for **S**NP between **pathogen strains**.  
> eg : Here there will be more alternate allele in sample that has strain A than strain B.

Similarly with viral side...  
**a_stratified** is the stratification you want to set for p**a**thogen variant between **pathogen strains**.  

**a_biased** is the stratification you want to set for p**a**thogen variant between **human sub-populations**.  

**associated_strains** in case of association it define the pathogen strain associated.  

**associated_populations** in case of association it define the the human sub-population associated.  

```R
G2G_setup = get_G2G_setup(rep = 1000, 
  s_stratified = c("P1","P2"), 
  s_biased = c("A","B"), 
  s_partial_bias = c("P1"), 
  a_stratified = c("A","B"))
```

#### Run G2G Setup
**fst_pop_strat** is the fixation coefficient that define strength for SNPs stratification between human sub-population stratification.   

**fst_pop_bias** is the fixation coefficient that define strength for SNPs stratification between pathogen strains.   

**fst_strain_strat** is the fixation coefficient that define strength for pathogen variant stratification between pathogen strains.  

**fst_strain_bias** is the fixation coefficient that define strength for pathogen variant stratification between human sub-population.  

**beta** if there is an association set the log of odd ratio.  

**tag** folder name to save results.  

```R
test_G2G_setup(study_design, G2G_setup, 
  fst_pop_strat = 0.2, 
  fst_pop_bias = 0.2, 
  fst_strain_strat = 0.2,
  tag = 'demo')
```

The results are automatically plotted in the **tag** folder.  

![test_G2G_setup()](doc/G2G_setup.png "G2G setup")

# Simulate G2G study

#### Define you data

Similarly to G2G setup simulation, define study design in term of population and strain distribution.   

```R
study_design =  to_study_design(list(
  `P1` = c(`A` = 250, `B` = 250), 
  `P2` = c(`A` = 250, `B`  = 250)))
```

Here we describe step-wise the data we want to generate as a composition of functions.  

```R
G2G_sample_data =	parse_G2G_config(
	study_design,
	G2G_conf(
		association(
			SNP(1),
			AA(1,
				associated_strains = "full",
				associated_populations = "full",
				stratified = "full",
				fst_strat = 0.2,
				beta = 0.5,
				bio_tag = 'tag'),
			replicate = 100),
		SNP(100, stratified = "full", fst_strat = 0.2),
		SNP(800)))
```

**parse_G2G_config**(**study_design**, **G2G_conf()**) the first argument has to be the study design element, then it will be a undefined number of G2G_conf() function calls (one is standard).  

**G2G_conf**(**...**, replicate = 1) combine **SNP()**, **AA()**, **association()** function calls, and a number of time to replicate the pattern.  

> eg here the described pattern will be repeated 5 times

**SNP**(**size**, stratified = NA, fst_strat=NA, biased = NA, fst_bias=NA, bio_tag=NA)  
* **size** : the number of SNPs
* **stratified**, similar to *s_stratified*
* **biased**, similar to *s_biased*
* **fst_strat**, similar to *fst_pop_strat*
* **fst_bias**, similar to *fst_pop_bias*

**AA**(**size**, stratified = NA, fst_strat=NA, biased = NA, fst_bias=NA, associated_strains = NA, associated_populations =NA, beta=NA, bio_tag=NA)  
* **size** : the number of pathogen variants  
* **stratified**, similar to *a_stratified*
* **biased**, similar to *a_biased*
* **fst_strat**, similar to *fst_strain_strat*
* **fst_bias**, similar to *fst_strain_bias*
* **associated_strains**, similar to  *associated_strains*  
* **associated_populations**, similar to *associated_strains*   
* **beta** in case of association if inside the association() function (see bellow)

**association(...)** combine **SNP()**, **AA()** that will be associated

**bio_tag** a marker that will be associated with the repetitions to make the different groups.

#### Run the G2G analysis
*NOTE:only logistic regression is still maintained*  

```R
analyse_G2G(G2G_sample_data,
  correction(WO_correction = T, W_human_PC = T, W_strain_group = T, W_strain_groups_human_PC = T), 
  analyse(logistic = T), 
  nb_cpu = 40)
```
**analyse** is a vector containing the different analytic methods to run on the dataset
* **logistic** : will run logistic regression
* **gt** : will collapse on human side and run global test
* **skat-L** : will collapse on human side and run skat
* **skato-L** : will collapse on human side and run skato
* **G2** : will collapse on both sides and run G2

**correction** is a vector containing the different correction methods to apply
* **WO_correction** : no correction
* **W_strain_group**: with strain groups
* **W_human_PC**: with 5 first PCs from SNPs data (imputed human groups)
* **W_strain_groups_human_PC**: with 5 first PCs from SNPs data and strain groups

**nb_cpu** : number of available CPU to use

See the script in paper/parse_paper_dataser.R to plot the result