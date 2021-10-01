# Introduction
This tool has been developed to simulate G2G analysis.  

G2G or genome-to-genome analysis is a joint analysis of host and pathogen genomes that study side by side correlation between host and pathogen systematic variation.  

#### Load tool

```R
source("G2G_simulator.R")
```

# OPTION 1: Simulate a full G2G study

Similarly to G2G simplified, a study design must be fined in term of host populations (population P1, P2...) and pathogen strains distribution (strain A, B,...).   

## A. Define G2G data structure

#### Description:  
**G2G_conf** defines the G2G data structure through a composition of SNP, AA and association function calls.  

#### Usage:  
G2G_conf(SNP, AA, association, ...)   

```R
G2G_conf =	get_G2G_conf(
		association(
			SNP(1),
			AA(1,
				stratified = c("A","B"),
				fst_strat = 0.2,
				beta = 0.5,
				bio_tag = 'tag'),
			replicate = 100),
		SNP(100, stratified = "full", fst_strat = 0.2),
		SNP(800)))
```

#### Arguments: 
* **SNP**, _fun_: SNP function call
  * **description:** defines (a) SNP(s). SNPs corresponding to the variations of the host side
  * **Usage:** SNP(size, stratified = NA, fst_strat=NA, biased = NA, fst_bias=NA, bio_tag=NA)
  * **Arguments:**
    * **size**, _int_: the number of SNPs
    * **stratified**, _vector of strings_: host populations groups order give the direction of the stratification from higher MAF to lower MAF
    * **fst_strat**, _int_: is the fixation coefficient that defines the stratification magnitude defined by **stratified**
    * **biased**, _vector of strings_: include a bias such that, the pathogen strains are associated with host stratification (regardless of the defined populations). The order gives the direction from higher MAF to lower MAF
    * **fst_bias**, _int_: is the fixation coefficient that defines the stratification magnitude defined by **biased**

* **AA**, _fun_:  AA function call
  * **description:** defines (an) AA(s). AAs for amino acids correspond to the variations on the pathogen side
  * **Usage:** AA(size, stratified = NA, fst_strat=NA, biased = NA, fst_bias=NA, beta=NA, bio_tag=NA)
  * **Arguments:**
    * **size**, _int_: the number of pathogen variant
    * **stratified**, _vector of strings_: pathogen strains order give the direction of the stratification from higher MAF to lower MAF 
    * **fst_strat**, _int_: is the fixation coefficient that defines the stratification magnitude defined by **stratified**
    * **biased**, _vector of strings_:  include a bias such that, the host populations are associated with pathogen stratification (regardless of the defined pathogen strains). The order gives the direction from higher MAF to lower MAF
    * **fst_bias**, _int_: is the fixation coefficient that defines the stratification magnitude defined by **biased**
    * **beta**, _int_: in case of association (and therefore inside the association() function call (see bellow)), the log of odd ratio.

* **association**, _fun_: association function call 
  * **description:** defines an association between (a) SNP(s) and (a) AA(s) 
  * **Usage:** association(SNP, AA, replicate)**
  * **Arguments:**
    * **SNP**, _fun_: is a SNP function call outcome, the number of SNP will define how many are associated with the AA function call
    * **AA**, _fun_: is a AA function call outcome, the number of AA will define how many are associated with the SNP function call
    * **replicate**, _int_: is the number of time such an association is added

* **...**, _fun_: other AA, SNP or association function calls

* **bio_tag**, _string_: a tag that will be added in the generated dataset.

## B. Define the host populations and pathogen strains distributions

#### Description:  
**get_study_design** defines the host populations and pathogen strains distributions

#### Usage:  
get_study_design(structure)  

```R
study_design =  get_study_design(list(
  `P1` = c(`A` = 250, `B` = 250), 
  `P2` = c(`A` = 250, `B`  = 250)))
```

#### Arguments:
**structure**, _list of nammed vector of nammed int_:  defines the study design with the host populations P1 and P2 and their respective proportion in pathogen strains A and B  
> eg : Here we have the same number of samples in each host population (500)  with each 250 with strain A and strain B  

## C. Generate G2G data

#### Description:  
**get_G2G_data** generates the G2G data

#### Usage:  
get_G2G_data(study_design, G2G_conf)  

```R
G2G_data =	get_G2G_data(
	study_design,
	G2G_conf)
```

<!---_ get_G2G_data(study_design, G2G_conf, ...)  
-->

#### Arguments: 

**study_design**, _fun_: get_study_design function call

**G2G_conf**, _fun_: G2G_conf function call

<!---_ **...** More G2G conf function outocmes can be passed. -->

## D. Analyse the G2G data

<!---_ *NOTE:only logistic regression is still maintained*  -->

<!---_ analyse(logistic = T), -->

#### Description:  
**analyse_G2G** runs the G2G analysis

#### Usage:  
analyse_G2G(G2G_data, correction, nb_cpu = 40)   

```R
analyse_G2G(G2G_data,
  get_correction(WO_correction = T, W_host_PC = T, W_pathogen_group = T, W_pathogen_groups_host_PC = T), 
  nb_cpu = 40)
```

#### Arguments: 
<!---_ **analyse** is a vector containing the different analytic methods to run on the dataset 
* **logistic** : will run logistic regression
* **gt** : will collapse on human side and run global test
* **skat-L** : will collapse on human side and run skat
* **skato-L** : will collapse on human side and run skato
* **G2** : will collapse on both sides and run G2
-->

* **G2G_data** is get_G2G_data function call

* **correction**, _fun_:  get_correction function call
  * **description:** defines the series of corrections to assess
  * **Usage:** get_correction(WO_correction = F, W_human_PC = F, W_pathogen_group = F, W_pathogen_groups_host_PC = F)
  * **Arguments:**
    * **WO_correction**, _bool_ : no correction
    * **W_pathogen_group**, _bool_: with pathogen strains
    * **W_host_PC**, _bool_: with 5 first PCs from SNPs data (imputed human groups)
    * **W_pathogen_groups_host_PC**, _bool_: with 5 first PCs from hosts data and pathogen strains

* **nb_cpu**, _int_ : number of available CPU to use

See the script in paper/parse_paper_dataser.R to plot the result


# OPTION 2: Simulate a single genome-to-genome (G2G) association case

## A. Define host's populations and pathogen's strains distribution

```R
study_design = get_study_design(structure = list(
	`P1` = c(`A` = 1500, `B` = 1000), 
	`P2` = c(`A` = 1000, `B`  = 1500)))
```

#### Description:  
**get_study_design** defines the host and pathogen structure

#### Usage:  
get_study_design(structure)  

#### Arguments:
**structure**, _list of nammed vector of nammed int_:  defines the study design with the host populations P1 and P2 and their respective proportion in strains A and B  
> eg : Here we have the same number of samples in each host population (2500) but in P1 1500 samples have strain A and 1000 strain B and conversely in P2, 1000 samples have strain A and 1500 strain B.   


## B. Define the correlation structure

```R
G2G_setup = get_G2G_setup(rep = 1000, 
  s_stratified = c("P1","P2"), 
  s_biased = c("A","B"),
  a_stratified = c("A","B"))
```
#### Description:  
**get_G2G_setup** allows to specify the stratification direction

#### Usage:  
get_G2G_setup(rep, s_stratified = NA, s_biased = NA, a_stratified = NA, a_biased = NA)  

#### Arguments: 

**rep**, _int_: is the number of repetition you want to execute to draw the pvalue distribution.  

**s_stratified**, _vector of strings_: host populations groups order give the direction of the stratification from higher MAF to lower MA
> eg : here there will be a higher minor allele frequencyi (MAF) in population P1 than in population P2.

**s_biased**, _vector of strings_: include a bias such that, the pathogen strains are associated with host stratification (regardless of the defined sub-populations groups). The order gives the direction from higher MAF to lower MAF
> eg : Here there will be a higher MAF for the hosts that have strain A than strain B.
> In conlcusion, the MAF decreases with a maximum for P1 with strain A (P1.A) to P1.B, P2.A and finally P2.B 

Similarly for the variants on the pathogen side...  
**a_stratified**,  _vector of strings_: pathogen strains order give the direction of the stratification from higher MAF to lower MAF

**a_biased**, _vector of strings_: include a bias such that, the host populations groups are associated with pathogen stratification (regardless of the defined pathiogen strains). The order gives the direction from higher MAF to lower MAF

<!---_
**associated_strains** in case of association it define the pathogen strain associated.  

**associated_populations** in case of association it define the the host populations associated.  
-->


## C. Run simplified G2G

```R
test_G2G_setup(study_design, G2G_setup, 
  fst_host_strat = 0.2, 
  fst_host_bias = 0.2, 
  fst_pathogen_strat = 0.2,
  tag = 'demo')
```
#### Description:  
**test_G2G_setup** runs the simplified G2G

#### Usage:  
test_G2G_setup(study_design, G2G_setup, fst_host_strat = NA, fst_host_bias = NA, fst_pathogen_strat = NA, fst_pathogen_bias=NA, tag = 'unnamed')  

#### Arguments: 

**study_design**, _fun_: get_study_design function call

**G2G_setup**, _fun_: get_G2G_setup function function call

**fst_host_strat**, _int_: is the fixation coefficient that defines the stratification magnitude defined by **s_stratified**

**fst_host_bias**, _int_: is the fixation coefficient that defines the stratification magnitude defined by **s_biased**

**fst_pathogen_strat**, _int_: is the fixation coefficient that defines the stratification magnitude defined by **a_stratified**

**fst_pathogen_bias**, _int_: is the fixation coefficient that defines the stratification magnitude defined by **a_biased**

**tag** folder name to save results

The results are automatically plotted in the **tag** folder.

![test_G2G_setup()](doc/G2G_setup.png "G2G setup")



# OPTION 3: Simulate Case-Control GWAS

#### Define populations

```R
my_population = generate_population_for_GWAS(list(
	`P1` = c(`case` = 200, `control` = 400), 
	`P2` = c(`case` = 400, `control`  = 200)))
```
Here we want two sub-populations P1 and P2. 
* From P1, 200 individuals are in case group and 400 in control group.
* From P2, 400 individuals are in case group and 200 in control group.

#### Define genotyping data & run analysis

```R
GWAS_result = GWAS_scenario(populations = my_population, 
 neutral = 100000, 
 neutral_S_rate = 0.05, 
 causal_NS = seq(1,2, by = 0.05), 
 causal_S = seq(1,2, by = 0.05), 
 fst_strat = 0.2)
```
Here we want 100,040 SNPs, neutral and causal, stratified or not
* Number of neutral SNP is 100,000
* On this 5% will be stratified
* 20 non stratified causal SNP will be added with R coefficient between 1 and 2
* 20 stratified causal SNP will be added with R coefficient between 1 and 2
* Fixation coefficient for making stratification strength is 0.2


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

