# ped_sim - Simulating Genetic Variance onto Family Pedigree

# Overview
`ped_sim` allows the use of simulating genomes onto complex family pedigrees. Using empirical genetic data to initialize the pedigree's 
starting founder genomes, we then will utilize a forward population simulator to create new offspring individuals based crossing between
bi-parental sexual reproduction. In addition to using empirical genetic data for the family's starting founders, 
this software allows you to simulate the founders genomes prior to the main simulation.
This software will the use of complex family pedigrees, with the capacity to handle implicit founders in the case 
of incomplete families.

This command line tool was developed with a python front end and uses SLiM as a forward genetic simulator to simulate genetic variance on non founders.
This software leverages family pedigrees as directed acyclic graph data structures where nodes represent individuals 
and directed edges represent parent-child relationships. While the input must be specified as [networkx](https://networkx.org/documentation/stable/tutorial.html#directed-graphs)
directed graph, you can additionally convert a standard ped file to a networkx based pedigree with this software.
 
**This software is still under development**, please feel free to reach out to me for any questions or suggestions 
I can create to make this code more useable! Miguel.Guardado@ucsf.edu
### Types of simulations preformed
The creation of this software is centered around founder genome initialization, and pedigree conversion. Most family pedigrees 
do not come with the generation of each individual explicitly. The first part of each of the 3 simulations comes with converting the 
pedigree into a slim readable pedigree file, where we know the parents for each descendants, and the generation in which they were conceived.
Currently, this software is able to preform 3 different versions of founder genome initialization. 

`load_founders` - This will be used to simulate an inputted family pedigree based on loading the founder genomes from a user supplied genetic file (vcf format).
This will require the user to input a vcf file that is identical to the number of founder that is inputted inside the family pedigree. In the case 
of implicit founders, it will randomly assign individuals to implicit and explicit founders for the simulation. 
If you want to specify the which founders should correspond to specific id from the vcf file, please refer to load_founder_exact function

`load_founder_exact` - This function is similar to load founder with the exception that we will not randomly assign individuals
from the inputted vcf as founders. We will instead require an additional text file that will specify individual id's from the vcffile and the individual pedigree
id to initialize the founders genome.

`sim_founders` - This function will instead simulate the founder genomes that will be inputted at the beginning of the pedigree simulation. The simulation
will also be done in SLiM, and will require many parameters from the user to define the burn in period of the desired founders genomes. 
Much care is needed for an founder burn in scenerio, for if the founder does not have enough genetic variance in it will 
lead to very similar founders being simulated.

## some other helper functions this software comes with.
`ped_to_networkx` - This is a helper function that is used to convert a standard pedigree (__.ped__) file into a networkx (__.nx__)  based family pedigree.

`networkx_to_ped` - This is a helper function that is used to convert a networkx based pedigree (__.nx__) file into a standard pedigree (__.ped__) file.

`check_founder` - This function is used return the meta information of the networkx pedigree, it will return the 
number of founders, num of implicit founders, and num of descendants inside the family.

### Basic Definitions
Founders - Founders are defined as individuals who do not have any known ascendants. Explicit founders are classified as known
founders who are defined by directed edges in the graph node based pedigree. Implicit founders are defined 
in the case of incomplete pedigrees, when only on one parent is known for a descendants. Since our simulations
need both parents to simulate the offspring, we refer to the missing parent as an implicit founder. My software is able to handle unknown parents.

## ped_sim tutorial
### Install Conda Environment
You will need to have conda installed on your system. To install conda/mini conda onto your system follow this 
[link](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

Once you have conda installed on your system, check its working by the simple command
```bash
conda info
```

This repo includes a `ped_sim.yml` file that includes a virtual conda environment of all the software that is needed. 
Once you clone the repo to you directory, create and load the conda environment as follows.

```
cd ped_sim
conda env create -f ped_sim_env.yml
conda activate ped_sim
python setup.py install

#check pipeline interface is working though run_ped_sim.py
python run_ped_sim.py -h
```
For each of the simulation preformed in here you must specify the type of simulation you want preformed. This will be done
by the `-t` parameter. This parameter is **necessary** for any simulation preformed, or else the software will break.

For this tutorial, we will simulate the genomes of a 3 generation family. This family will have 20 individuals, 
7 of which are the family's founders. This means we will need to initialize the genomes of the 7 founders so we can be 
able to simulate the genetic architecture of the 13 remaining descendants. A visual representation of this family can
be found in test_data/test_fam_ped.png

The data for the family tree will be represented in a standard 
[ped](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format) format. For ped_sim to be able
to read in a user family pedigree file, we will first need to convert the .ped format pedigree into a networkx based 
pedigree. To do this, we will use ped_sim ped_to_networkx functions. After creating the networkx based pedigree, we 
will use the check_founders function to inspect the number of founders and descendants.
```
cd test_data
python ../run_ped_sim.py -t ped_to_networkx -p test_fam.ped -o test_fam

python ../run_ped_sim.py -t check_founders -n test_fam.nx
```
Now that we have a networkx pedigree of the family, we will start with the first type of founder initialization;
simulating the founders genomes before running the pedigree simulations. 

### sim_founders

```bash
python ../run_ped_sim.py -t sim_founders -n test_fam.nx  -o testfam_sim
```
Suppose we have genetic data(.vcf file) of genetic data to initialize the founders with, that is where load_founders 
can be initialized. Reminder that this will assign the individuals in the vcf file to founders randomly. 
### load_founders

```bash
`python ../run_ped_sim.py -t load_founders -n test_fam.nx -v load_founders/lwk_1kg_toydata.vcf -o testfam_load`
```

### load_founders nucleotide explicit simulation
One limitation of this software is that the chromosome and position location will be lost inside the simulations, such that if 
you are simulating chromosome 22, the outputting genetic file will be stored as chr 1 and the exact positions will be lost. 
This might be unideal if you want to compare your simulations with external genetic resources.
`load_founder` and `load_founder_exact` comes with the ability for nucleotide specific simulations.
You will additionally need to provide a fasta sequence of the vcf file you are tyring to reference. 

```bash
python ../run_ped_sim.py -t load_founders -n test_fam.nx -v load_founders/lwk_1kg_toydata.vcf -f load_founders/test_fam_fasta.fa -o testfam_load_wnuc
```
***note that if want family simulations across multiple chromosomes, you should create those simulations independantly/
in parallel, this software is not flexible to distinguish 22 chromosomes inside a single simulations. 

Now we want to be specific about who the founders genomes who get initialize before the simulations. load_founder_exact 
will allow you to provide an additional 2 column text file, for each row there is a founder id and the individual id of the vcf file
### load_founder_exact
```bash
python ../run_ped_sim.py -t load_founders_exact -n test_fam.nx -e load_founders/exact_founder_input.txt -v load_founders/lwk_1kg_toydata.vcf -o testfam_load_exact
```
This function will convert a traditional pedigree (.ped) file into a networkx pedigree (.nx)
### ped_to_networkx
```bash
python ../run_ped_sim.py -t ped_to_networkx -p test_fam.ped -o test_fam
```
This function will convert a networkx pedigree (.nx) to a traditional pedigree (.ped) file
### networkx_to_ped
```bash
python ../run_ped_sim.py -t networkx_to_ped -n test_fam.nx -o test_fam
```
## Output
The user will specify the output prefix of the simulations results via the `-o` user parameter.

`_genomes.vcf` - Genetic file of the family simulation. Input is presented as VCF format, information about vcf_file can
be found [here](https://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/) 

`_slim_pedigree.txt` - This will output the slim readable pedigree output that is read into SLiM for the genetic simulations.
The file looks similar to a ped file, but the columns for each file will be [Generation, Offsrping, Parent1, Parent2].
The file is created internally inside ped_sim, and can be found as a python module inside `scripts/convert_pedigree.py`.

## Full parameters for each simulation

### sim_founder 
This feature will initialize the founders genome by running an additional genetic simulation prior to the main
family simulation. 
##### Required input parameters 

`-t sim_founder`

`-n` - networkx represented family pedigree

`-o` - output prefix of the file

##### Additional parameters  

`-r` - recombination rate to use for founder initialization + simulations (default = 1e-8)

`-mu` - mutation rate to use founder initialization + simulations (default = 1e-8)

`-n_gen` - Number of generations to use for founder genome burn in simulation (default = 12000)

`-n_indiv` - Number of individuals to use for founder genomes burn in simulation (default = 1000)

`-l` - length of the genome to simulation for founder burn in and pedigree simulation.

`-s` - Seed number to use for the genetic simulation (default = 1)

### load_founders

This simulation will initialize the founders genomes randomly from the inputted vcf file. 
#####Required input parameters
`-t load_founders`

`-n` - networkx represented family pedigree

`-v` - vcf file for founder genome initialization

`-o` - output prefix of the file

##### Additional parameters
`-f` - fasta file for the inputted vcf file, this will activate nucleotide specifc simulations

`-r` - recombination rate to use for founder initialization + simulations (default = 1e-8)

`-mu` - mutation rate to use founder initialization + simulations (default = 1e-8)

`-s` - Seed number to use for the genetic simulation (default = 1)


### load_founders_exact

This simulation will initialize the founder genomes by providing an additional file to exactly map each founder to a 
individual inputted from a vcf file. 

##### Required input parameters:

`-t load_founders_exact`

`-n` - networkx represented family pedigree

`-v` - vcf file for founder genome initalization

`-e` - 2 column text file matching the pedigress founder to individuals in vcf file to initalizatine the genomes non randomly.
`[vcf_sample_id, networkx_id]`

`-o` - output prefix of the file


##### Additional parameters
`-f` - fasta file for the inputted vcf file, this will activate nucleotide specifc simulations

`-r` - recombination rate to use for founder initialization + simulations (default = 1e-8)

`-mu` - mutation rate to use founder initialization + simulations (default = 1e-8)

`-s` - Seed number to use for the genetic simulation (default = 1)


### networkx_to_ped
This feature will convert networkx represented pedigree files into traditional family pedigrees files. 

##### Required input parameters:
`-n` - networkx represented family pedigree
`-o` - output prefix of the file

##### Output
{output_prefix}.ped - traditional pedigree file

### ped_to_networkx
This feature will convert traditional pedigree files into networkx represented family pedigrees. 

Required input parameters:

`-p` - Traditional family pedigree file

`-o` - output prefix of the file

##### Files that get outputted
{output_prefix}.nx - traditional networkx file

### check_founders
This feature is used for checking the number of explicit, implicit, and descendants created in the simulation.

##### Required input parameters

`-p` - Traditional family pedigree file

output of this function is outputted on the terminal command line

### convert_pedigree
Internal debugging function that can be used to convert networkx represented family pedigree to a slim readable family pedigree.

Required input parameters:
`-n` - networkx represented family pedigree
`-o` - output prefix of the file

##### Files that get outputted
{output_prefix}_ - traditional networkx file





## References
Haller, B.C., and Messer, P.W. (2019). SLiM 3: Forward genetic simulations beyond the Wright–Fisher
model. Molecular Biology and Evolution 36(3), 632–637. DOI: https://doi.org/10.1093/molbev/
msy228

Aric A. Hagberg, Daniel A. Schult and Pieter J. Swart, “Exploring network structure, dynamics, and function using NetworkX”, 
in Proceedings of the 7th Python in Science Conference (SciPy2008), Gäel Varoquaux, 
Travis Vaught, and Jarrod Millman (Eds), (Pasadena, CA USA), pp. 11–15, Aug 200
