# Bayes Factor Species Delimation with SNP Data

A tutorial on Bayes factor species delimitation with SNP data

## Summary

The term "species delimitation" describes statistical methods that use genetic data to quantify whether certain groups of individuals are best considered at the level of "species" or at lower taxonomic levels (e.g., subspecies, populations). To qualify as species, these groups of individuals are required to show a degree of genetic separation from other groups of individuals, but how that degree is measured differs among the methods developed for this purpose. Perhaps a shortcoming of most of these methods is that they require alignment data or sets of gene trees as input, even though most genomic data is in the form of VCF files. In contrast, BFD\* (Bayes Factor Delimitation with SNPs) ([Leaché et al., 2014](https://doi.org/10.1093/sysbio/syu018)) is a method that is implemented in [SNAPP](http://beast2.org/snapp/) ([Bryant et al., 2012](https://doi.org/10.1093/molbev/mss086)) and can be applied to SNP data. Most tutorials on BFD\* analyses use BEAUti to set up the required BFD\* input files; however, this is inconvenient because BEAUti does not handle VCF files. This tutorial therefore demonstrates how BFD\* analyses can be prepared with a script that accepts VCF input files.

## Table of contents

* [Outline](#outline)
* [Dataset](#dataset)
* [Requirements](#requirements)
* [BFD\* with SNAPP](#snapp)

<a name="outline"></a>
## Outline

In this tutorial I am going to present how the BEAST2 add-on package SNAPP can be used for species delimitation with SNP data. To set up this analysis for a SNP dataset in VCF format, the Ruby script [`bfd_prep.rb`](https://github.com/mmatschiner/bfd_prep) will be used. The files generated in this way will be used as input for SNAPP, and SNAPP will estimate Bayes factor support for two contrasting configurations of species delimitation.


<a name="dataset"></a>
## Dataset

The SNP data used in this tutorial are the filtered dataset used for species-tree inference with SVDQuartets in tutorial [Species-Tree Inference with SNP Data](../species_tree_inference_with_snp_data/README.md). You can find more information about the origin of this dataset in the Dataset section of this other tutorial. In brief, the dataset has been filtered to include only bi-allelic SNPs with a low proportion of missing data, for the 26 individuals from 13 cichlid species listed in the table below. Only SNPs mapping to chromosome 5 of the Nile tilapia genome assembly ([Conte et al. 2017](https://doi.org/10.1186/s12864-017-3723-5)) are included in the dataset, and these have been thinned so that no pair of SNPs is closer to each other than 100 bp.

<center>

| Individual ID | Species ID | Species name                  | Tribe         |
|---------------|------------|-------------------------------|---------------|
| IZA1          | astbur     | *Astatotilapia burtoni*       | Haplochromini |
| IZC5          | astbur     | *Astatotilapia burtoni*       | Haplochromini |
| AUE7          | altfas     | *Altolamprologus fasciatus*   | Lamprologini  |
| AXD5          | altfas     | *Altolamprologus fasciatus*   | Lamprologini  |
| JBD5          | telvit     | *Telmatochromis vittatus*     | Lamprologini  |
| JBD6          | telvit     | *Telmatochromis vittatus*     | Lamprologini  |
| JUH9          | neobri     | *Neolamprologus brichardi*    | Lamprologini  |
| JUI1          | neobri     | *Neolamprologus brichardi*    | Lamprologini  |
| KHA7          | neochi     | *Neolamprologus chitamwebwai* | Lamprologini  |
| KHA9          | neochi     | *Neolamprologus chitamwebwai* | Lamprologini  |
| IVE8          | neocra     | *Neolamprologus crassus*      | Lamprologini  |
| IVF1          | neocra     | *Neolamprologus crassus*      | Lamprologini  |
| JWH1          | neogra     | *Neolamprologus gracilis*     | Lamprologini  |
| JWH2          | neogra     | *Neolamprologus gracilis*     | Lamprologini  |
| JWG8          | neohel     | *Neolamprologus helianthus*   | Lamprologini  |
| JWG9          | neohel     | *Neolamprologus helianthus*   | Lamprologini  |
| JWH3          | neomar     | *Neolamprologus marunguensis* | Lamprologini  |
| JWH4          | neomar     | *Neolamprologus marunguensis* | Lamprologini  |
| JWH5          | neooli     | *Neolamprologus olivaceous*   | Lamprologini  |
| JWH6          | neooli     | *Neolamprologus olivaceous*   | Lamprologini  |
| ISA6          | neopul     | *Neolamprologus pulcher*      | Lamprologini  |
| ISB3          | neopul     | *Neolamprologus pulcher*      | Lamprologini  |
| ISA8          | neosav     | *Neolamprologus savoryi*      | Lamprologini  |
| IYA4          | neosav     | *Neolamprologus savoryi*      | Lamprologini  |
| KFD2          | neowal     | *Neolamprologus walteri*      | Lamprologini  |
| KFD4          | neowal     | *Neolamprologus walteri*      | Lamprologini  |

</center>


<a name="requirements"></a>
## Requirements

* **BEAST2:** BEAST2 is a program for Bayesian phylogenetic analyses, that comes bundled with other tools, such as BEAUti (used to prepare BEAST2 input) and TreeAnnotator (used to process BEAST2 output). For this tutorial, the BEAST2 program package is only required on lynx and will not need to be installed on your local computer.

* **SNAPP:** The [SNAPP method](https://www.beast2.org/snapp/) ([Bryant et al. 2012](https://doi.org/10.1093/molbev/mss086)) add-on package for BEAST2 implements a version of the multi-species coalescent model that mathematically integrates over all possible trees at biallelic loci, and is therefore particularly well suited for SNP data. SNAPP can be installed through the BEAST2 Package Manager; however, since BEAUti will not be used to prepare SNAPP input files, the SNAPP package does not need to be installed on local computers, but only on lynx, as described in the tutorial.

* **MODEL_SELECTION:** The MODEL_SELECTION add-on package for BEAST2 implements two approaches to estimate marginal likelihoods, called Path Sampling and Stepping Stone. This  package does not need to be installed on local computers, but only on lynx, as described in the tutorial.


<a name="snapp"></a>
## BFD\* with SNAPP

As the name indicates, Bayes factor species delimitation is based on a statistic called the [Bayes factor](https://en.wikipedia.org/wiki/Bayes_factor). This Bayes factor compares two different models on the basis of the ["marginal likelihood"](https://en.wikipedia.org/wiki/Marginal_likelihood) that is estimated separately for both of these models, in a Bayesian analysis. While the calculation of the Bayes factor from these marginal likelihoods is easy (namely as their ratio), and there are clear rules for the interpretation of Bayes factors (laid out by [Kass and Raftery, 1995](https://doi.org/10.1080/01621459.1995.10476572); see below), the calculation of the marginal likelihoods themselves is tricky. The marginal likelihood can not simply be calculated from the output of any standard BEAST analysis ([Lartillot and Philippe, 2006](https://doi.org/10.1080/10635150500433722)), but requires specific approaches, with multiple separate analyses that factor in the prior and the posterior to different degrees. Two such approaches, similar to each other, are called "Stepping Stone" and "Path Sampling" ([Xie et al., 2011](https://doi.org/10.1093/sysbio/syq085); [Baele et al., 2012](https://doi.org/10.1093/molbev/msl161)). Both of these can be applied with BEAST2, through an implementation that is - perhaps confusingly - called [Path Sampling](https://www.beast2.org/path-sampling/). Here, we are going to use this implementation to perform the Stepping Stone approach, to estimate the marginal likelihood of two alternative models. In one of these models, two groups of individuals will be considered one and the same species, and in the alternative model, both will be considered separate species. The Bayes factor calculated from the two marginal likelihoods will then quantify the relative support for these two models.

Before we apply BFD\*, we will prepare a subset of the SNP data; this will speed up the computationally challenging analysis. The two most closely related species of our dataset seem to be *Neolamprologus chitamwebwai* ("neochi") and *Neolamprologus walteri* ("neowal"). In tutorial [Divergence-Time Estimation with SNP Data](../divergence_time_estimation_with_snp_data/README.md), the divergence of these two taxa was estimated to have occurred only around 30,000 years ago. It might therefore be questioned whether these taxa should in fact be considered species or not. Both taxa have been described rather recently, by [Verburg and Bills (2007)](https://www.mapress.com/zootaxa/2007f/z01612p044f.pdf). They occur sympatrically at the same location (Bangwe Peninsula in Lake Tanganyika, Tanzania), but are said to differ "in
the linear relations between cheek depth and head length and between body depth and standard length", in the number of scales in the lower lateral line, and in habitat preference - *Neolamprologus chitamwebwai* ("neochi") being associated with "large boulders with sand" while *Neolamprologus walteri* ("neowal") appears to prefer "rubble substrate with fine sediment" ([Verburg and Bills, 2007](https://www.mapress.com/zootaxa/2007f/z01612p044f.pdf)). Besides these two species, we'll also retain *Neolamprologus brichardi* ("neobri"), *Altolamprologus fasciatus* ("altfas"), and *Astatotilapia burtoni* ("astbur") as outgroups, but we'll exclude all other species from the BFD\* analysis. Of each species, we have SNP data for two diploid individuals.

* If you no longer have file `NC_031969.f5.sub4.vcf` (it was produced in tutorial [Species-Tree Inference with SNP Data](../species_tree_inference_with_snp_data/README.md)), download it from GitHub:
	
		wget https://github.com/mmatschiner/lynx/raw/refs/heads/main/divergence_time_estimation_with_snp_data/data/NC_031969.f5.sub4.vcf

* Produce a subset that only includes the ten individuals of the five species *Neolamprologus chitamwebwai* ("neochi"), *Neolamprologus walteri* ("neowal"), *Neolamprologus brichardi* ("neobri"), *Altolamprologus fasciatus* ("altfas"), and *Astatotilapia burtoni* ("astbur"). These individuals are: "KHA7", "KHA9", "KFD2", "KFD4", "JUH9", "JUI1", "AUE7", "AXD5", "IZA1", and "IZC5". To write a new VCF file with data only for these individuals, use the following command:

		bcftools view -s KHA7,KHA9,KFD2,KFD4,JUH9,JUI1,AUE7,AXD5,IZA1,IZC5 -o NC_031969.f5.sub6.vcf NC_031969.f5.sub4.vcf
		
	The reduction of individuals might have led to some sites becoming monomorphic for either the reference or the alternate allele. These could now be excluded these sites again with BCFtools, for example with the `bcftools view` option `-e 'AC==0 || AC==AN || F_MISSING > 0.0'`. However, this is not necessary because the `bfd_prep.rb` script, which will be used to write the XML file for SNAPP, will automatically exclude these sites anyway.
	
* Download the script `bfd_prep.rb` from GitHub:

		wget https://raw.githubusercontent.com/mmatschiner/bfd_prep/refs/heads/main/bfd_prep.rb
		
* To see the options available for generating SNAPP input files with `bfd_prep.rb`, have a look at the help text of the script by using the following commands:

		ruby bfd_prep.rb -h
		
	You'll see that `bfd_prep.rb` accepts input either in Phylip (with option `-p` or `--phylip`) or VCF (option `-v` or `--vcf`) format. In addition, the script requires a table file (option `-t` or `--table`) in which individuals are assigned to species, and a directory needs to be specified (option `-d` or `--directory`). SNAPP will then write a series of subdirectories into this directory, each with another XML file. These XML will differ in the specification of how prior and posterior are factored into the analysis, as required for the Stepping Stone approach. Another important variable is the number of steps to be used for the Stepping Stone approach, corresponding to the number of subdirectories and XML files that will be written. With more steps, the analysis will be computationally more demanding, but the estimate of the marginal likelihood will be more accurate. A good number of steps may be 16. Further parameters to specify are the length of each individual analysis in MCMC iterations (`-l` or `--length`), the percentage that should be removed from each MCMC chain as burnin (`-b` or `--burnin`), the maximum number of SNPs to be included in the analysis (`-m` or `--max-snps`), the minimum distance between these SNPs on the chromosome or scaffold (`-q` or `--min-dist`), and the name of the XML output file that is to be written (`-x` or `--xml`). It is also possible to specify that SNAPP should quit after setting up the series of subdirectories for the Stepping Stone approach (`-c` or `--quit`). This can be useful when the analysis is performed on a remote server, and we will use this option here.
	
* To prepare the table file assigning individuals to species, open a new file named `individuals.txt` with a text editor available on lynx, and write the following text to this new file:

		species	individual
		neochi	KHA7
		neochi	KHA9
		neowal	KFD2
		neowal	KFD4
		neobri	JUH9
		neobri	JUI1
		altfas	AUE7
		altfas	AXD5
		astbur	IZA1
		astbur	IZC5

* With the VCF file `NC_031969.f5.sub5.vcf` and the table file `individuals.txt `, we can now prepare the input file for SNAPP with the script `bfd_prep.rb`. We'll also name the analysis directory "bfd_steps", and we'll specify 16 steps for the Stepping Stone approach, a length of 10,000 iterations for each MCMC, a burnin percentage of 10%, a maximum number of 50,000 SNPs, and a minimum distance between SNPs of 1,000 bp. Importantly, we'll specify that SNAPP should quit after setting up the series of subdirectories for the Stepping Stone approach. The default names for the output XML file and the prefix ("BFD.xml" and "BFD") can be kept. Thus, use the following command to generate the XML input file for SNAPP:

		ruby bfd_prep.rb \
			--vcf NC_031969.f5.sub6.vcf \
			--table individuals.txt \
			--directory bfd_steps \
			--steps 16 \
			--length 10000 \
			--burnin 10 \
			--max-snps 10000 \
			--min-dist 1000 \
			--quit

		
* The script `bfd_prep.rb` should have written a file named `BFD.xml`. You could open that, for example with `less -S BFD.xml` if you'ld like to know more about the settings that we are about to use with SNAPP. If you do so, you'll notice a block near the top of the file that starts with `<run id="pathsampler" ...` and ends with `java -cp ...`. This block specifies how the series of subdirectories should be written.
		
* Next, we need to run SNAPP with this XML file as input to write the series of subdirectories. However, to "run SNAPP", we actually run BEAST2, given that SNAPP is implemented into BEAST2 as an add-on package. Thus we need to make sure that SNAPP is installed as an add-on package. This can be done with BEAST2's utility Package Manager. To see which add-on packages are installed and available for you, use the following command:

		packagemanager -list

* If the SNAPP add-on package should not be installed, use this command to install it for you as a user:

		packagemanager -add SNAPP
		
* In the same way, make sure that the add-on package "MODEL_SELECTION" is installed. This package is needed to perform the Stepping Stone approach. If it is not yet installed and available for you, you can change that with the following command:

		packagemanager -add MODEL_SELECTION

* With both required add-on packages installed, run SNAPP with the XML file as input:

		beast2 BFD.xml

	This should produce a long output on the screen, but nevertheless finish very quickly.
	
* After the above command finishes, check that the analysis directory named `bfd_steps` has been created, and have a look inside of that directory, e.g., with `ls`. You should see that it contains 16 subdirectories named `step0` to `step15`, as well as a single file with ending `.sh`. Each of the 16 step directories should contain two files with ending `.sh` and one XML file named `beast.xml`.

We now need to set up scripts to run SNAPP with each of the 16 files named `beast.xml` as input - one Slurm script and one script to launch the Slurm script for each subdirectory.

* Write a Slurm script named `run_snapp.slurm` with the following content to prepare the SNAPP analysis:

		#!/bin/bash

		# Job name:
		#SBATCH --job-name=snapp
		#
		# Wall clock limit:
		#SBATCH --time=6:00:00
		#
		# Processor and memory usage:
		#SBATCH --ntasks=1
		#SBATCH --cpus-per-task=1
		#SBATCH --mem-per-cpu=1G
		
		# Run snapp.
		beast2 ${xml}

* Write the second script, named `run_snapp.sh`, with the following content:

		# Get the initial directory.
		init_dir=`readlink -f .`
		
		# Change to the step directory.
		cd bfd_steps
		
		# Repeat for each step directory.
		for step_dir in step*
		do
			# Change to the step directory.
			cd ${step_dir}
			
			# Copy the slurm script to the step directory.
			cp ${init_dir}/run_snapp.slurm .
			
			# Set the log files.
			log_out=snapp.out
			log_err=snapp.err
			rm -f ${log_out}
			rm -f ${log_err}
			
			# Run snapp to perform this step of bfd path sampling.
			sbatch -o ${log_out} -e ${log_err} run_snapp.slurm beast.xml
			
			# Return to the bfd directory.
			cd ..
		done

* Launch the second script, which in turn triggers the Slurm script:

		bash run_snapp.sh