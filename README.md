# ARGON

**Current version:** 0.1.160415 (April 15, 2016).<br>
**Reference:** P.F. Palamara, “ARGON: fast, whole-genome simulation of the discrete time Wright-Fisher process”, [Bioinformatics 2016](http://bioinformatics.oxfordjournals.org/content/early/2016/06/15/bioinformatics.btw355.abstract?keytype=ref&ijkey=Gvb31GJrsCKzeyD). [Preprint (with supplements), doi: http://dx.doi.org/10.1101/036376](http://biorxiv.org/content/early/2016/01/12/036376).

ARGON simulates the discrete time Wright Fisher process (DTWF) backwards in time. The coalescent is equivalent to the DTWF process if the sample size is small compared to the effective population size, but will deviate from it as the sample size increases ([Wakeley and Takahashi, MBE 2003](http://www.ncbi.nlm.nih.gov/pubmed/12598687); [Bhaskar, Clark and Song, PNAS 2014](http://www.ncbi.nlm.nih.gov/pubmed/24469801)). ARGON supports arbitrary demographic history, migration, variable mutation/recombination rates and gene conversion, and efficiently outputs pairwise identical-by-descent (IBD) sharing data.

## Usage

To run ARGON, you can type

		java -jar ARGON.0.1.jar [options]

If no option is specified, the program will generate output using default parameters (10 Mb region for 1,000 samples from a population of 1,000 individuals, with rec=1E-8 and mut=1.65E-8). Command line options are described below.

The input format of ARGON is similar to that of the [GENOME](http://csg.sph.umich.edu/liang/genome/) simulator ([Liang et al., Bioinformatics 2007](http://bioinformatics.oxfordjournals.org/content/23/12/1565.full.pdf)).

### Overview of commands

###### demographic model, sample size, chromosome length
	-N	        Population size or file.
			    (Default = 1000; example: "-N 10000" or “-N model.txt”)
	-pop		Number of samples.
		    	(Default = 1 pop, 1000 samples; example: "-pop 1 1000" or "-pop 2 1000 2000")
	-size		Chromosome length (Mb).
			    (Default = 10 cM, 10 Mb; example: "-size 10")
###### mutation, recombination, and gene conversion rates
	-rec		Recombination rate per base pair.
		    	(Default = 1.0E-8; example: "-rec 1E-8")
	-mut		Mutation rate per base pair.
		    	(Default = 1.65E-8; example: "-mut 1.65E-8")
	-map		Recombination/mutation map file.
	    		(Default = no map; example: "-map map.txt")
	-GC	  	    Non-crossover gene conversion (NCOGC).
		    	(Default = non-crossover GC disabled; example: "-GC 1.0")
	-meanGC	    Mean NCOGC tract length.
		    	(Default = 300; example: "-meanGC 300")
	-minGC	    Minimum NCOGC tract length.
		    	(Default = 0; example: "-minGC 0")
###### output options
	-out		Output files prefix.
				(Default = output to files ARGON.*; example: "-out ARGON")
	-screen		Activate to write output to screen.
				(Default = output to files; example: "-screen")
	-seq		Minimum allele frequency for sequence output.
				(Default = print sequence, MAF = 0.0; example: "-seq 0.0")
	-seq-out	Activate seq output format
				(Default = use VCF file format; example: "-seq-out")
	-haps-out	Activate haps/samples output format
				(Default = use VCF file format; example: "-haps-out")
	-IBD		Output IBD longer than specified cM threshold.
				(Default = do not print IBD; example: "-IBD 1.0")
	-shrink		Output sequence of mutations with small DAF as list (only for seq format).
				(Default = do not shrink sequence; example: "-shrink")
	-quiet		Suppress progress details.
				(Default = print progress; example: "-quiet")
	-age		Output age of alleles.
				(Default = do not output age of alleles; example: "-age")
	-gz			Compress output.
				(Default = do not compress output; example: "-gz")
	-trees		Output Newick trees.
				(Default = do not output trees; example: "-trees")
###### approximation
	-len		Minimum block size.
			    (Default = 1; example: "-len 1" or "-len 10000")
###### other options
	-seed		Random seed.
	    		(Default = random; example: "-seed 1234")
	-help		Print parameter defaults and examples.
		        (Default = do not print parameter defaults and examples; example: "-help")

### Description of command line options

**Demographic model**

-N&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Default = 1000; example: "-N 10000" or “-N model.txt”)

A demographic model is specified using the “-N” flag, followed by one argument. If the argument is an integer number, a single population of constant corresponding haploid size is assumed. If the argument is a file name, the demographic model will be read from the specified file. The format of the demographic model in the file is similar to that of the [GENOME](http://csg.sph.umich.edu/liang/genome/) simulator. *Odd* lines specify haploid population sizes for each generation, starting from generation 0 (present) and proceeding backwards in time. All generations need to be specified from discrete time 0 to the time at which no more population size or migration rate changes occur. The last population sizes specified in the file are assumed to remain constant from the last generation on. *Even* lines specify population merge-split events, and migration rates across populations (backwards in time). Each entry in even lines has format “pop1-pop2_rate”, meaning that individuals will flow from population 1 to population 2 at “rate” individuals per generation. Populations are indexed as they appear in odd lines, using integers starting from 1. If only “pop1-pop2” is specified, “_1” is implied, i.e. all individuals will flow from population 1 to 2 (this may be used to specify population merge/split events). Migration rates are backwards in time. For instance, 1-2_0.1 means that individuals in population 1 have a probability 0.2 to sample a parent in population 2 at the previous generation.

Example:

    0	10000	20000
    1-1_0.998	1-2_0.002	2-2_0.999	2-1_0.001
    1	1000	1000
    1-1	2-1
    2	2000

This file specifies that two populations exist at generation 0. Population 1 has size 10,000, and population 2 has size 20,000. An individual from population 1 at generation 0 has a chance of 0.998 of having an ancestor in population 1, and a chance of 0.002 of having an ancestor in population 2 (this specifies migration from population 2 to 1 forward in time). Note that if the sum of outgoing rates for a population does not sum to 1, it will be automatically normalized so that it does (a warning will be printed). Similarly, an individual from population 2 has a chance 0.999 of having a parent in population 2, and 0.001 of having a parent in population 1. At generation 1, both populations have size 1,000. All individuals at generation 1 will have parents from population 1, which at generation 2 will have size 2,000 individuals. The fourth line of the file is equivalent to “1-1_1	2-1_1”, but “_1” can be omitted.

A few examples can be found in the [FILES/DEMOGRAPHIC_MODELS](https://github.com/pierpal/ARGON/tree/master/FILES/DEMOGRAPHIC_MODELS) folder, including the demographic model used in [Gravel et al. PNAS 2011](http://www.pnas.org/content/108/29/11983.abstract). Additional examples of demographic models can be found in the [GENOME simulator manual](http://csg.sph.umich.edu/liang/genome/). Note, however, that ARGON 0.1 requires all generations to be specified (i.e. there can be no gap between generations), and that the GENOME format does not allow migration rates (introduced with the “_” character in the ARGON format as described above). The demographic model parsed by ARGON can be verified when printed in output before simulation begins.

**Sample sizes**

-pop&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Default = 1 pop, 1000 samples; example: "-pop 1 1000" or "-pop 2 1000 2000")

The number of samples from each population is specified using the “-pop” flag, followed by the number of populations at generation 0, and the number of samples to be obtained from each population. For instance "-pop 1 1000" indicates there will be 1 population at present time, and that 1,000 samples are to be obtained from population 1. "-pop 2 1000 2000” indicates there are 2 populations at present time, 1,000 samples are to be obtained from population 1, and 2,000 samples are to be obtained from population 2. All sizes are haploid.

**Block size approximation**

-len&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Default = 1; example: "-len 1" or "-len 10000")

ARGON can run in approximate mode, which results in improved speed and memory usage. The approximation is similar to that adopted by the GENOME simulator. In approximate mode, the location of recombination events is rounded so that all events occur at the edges of discrete blocks of fixed genetic length. The block length is specified using the “-len” flag. For instance, “-len 10000” requires recombination to occur at the edges of blocks of length 10 micromorgans (1 micromorgan = 0.000001 Morgans). The approximation results in increased correlation for nearby markers (benchmarks for “-len 10000” and “-len 50000” are provided in the supplementary materials of the ARGON paper). 

**Chromosome size**

-size&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Default = 10 cM, 10 Mb; example: "-size 10")

The size (in Mb) of the simulated region is specified using the “-size” flag. For instance, “-size 10” indicates a region of 10 Mb is to be simulated.

**Recombination rate**

-rec&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Default = 1.0E-8; example: "-rec 1E-8")

The recombination rate (per bp, per generation) is specified using the “-rec” flag. For instance, "-rec 1E-8" indicates recombination will happen with probability 1E-8 between two adjacent base pairs (i.e. 1 cM per Mb). The recombination rate can be variable along the genome, in which case the “-map” flag should be used (see below).

**Mutation rate**

-mut&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Default = 1.65E-8; example: "-mut 1.65E-8")

The mutation rate (per bp, per generation) is specified using the “-mut” flag. For instance, "-mut 1E-8" indicates a base pair has a probability of 1E-8 to mutate during one generation. The mutation rate can be variable along the genome, in which case the “-map” flag should be used (see below).

**Variable recombination and/or mutation rate**

-map&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Default = no map; example: "-map map.txt")

The user can specify a variable rate for the recombination and/or mutation process along the genome. To do so, a file should be input using the “-map” flag (e.g. “-map map.txt”). Each line of the file should have format

    Position  recRate  cumRecRate  mutRate

Where “Position” is specified in base pairs, and determines the beginning of a physical interval in which the specified recombination and mutation rates will hold. “recRate” specifies the recombination rate, expressed in centimorgans per megabase (e.g. "-rec 1E-8” is equivalent to 1 cM/Mb). The “cumRecRate” field indicates the genetic position corresponding to the previously specified "Position" field (i.e. the cumulative recombination rate), and “mutRate” specifies the mutation rate for the interval, expressed in chance of mutation per bp per generation (same scale as the “-mut” flag). The [FILES/GENETIC_MAP/](https://github.com/pierpal/ARGON/tree/master/FILES/GENETIC_MAP) folder contains maps from the 1KG project. Note the format is similar to commonly available genetic maps, except the mutation rate field needs to be added. Also see maps available [here](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#gmap).

**Random seed**

-seed&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Default = random; example: "-seed 1234")

A random seed is specified using the “-seed” flag. If not specified, a seed is randomly sampled and output before the simulation begins.

**Sequence output**

-seq&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Default = print sequence, MAF = 0.0; example: "-seq 0.0")

The minimum derived allele frequency of the output sequence is specified using the “-seq” flag, followed by the minimum frequency threshold. Note that this is a derived allele frequency and not a minimum allele frequency.

**Output location**

-out&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Default = output to files ARGON.*; example: "-out ARGON")<br />
-screen&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Default = Not active).

The user can specify the output files prefix using the “-out” flag (e.g. “-out file” will output the sequence data in “file.mut”, see details below). All output can be redirected to screen using the “-screen” flag, in which case files will not be written.

**identical-by-descent (IBD) segments output**

-IBD&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Default = do not print IBD; example: "-IBD 1.0")

IBD segments may be output using the “-IBD” flag, which takes the minimum length (in cM) of output segments as argument. For instance, “-IBD 1.0” outputs all segments longer than 1 centimorgan. Note that small values of this threshold may result in very large output.

**Non-crossover gene conversion**

-GC&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Default = non-crossover GC disabled; example: "-GC 1.0")<br />
-meanGC&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Default = 300; example: "-meanGC 300")<br />
-minGC&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Default = 0; example: "-minGC 0")

Non-crossover gene conversion is activated using the “-GC” flag. The argument of the "-GC" flag indicates the relative rate of non-crossover events. Given a recombination event has been sampled (with rate given in "-rec"), it will be a crossover with probability 1/(1+GC), and non-crossover with probability GC/(1+GC). For instance, If you set "-GC 1", non-crossover events will be as likely as crossover events. Example: "-rec 1E-8 -GC 1.0" --> recombination events are sampled with rate 1E-8, and when a recombination happens, it is crossover with probability 1/(1+1)=1/2, or non-crossover with probability 1/(1+1)=1/2. If "-rec 1E-8 -GC 6.0" --> when a recombination happens, it is crossover with probability 1/(1+6)=1/7, or non-crossover with probability 6/(1+6)=6/7.

The flags “-meanGC” and “minGC” can be used to specify the mean and minimum length of non-crossover tracts (e.g. “-meanGC 300 –minGC 100” will result in non-crossover tracts with average length 300 bp, and all tracts will have length at least 100 bp). Note: a non-uniform genetic map input using the "-map" flag is used to determine the probability that a region hosts a recombination event. Given a recombination event occurs, the specific subtype (crossover or non-crossover) is sampled based on the ratio provided with the “-GC” flag. Therefore, the genetic map does not specify the intensity of crossover recombination events alone.

**Shrinking sequence output**

-shrink&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Default = do not shrink sequence; example: "-shrink")

The "-shrink" flag causes mutations with small derived allele frequency to be output using a different format. This results in smaller sequencing data files. It can only be used with the "-seq-out" flag. See below for format details.

**Quiet output**

-quiet&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Default = print progress; example: "-quiet")

The “-quiet” flag can be activated to suppress printing of simulation progress.

**Compressing the output**

-gz&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Default = do not compress output; example: "-gz")

If the “-gz” flag is activated, most output files will be compressed using gzip.

**Printing command line examples**

-help&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Default = do not print parameter defaults and examples; example: "-help")

If the “-help” flag is activated, default values and a few examples of the available command line arguments are printed.

###Output format

**VCF output**

*NOTE: variants in Haps/samples output are currently not sorted by physical position. Will add an option to have sorted output. If you need sorted output, [Plink](https://www.cog-genomics.org/plink2) should be able to parse an unsorted VCF and output a sorted VCF.* <br>
ARGON outputs VCF files by default. Individual IDs in the VCF have format "popID_indID", where the ID of populations and individuals are integer numbers starting from 1. Diploid individuals are created by merging two haploid individuals, and phase information is maintained (using the "|" separator for the alleles). The full specification of the VCF 4.2 format is available [here](http://samtools.github.io/hts-specs/VCFv4.2.pdf). The VCF format is accepted in input by [Plink 2](https://www.cog-genomics.org/plink2).

**Haps/samples output**

*NOTE: variants in Haps/samples output are currently not sorted by physical position. Will add an option to have sorted output. If you need sorted output, [Plink](https://www.cog-genomics.org/plink2) should be able to parse an unsorted VCF and output a sorted VCF.* <br>
If the -haps-out flag is used the output will be written in haps/samples format (see also [this](http://www.shapeit.fr/pages/m02_formats/hapssample.html) format reference). Sequence data will be written to a “.hap” file (or to screen, with the “MUT” string prepended to each relevant line). Each line will have format

     SNPID1 SNPID2 POS 1 2 ID1_0 ID1_1 ID2_0 ID2_1 ...

where SNPID1 and SNPID2 are identifiers for the variant; POS is the physical position; 1 2 represents the allele values (1=ancestral, 2=derived); and ID1_0 ID1_1 ID2_0 ID2_1 is a list of 1s and 2s corresponding to each (haploid/phased) individual chromosome (0=carries ancestral allele, 1=carries derived allele). The samples file has format

     ID_1 ID_2 missing
     0 0 0
     FAM_ID_1 IND_ID_1 0
     FAM_ID_2 IND_ID_2 0
     ...

where the first two line are header information, and each following line contains family ID, individual ID, fraction missing (always 0 in simulation output). Individuals are diploid (i.e. merge of two haploids). You can convert output from the haps/samples format to vcf using the "--hapsample2vcf" in [BCFTools](https://samtools.github.io/bcftools/bcftools.html). The VCF format is accepted in input by [Plink 2](https://www.cog-genomics.org/plink2).

**Seq output**

If the -seq-out flag is used, Sequence data will be written to a “.mut” file (or to screen, with the “MUT” string prepended to each relevant line). Each line will have format

    posFrom  posTo  numMut  DAF  sequence

The “posFrom” and “posTo” fields specify the physical range where the sampled mutation(s) occurred. This corresponds to the physical range for the ARG branch where the mutation(s) occurred. A specific location for the mutation(s) can be obtained by uniform sampling within the range. The “numMut” field specifies the number of mutations that affected the ARG branch. If this number if greater than 1, it is sufficient to uniformly sample additional mutations within the physical range. The frequency and sequence for all sampled mutations from the same line will be identical. The “DAF” field specifies the derived allele frequency for the mutation(s) (i.e. the number of samples carrying a derived allele). The “sequence” field contains the sequence data. A “0” indicates the corresponding haploid individual does not carry a derived allele, while “1” indicates the individual is a carrier for the derived allele. If the "-shrink" flag was activated, the sequence field for mutations with a small DAF will be output as a list of individuals carrying the derived allele instead of a single string of 0's and 1's. Because only a few individual IDs will be included in the list, this will result in smaller output.

**Newick trees output**

The -trees flag causes marginal [Newick trees](http://evolution.genetics.washington.edu/phylip/newicktree.html) to be output in a ".trees" file. Each line has format

    treeHeight  posFrom  posTo NewickTree

Note that posFrom and posTo are in microMorgans.

**Marginal tree intervals**

The “.map” file contains a list of “breakpoints”. These indicate the starting position of each marginal tree of the sampled ARG. If the output is written to screen, the string “BREAKPOINTS” will introduce this list on stdout.

**Allele age**

If the -age flag is used, allele ages are written in a ".age" file. Each line of the ".age" file contains one allele age, with format

    chr SNP_name SNP_pos branchPhysStart branchPhysEnd mutAge branchGenStart branchGenEnd

where mutAge is the age of the allele, and branchPhysStart, branchPhysEnd, branchGenStart, branchGenEnd are the physical and generation start/end coordinates of the ARG branch that was hit by the mutation. SNP_name SNP_pos are omitted if the seq output is used.

**IBD output**

The IBD output will be written to a “.ibd” file (or to screen, if the “-screen” flag is active, with the “IBD” string prepended to each relevant line). Lines will have format

    ID1  ID2  physFrom  physTo  genFrom  genTo  lengthCM  ancestorAge  ancestorID,

where “ID1” and “ID2” specify the IDs of the individuals sharing an IBD segment, “physFrom” and “physTo” specify the physical range of the IBD segment, “genFrom” and “genTo” specify the range expressed using genetic positions. The “lengthCM” field contains the segment length (in centimorgans). The “ancestorAge” and “ancestorID” fields contain the generation and the ID of the ancestor transmitting the IBD segment.

###Known issues and planned updates

**TODO list**

-	Use physical position for tree from/to output.
-	Enable specifying piecewise constant/exponential periods in demographic file format (fewer lines). May integrate with the [Demographic Language Parser](https://github.com/pierpal/DemographicLanguageParser).

**Known issues**
-	Output is not sorted by phisical position. Will add an option to have sorted output.
-	If the recombination map in input is extremely fine grained, minor rounding issues may occur. This should not be a problem for most analyses.

###Contact
For bug reports or suggestions, please contact me at ppalama AT hsph DOT harvard DOTAGAIN edu

###Version history

- June 13 2016&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Uploading release code and binaries. Added Newick, VCF, Haps/Samples, allele age (160415).
- Jan 13 2016&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Fixed recombination problem. Back to normal speed/memory usage in 160113.
- Jan 10 2016&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ARGON 0.1.160101 released.
