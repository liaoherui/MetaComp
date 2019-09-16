MetaComp --- A flexible pipeline for evaluating metagenome assemblies by integrating different algorithms.
==============

### E-mail: heruiliao2-c@my.cityu.edu.hk
### Version: V1.0

--------------
### Dependencies
* Python >=2.6
* MaxBin 2.2.4
* BWA 0.7.17
* Quast v5.0.0
* Checkm v1.0.12
* Kraken 0.10.6
* Aragorn v1.2.38
* Barrnap 0.9
Make sure these programs have been added in path

### Install(Only for linux)

####
Install required python libraries:  `pip install  pandas numpy rpy2 matplotlib seaborn pyecharts==0.5.11`

####
`git clone https://github.com/liaoherui/MetaComp.git`<BR/>
 

### Usage

**python Lazy_MAEP.py -l data.list -s samplename -o out_path** <BR/>
  
optional arguments:  

**-l** : <BR/>
We need a list file including five columns with **(tab as delimiter)**.<BR/>
**column 1: Sample id**<BR/>
**column 2: Assembly id**<BR/>
**column 3: Assembly sequences (.fasta)**<BR/>
**column 4: Forwared reads (.fastq for short-reads sequencing)/Long reads (.fastq)
**column 5 (Option): Reverse reads (.fastq for short-reads sequencing)
**This example provides four assemblies for the sample of "hlj" ** <BR/>
 ```
 list/contig_raw_reads.list
 ```
****<BR/>
**-s** : <BR/>

A file contains sample names, each line per sample.

**-o** : <BR/>
Output path.
 
### Output
MetaComp genetats both html report and standalone png figures <BR/>
* 1.Bin Quality Bar
![Bin Quality Bar](https://github.com/liaoherui/MetaComp/blob/master/Figure/1.1_bar.png)
* 2.Bin Quality Stack Bar
![Bin Quality Stack Bar](https://github.com/liaoherui/MetaComp/blob/master/Figure/1.2_stack_bar.png)
* 3.Bin Completeness/Contamination Scatter
![Bin CC Scatter](https://github.com/liaoherui/MetaComp/blob/master/Figure/3.scatter.png)
* 4.Taxomony UpSet
![Bin Taxomony UpSet](https://github.com/liaoherui/MetaComp/blob/master/Figure/4_upset_species.png)
* 5.Bin N50 and Coverage Boxplot
![Bin N50 and Coverage Boxplot](https://github.com/liaoherui/MetaComp/blob/master/Figure/5_boxplot_cov_n50_overall.png)
* 6.Taxonomy Relative Abundance
![Taxonomy Relative Abundance](https://github.com/liaoherui/MetaComp/blob/master/Figure/6_2_stack_bar_species.png)


