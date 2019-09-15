MetaComp --- A flexible pipeline for comparing bins'(from different metagenomic data ) quality.
==============
<b> This is not a production-ready software repository and is still under active development. Bugs and feature requests will not addressed. Please use with caution.</b>

### E-mail: liaoherui@mail.dlut.edu.cn
### Version: V1.0

--------------

### Abstract
MetaComp evaluates bins(draft genome) from different metagenomic data.These metagenomic data could come from different sequencing platforms or handled by different analysis pipelines.Anyway,if you want to know which platform or assembly tool or even their combinations can bring you the best binning result,just input the fastq files(after quality control) and contig/scaffold(.fa/.fasta) files, MetaComp will help you binning and compare the quality of bins from these data automatically. <BR/>

### Dependencies
* Python >=2.6
* MaxBin 2.2.4
* BWA 0.7.17
* Quast v5.0.0
* Checkm v1.0.12
* Kraken 0.10.6
* Aragorn v1.2.38
* Barrnap 0.9
* Make sure these programs are located on your PATH

### Install(Only for linux)
Install required python libraries:  `pip install  pandas numpy rpy2 matplotlib seaborn pyecharts==0.5.11`

clone software: `git clone https://github.com/liaoherui/MetaComp.git`<BR/>
 

### Usage

**python Lazy_MAEP.py -l DIR -s DIR -o DIR** <BR/>
  
optional arguments:  

**-l** : <BR/>
This option refers to the directory of input list.The list **(tab seperated)** is composed of 4 parts.<BR/>
**column 1: sample name**<BR/>
**column 2: prefix (Usuallly refer to different sequencing platforms or assembly tools' name)**<BR/>
**column 3: assembly result( the directory of  contig/scaffold .fa/.fasta file ->) dir**<BR/>
For PE:<BR/>
**column 4 and column 5: the directory of PE fastq files**<BR/>
(Or,for SE/Long Reads:<BR/>
**column 4 : the directory of SE fastq files**)<BR/>

**Example List:(One sample ,two sequencing platforms)** <BR/>
 ```
 list/contig_raw_reads.list
 ```
  **-s** : <BR/>
 This option refers to the directory of sample name list.For example,if your input data only refers to one sample(suppose the sample name is 'sample_A'), then your sample list should look like:<BR/>
 ```
 sample_A
 ```
 Or, you have two samples(sample_A and sample_B),then,your list should look like:<BR/>
 ```
 sample_A
 sample_B
 ```
Then, for multiple samples ,the form is similar like:<BR/>
 ```
 sample_A
 sample_B
 ...
 sample_N
 ```
 
  **-o** : <BR/>
  This option refers to output directory.
 
### Example
Evaluate 4 different metagenomic data from 4 sequencing platforms.(These data come from the same sample named 'hlj'.) <BR/>
```sh lazy_test.sh```

### Output
There are two parts of output.One is report(.html),another is literature figure.<BR/>
The results can be downloaded in  Report folder(.html) or Figure  folder(literature figure) 



