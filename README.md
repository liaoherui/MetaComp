MetaComp(Metagenomics Assembly Compare Pipeline)
==============
<b> This is not a production-ready software repository and is still under active development. Bugs and feature requests will not addressed. Please use with caution.</b>

### E-mail: liaoherui@mail.dlut.edu.cn
### Version: V1.0

--------------

### Abstract
MetaComp uses a combination of different features and algorithms to compare the qualities of different metagenomic datasets. The datasets are defined as the de novo assembly results(contigs/scaffolds) from diffenrent types of raw reads which are processed variously.<BR/>



### Install(Only for linux)

Install 3rd party programs:

* [MaxBin 2.2.4](https://downloads.jbei.org/data/microbial_communities/MaxBin/MaxBin.html)
* [BWA 0.7.17](http://bio-bwa.sourceforge.net/)
* [Quast v5.0.0](https://sourceforge.net/projects/quast/files/)
* [Checkm v1.0.12](https://github.com/Ecogenomics/CheckM/wiki/Installation)
* [Kraken 0.10.6](http://ccb.jhu.edu/software/kraken/)
* [Aragorn v1.2.38](https://anaconda.org/bioconda/aragorn)
* [Barrnap 0.9](http://www.vicbioinformatics.com/software.barrnap.shtml)
* Make sure these programs are located on your PATH
* Tested versions indicated, but other versions might also work
* Quick Start:<BR/>
  Run the command:
  `git clone https://github.com/liaoherui/MetaComp.git`<BR/>
 

### Manuals
* Quick Start:<BR/>
  `python Lazy_MAEP.py -l list/contig_raw_reads.list -s list/sample.list -o zxy_p4_test `<BR/>
  
* Option Illustration:<BR/>

**-l** : <BR/>
This option refers to the input list.The list **(tab seperated)** is composed of 4 parts.<BR/>
**column 1: sample name**<BR/>
**column 2: prefix (Usuallly refer to different sequencing platforms or assembly strtegies)**<BR/>
**column 3: assembly result(.fasta file with multi contigs/scaffolds) dir**<BR/>
**column 4 and column 5: PE raw reads dir**<BR/>
**Example List:(One sample ,two sequencing platforms)** <BR/>
 ```
 sample_ID Athena_10X  /mnt/10X_athena.fasta /mnt/osf1/_10X_R1.fq.gz  /mnt/osf1/zxy_10X_R2.fq.gz 
 sample_ID Illumina /mnt/Illumina.fasta  /mnt/osf1/zxy_1.fq.gz /mnt/osf1/zxy_2.fq.gz
 ```
  **-s** : <BR/>
 This option refers to the sample name list.For example,if your input data only refers to one sample(suppose the sample name is 'zxy'), then your sample list should be:<BR/>
 ```
 zxy
 ```
 Or, you have two samples(zxy and hlj),then,your list should be:<BR/>
 ```
 zxy
 hlj
 ```
Then,multiple samples' condition is similar.<BR/><BR/>

  **-o** : <BR/>
  This option refers to output dir.
 
 


### Output
There are two parts of output.One is report(.html),another is literature figure.<BR/>
The results can be downloaded in  Report folder(.html) or Figure  folder(literature figure) 



