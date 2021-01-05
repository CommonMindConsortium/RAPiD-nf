# RAPiD-nf

RAPiD is an RNA-Seq pipeline used and maintained by Ying-Chih Wang and Hardik Shah in the [Department of Genetics and Genomic Sciences](https://icahn.mssm.edu/research/genomics), [Icahn Institute for Data Science and Genomic Technology](http://datascience.icahn.mssm.edu).  RAPiD is implemented using the [Nextflow](https://www.nextflow.io) pipeline engine and the [LSF](https://www.ibm.com/products/hpc-workload-management) job manager to enabl processing of raw RNA-seq FASTQ files for analysing using:

 - [trimmomatic v0.36](http://www.usadellab.org/cms/index.php?page=trimmomatic)
 - [picard v2.20.0](https://broadinstitute.github.io/picard/)
 - [fastqc v0.11.8](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
 - [STAR v2.7.2a](https://github.com/alexdobin/STAR)
 - [STAR/WASP](https://github.com/alexdobin/STAR/blob/master/RELEASEnotes.md#star-260a-20180423)
 - [featurecounts v1.6.3](http://subread.sourceforge.net)
 - [salmon v0.13.1](https://combine-lab.github.io/salmon/)
 - [kallisto v0.45.0](https://pachterlab.github.io/kallisto/)
 - [rsem v1.3.1](https://deweylab.github.io/RSEM/)
 - [leafcutter v0.2.8](https://davidaknowles.github.io/leafcutter/)

# Running RAPiD
See [tutorial](https://github.com/CommonMindConsortium/RAPiD-nf/blob/master/tutorial.md) on running RAPiD using Nextflow v19.10.0 and LSF v10.1.0.0.

Authors
-------
 - Ying-Chih Wang
 - Hardik Shah

Collaborators
-------------
 - Kaur Alasoo
 - Gabriel Hoffman
 - Veera Manikandan
 - Panos Roussos
 - Solly Sieberts
 - Laura Sloofman
 - Georgios Voloudakis
 
