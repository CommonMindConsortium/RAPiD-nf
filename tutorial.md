
# Tutorial on RAPiD Nextflow pipeline for RNA-seq data
Gabriel Hoffman

September 7, 2020

Start by download RAPID to a local directory from [https://github.com/CommonMindConsortium/RAPiD-nf](https://github.com/CommonMindConsortium/RAPiD-nf).

Paths here refer to Mount Sinai's Minerva HPC system, but should be changed for a different system. The reference genome (GRCh38) and transcripome (Gencode.v30) used here are specified by `--genome GRCh38.Gencode.v30`.  Note that the `config/chimera.config` specifies the locations of annotations files on the local system.

#### Set up variables and paths

```bash
# load NextFlow
ml nextflow/19.10.0

# where to store temporary work files
SCRATCH=/sc/arion/scratch/hoffmg01/

WORK=/sc/arion/scratch/hoffmg01/rapid_run
cd $WORK

# Directory for batch
# each subfolder contains the FASTQ's from a single sample:
# The format of the FASTQ should be
# NAME_GTCCGC_L001_R2_001.C4HALACXX.fastq.gz
# NOTE: path *cannot* end with a backslash
BATCH=/sc/arion/scratch/hoffmg01/test_rapid

# Location of RAPiD pipeline
NF=/hpc/users/hoffmg01/build2/RAPiD-nf/RAPiD.nf

# Be sure to specify the correct strand
# Most kits are "reverse", but check with the wet lab team
STRAND='none'
```
If there is a single sample VCF in the path specified by `--vcfPath`, `STAR` will run `WASP` to identify reads overlaping variants that are succesptible to reference bias (i.e. reads with the reference allel are most likely to map.). It will then folder out biased reads in to a `*wasp.bam` file

In these examples, most methods are enabled.  You can remove `--rsem` to same time.


#### Test run
Note that this just checks the files: `--dryRun` means that jobs won't be submitted to queue.

```shell
nextflow run $NF \
	--run $BATCH \
	-work-dir $SCRATCH \
	--genome GRCh38.Gencode.v30 \
	--stranded ${STRAND} \
	-profile chimera \
	--trimAdapter TruSeq3-PE \
	--rawPath . \
	--vcfPath . \
	--wasp \
	--continueWithoutVcf \
	--fastqc \
	--kallisto \
	--leafcutter \
	--featureCounts \
	--qc \
	--salmon \
	--pathogen \
	--dryRun 
	
# --rsem \ Omit RSEM, it can be slow
```

#### Real run
Run this in a `screen` session or submit as a job.

```shell
nextflow run $NF \
	--run $BATCH \
	-work-dir $SCRATCH \
	--genome GRCh38.Gencode.v30 \
	--stranded ${STRAND} \
	-profile chimera \
	--trimAdapter TruSeq3-PE \
	--rawPath . \
	--vcfPath . \
	--wasp \
	--continueWithoutVcf \
	--fastqc \
	--kallisto \
	--leafcutter \
	--featureCounts \
	--qc \
	--salmon \
	--pathogen \
	--txrevise 
```

When jobs are submitted to queue, you can use `bjobs` to look at all your jobs.


#### Resume run
If a run crashes, you can resume using cached files

```shell
nextflow run $NF \
	--run $BATCH \
	-work-dir $SCRATCH \
	-resume \ # flag to resume using cached data
	--genome GRCh38.Gencode.v30 \
	--stranded ${STRAND} \
	-profile chimera \
	--trimAdapter TruSeq3-PE \
	--rawPath . \
	--vcfPath . \
	--wasp \
	--continueWithoutVcf \
	--fastqc \
	--kallisto \
	--leafcutter \
	--featureCounts \
	--qc \
	--salmon \
	--pathogen \
	--txrevise 
```



<!---
setfacl -m u:wangy33:rwx /sc/arion/scratch/hoffmg01/test_rapid/
--->



