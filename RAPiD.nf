/*
===============================================================
 RAPiD pipeline wrapped in Nextflow
===============================================================
*/

def helpMessage() {
    println """
    ===================================
     ${workflow.manifest.name} v${workflow.manifest.version}
    ===================================

    Usage:

    nextflow run RAPiD-nf --run [Run_Folder] --genome [Genome] --stranded [Strandness] -profile [Profile]

    Mandatory arguments:
      --run                         Run directory
      --genome                      Reference genome
                                      Human:GRCh38.Gencode.v30
                                      Mouse:GRCm38.Gencode.vM21
      --stranded                    Strandness of the sample [none|forward|reverse]
      -profile                      Configuration profile to use. [local|chimera|liteQC]
                                      liteQC:STAR+kallisto+featureCounts+QC

    Optional arguments:
      Main pipeline arguments:
        Path:
        --startFrom                   Run pipeline from FASTQ or BAM(assuming BAI exist) [fastq(default)|bam]
        --rawPath                     Relative path to sample dir to search for FASTQ files [Raw/Illumina]
        --outPath                     Relative path to sample dir to store outputs [Processed/RAPiD]
        --bamPath                     Relative path to sample dir to search for BAM files [Processed/RAPiD/bams]
                                        (at most 1 bam file per sample, *.wasp.bam is automatically ignored)
        --singleEnd                   For single-end dataset
        --trimAdapter                 Adaptors to be trimmed
        --noMarkDup                   Disable MarkDuplicates
        --wasp                        Enable WASP mode
          --vcfPath                   Relative path to sample dir to search for VCF files [VCF]
          --continueWithoutVcf        Continue running STAR without WASP mode for samples without VCF files
        --qc                          Enable QC pipeline
        --target                      Provide a target interval list for picard hsmetric
        --bait                        Provide a bait interval list for picard hsmetric [use --target if not provided]
        --fastqc                      Enable FASTQC
        --pathogen                    Use STAR to align unmapped reads to pathogen sequences
        --rsem                        Enable rsem
        --featureCounts               Enable featureCounts
        --enhancer                    Enable featureCounts for enhancer RNA
        --kallisto                    Enable kallisto
        --salmon                      Enable salmon
        --txrevise                    Enable txrevise_salmon [GRCh38.ensembl.v96]

      Sorted-BAM-downstream arguments:
        --leafcutter                  Enable regtools+Leafcutter
        --bigwig                      Generate BigWig with bedGraphToBigWig
        --bamCoverage                 Generate BigWig with bamCoverage
        --globalVCF                   Set a global VCF file
          --verifyBamID               Enable verifyBamID
          --QTLtools_mbv              Enable QTLtools mbv

      Dry run:
        --dryRun                      Show pipeline configurations and detected input files

      Nextflow:
        -w/-work-dir                  Working dir ["work"]
        -resume                       Resuming previous run
    """.stripIndent()
}

// Show help emssage
if (params.help || !params.genome || !params.stranded || !params.run){
    helpMessage()
    exit 0
}

if (params.run.startsWith("/")) {
    run_path = params.run
    run_name = params.run.replaceAll("(.*/)", "")
} else {
    run_path = "$PWD/${params.run}"
    run_name = run_path.replaceAll("(.*/)", "")
}

if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "Genome: ${params.genome} is not available in the configuration file. Available genomes: [${params.genomes.keySet().join(", ")}]"
}

if (params.adapters && params.trimAdapter && !params.adapters.containsKey(params.trimAdapter)) {
    exit 1, "Adapter: ${params.trimAdapter}' is not available in the configuration file. Available adapter: [${params.adapters.keySet().join(", ")}]"
}

// Reference index path loaded from configuration
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.star_index = params.genome ? params.genomes[ params.genome ].star ?: false : false
params.kallisto_index = params.genome ? params.genomes[ params.genome ].kallisto ?: false : false
params.salmon_index = params.genome ? params.genomes[ params.genome ].salmon ?: false : false
params.rsem_index = params.genome ? params.genomes[ params.genome ].rsem ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.rrna = params.genome ? params.genomes[ params.genome ].rrna ?: false : false
params.reflat = params.genome ? params.genomes[ params.genome ].reflat ?: false : false
params.enhancer_bed = params.genome ? params.genomes[ params.genome ].enhancer ?: false : false
params.txrevise_indices = params.containsKey('txrevise_indices')? params.txrevise_indices : false

// set params.bam flag from --startFrom
params.bam = params.startFrom == "bam"

if(params.trimAdapter) {
    params.adapter = params.adapters[ params.trimAdapter ].fasta
}

if (!params.star_index) {
    exit 1, "STAR index not defined."
}

if (!params.gtf) {
    exit 1, "GTF annotation file not defined."
}

if (params.qc && !params.reflat) {
    exit 1, "REFLAT annotation file not defined."
}

if (params.qc && !params.rrna) {
    exit 1, "rRNA interval list not defined."
}

if (params.kallisto && !params.kallisto_index) {
    exit 1, "kallisto index not defined."
}

if (params.salmon && !params.salmon_index) {
    exit 1, "salmon index not defined."
}

Channel.fromPath("${params.txrevise_indices}/*.index", type: 'dir')
.ifEmpty { exit 1, "no txrevise indices found."}
.set { txrevise_salmon_indices }
if (params.txrevise) {
    if (!params.txrevise_indices) {
        exit 1, "txrevise indices folder not defined"
    }
}

if (params.rsem && !params.rsem_index) {
    exit 1, "rsem index not defined."
}

if (params.pathogen && !params.containsKey('pathogen_star_index')) {
    exit 1, "Pathogen index not defined."
}

if (params.enhancer &&  !params.enhancer_bed) {
    exit 1, "enhancer RNA bed file not defined."
}

println "========================================="
println "${workflow.manifest.name} v${workflow.manifest.version}"
println "========================================="
println "Run Folder: $run_path"
println "Genome: ${params.genome}"
println "Stranded: ${params.stranded}"
println "PairedEnd: ${!params.singleEnd}"
println "Adapter: ${params.trimAdapter}"
println "WASP: ${params.wasp} ${params.continueWithoutVcf?'- continue without vcf':''}"
if (params.target){ println "Target: ${params.target}"}
println "============  Main Pipeline  ============"
if (!params.bam) {
    print "Trimmomatic[${params.trimAdapter?:'-'}]\n"
    print "STAR[*]${params.wasp?'[WASP]':''}\n"
    print "FASTQC[${params.fastqc?'+':'-'}]\n" 
    print "Filter[${params.wasp?'WASP':'-'}]/MarkDup[${params.noMarkDup?'-':'+'}]\n"
    print "QC[${params.qc?'+':'-'}]\n"
    print "Pathogen[${params.pathogen?'+':'-'}]\n"
    print "featureCounts[${params.featureCounts?'+':'-'}]${params.enhancer?'[enhancer]':''} "
    print "RSEM[${params.rsem?'+':'-'}] "
    print "Kallisto[${params.kallisto?'+':'-'}] "
    print "Salmon${params.salmon?'[+]':'[-]'}${params.txrevise?'[txrevise]':''}\n"
} else {
    print "  Omitted - pipeline starts from BAM\n"
}
println "============  BAM downstream  ==========="
print "LeafCutter[${params.leafcutter?'+':'-'}]\n"
print "BigWig[${params.bigwig?'bedGraphToBigWig':'-'}][${params.bamCoverage?'bamCoverage':'-'}]\n"
print "globalVCF:${params.globalVCF} [${params.verifyBamID?'verifyBamID':'-'}][${params.QTLtools_mbv?'QTLtools_mbv':'-'}]\n"

/*
 * Create a channel for input read files
 */

if (!params.bam) {
    if (params.rawPath==".") {
        raw_search_path = "${params.run}/*/*.fastq.gz"
    } else {
        raw_search_path = "${params.run}/*/${params.rawPath}/*.fastq.gz"
    }
    Channel
        .fromPath(raw_search_path)
        .ifEmpty { exit 1, "${params.run}/[SAMPLE_FOLDER]/${params.rawPath} - no fastq.gz file detected." }
        .toSortedList()
        .flatten()
        .map { it ->
                   def fastq = it
                   while(it.getParent().toString() != run_path) {it = it.getParent()}
                   def sample_fastq = [it.getName(), fastq]
             }
        .groupTuple()
        .set { sample_fastq }
} else {
    // Set an empty channel
    Channel.empty().set { sample_fastq }

    // Set input channel for BAM and BAI
    if (params.bamPath==".") {
        bam_search_path = "${params.run}/*/*.bam"
        bai_search_path = "${params.run}/*/*.bai"
    } else {
        bam_search_path = "${params.run}/*/${params.bamPath}/*.bam"
        bai_search_path = "${params.run}/*/${params.bamPath}/*.bai"
    }
    Channel
        .fromPath(bam_search_path)
        .ifEmpty { exit 1, "${params.run}/[SAMPLE_FOLDER]/${params.bamPath} - no bam file detected." }
        .toSortedList()
        .flatten()
        .map { it ->
                   def bam = it
                   while(it.getParent().toString() != run_path) {it = it.getParent()}
                   def sample_bam = [it.getName(), bam]
             }
        .groupTuple()
        .set { sample_bam }
    Channel
        .fromPath(bai_search_path)
        .ifEmpty { exit 1, "${params.run}/[SAMPLE_FOLDER]/${params.bamPath} - no bai file detected." }
        .toSortedList()
        .flatten()
        .map { it ->
                   def bai = it
                   while(it.getParent().toString() != run_path) {it = it.getParent()}
                   def sample_bai = [it.getName(), bai]
             }
        .groupTuple()
        .set { sample_bai }

    sample_bam.join(sample_bai).into{sample_bam_bai_1; sample_bam_bai_2}
    // seperate the channel and remove .wasp.bam from the matched files
    sample_bam_bai_1.map{ch -> [ch[0], ch[1].findAll{ it -> !it.getName().matches(".*[.]wasp[.]bam") }]}.set{sample_bam}
    sample_bam_bai_2.map{ch -> ch[2]}.set{sample_bai}
}

sample_fastq.into{sample_fastq; sample_fastq_count; sample_fastq_se_split; sample_fastq_pe_split; sample_fastq_dry}

sample_fastq_se_split.map {
    ch -> [ ch[0],
            ch[1],
            new nextflow.util.ArrayBag() ]
    }.set{sample_fastq_se_r1}

sample_fastq_pe_split.map {
    ch -> [ ch[0],
            ch[1].findAll{ it.getName().matches(".*[._]R1[._].*fastq.gz") },
            ch[1].findAll{ it.getName().matches(".*[._]R2[._].*fastq.gz") } ]
    }.set{sample_fastq_pe_r1_r2}
sample_fastq_pe_r1_r2.into{sample_fastq_pe_r1_r2;merge_pe_r1_r2}

if (!params.bam) {
    sample_count = sample_fastq_count.count().get()
} else {
    sample_bam.into {sample_bam; sample_bam_count}
    sample_count = sample_bam_count.count().get()
}
println "Sample Detected: ${sample_count}"

if(!params.bam && params.wasp){
    if (params.vcfPath==".") {
        vcf_search_path="${params.run}/*/*.vcf.gz"
    } else {
        vcf_search_path="${params.run}/*/${params.vcfPath}/*.vcf.gz"
    }
    Channel
        .fromPath(vcf_search_path)
        .ifEmpty { exit 1, "${params.run}/[SAMPLE_FOLDER]/${params.vcfPath} - no vcf.gz file detected." }
        .toSortedList()
        .flatten()
        .map { it ->
                   def vcf = it
                   while(it.getParent().toString() != run_path) {it = it.getParent()}
                   def sample_vcf = [it.getName(), vcf]
             }
        .into { sample_vcf; sample_vcf_count}
    vcf_count = sample_vcf_count.count().get()
    println "VCF Detected: ${vcf_count}"
    if (sample_count != vcf_count && !params.continueWithoutVcf) {
        // Added Sanity Check for WASP mode
        println "Error: Missing VCF.GZ."
        exit(1)
    }
}

/*
 * Print Dry Run Stats
 */
if (params.dryRun) {
    if(!params.bam) {
        if(params.wasp) {
            sample_fastq_dry.join( sample_vcf, remainder:true ).set { sample_fastq_dry }
        }
        sample_fastq_dry.subscribe { println "$it" }
    } else {
        sample_bam.into {sample_bam; sample_bam_dry}
        sample_bam_dry.subscribe { println "$it" }
    }
    exit(0)
}

/*
 * STEP 1 - FASTQ Merging and Trimmomatic
 *   (optional) --trimAdapter
 */

sample_fastq.into { merge_se_fastq; merge_pe_fastq; no_merge_fastq }

process merge_se {
    tag "${sample_name}"

    input:
    set val(sample_name), file(fastq) from merge_se_fastq

    output:
    set val(sample_name), file("${sample_name}.merged.fastq.gz") into merged_se_fastq

    when:
    !params.bam && params.trimAdapter && params.singleEnd && !(fastq instanceof nextflow.processor.TaskPath) && fastq.size()>1

    """
    zcat $fastq | pigz -c -p ${task.cpus} - > ${sample_name}.merged.fastq.gz
    """
}

process merge_pe {
    tag "${sample_name}"

    input:
    set val(sample_name), file(r1_fastq), file(r2_fastq) from merge_pe_r1_r2

    output:
    set val(sample_name), file("${sample_name}.{R1.,R2.}merged.fastq.gz") into merged_pe_fastq

    when:
    !params.bam && (params.trimAdapter &&
     !params.singleEnd &&
     !(r1_fastq instanceof nextflow.processor.TaskPath) &&
     !(r2_fastq instanceof nextflow.processor.TaskPath) &&
     r1_fastq.size()>1 &&
     r2_fastq.size()>1 )

    script:
        size = (r1_fastq instanceof nextflow.processor.TaskPath)
    """
    zcat ${r1_fastq} | pigz -c -p ${task.cpus} - > ${sample_name}.R1.merged.fastq.gz && \
    zcat ${r2_fastq} | pigz -c -p ${task.cpus} - > ${sample_name}.R2.merged.fastq.gz
    """
}

process no_merge {
    tag "${sample_name}"

    input:
    set val(sample_name), file(fastq) from no_merge_fastq

    output:
    set val(sample_name), file(fastq) into not_merged_fastq

    when:
    !params.bam && params.trimAdapter && (( params.singleEnd && (fastq instanceof nextflow.processor.TaskPath)) || (!params.singleEnd && fastq.size()==2))

    """
    """
}
merged_se_fastq.mix(merged_pe_fastq).mix(not_merged_fastq).set{ merged_fastq }

process trimmomatic {
    tag "${sample_name}"
    publishDir "${run_path}/${sample_name}/${params.outPath}/trimmomatic", mode: 'copy',
               saveAs: {filename -> if (filename.indexOf("trimmed.fastq.gz") > 0) null
               else if (filename.indexOf("trimmomatic.log") > 0) null
               else filename
               }

    input:
    set val(sample_name), file(fastq) from merged_fastq

    output:
    set val(sample_name), file("${sample_name}.{R1.,R2.,}trimmed.fastq.gz") into trimmed_fastq
    file "*trimmomatic.log" into trimmomatic_log
    file "*trimmomatic.stdout"
    file "*trimmomatic.stderr"

    when:
    !params.bam && params.trimAdapter

    script:
    paired_end = !params.singleEnd ? "PE" : "SE"
    threads = "${task.cpus}"
    outfiles = !params.singleEnd ? "${sample_name}.R1.trimmed.fastq.gz ${sample_name}.R1.unpair.trimmed.fastq.gz ${sample_name}.R2.trimmed.fastq.gz ${sample_name}.R2.unpair.trimmed.fastq.gz" : "${sample_name}.trimmed.fastq.gz"

    """
    trimmomatic \\
        $paired_end \\
        -threads $threads \\
        ${params.trimmomatic.quality_encoding} \\
        $fastq \\
        $outfiles \\
        ILLUMINACLIP:${params.adapter}:${params.trimmomatic.seed_mismatches}:${params.trimmomatic.palindrome_clip_threshold}:${params.trimmomatic.simple_clip_threshold}:${params.trimmomatic.min_adapter_length}:${params.trimmomatic.keep_both_reads} \\
        LEADING:${params.trimmomatic.leading} \\
        TRAILING:${params.trimmomatic.trailing} \\
        SLIDINGWINDOW:${params.trimmomatic.slidingwindow_size}:${params.trimmomatic.slidingwindow_quality} \\
        MINLEN:${params.trimmomatic.minlen} \\
        >  ${sample_name}.trimmomatic.stdout \\
        2> ${sample_name}.trimmomatic.stderr && \\
        sed "s/[ ][^ ]\\+fastq.gz/ ${sample_name}.fastq.gz/" ${sample_name}.trimmomatic.stderr > ${sample_name}.trimmomatic.log
    """
}

if (params.trimAdapter) {
    trimmed_fastq.into { trimmed_fastq_se_split; trimmed_fastq_pe_split; kallisto_fastq}
    trimmed_fastq_se_split.map {
        ch -> [ ch[0],
                ch[1],
                new nextflow.util.ArrayBag() ]
        }.set{trimmed_fastq_se_r1}

    trimmed_fastq_pe_split.map {
        ch -> [ ch[0],
                ch[1].findAll{ it.getName().matches(".*[._]R1[._].*fastq.gz") },
                ch[1].findAll{ it.getName().matches(".*[._]R2[._].*fastq.gz") } ]
        }.set{trimmed_fastq_pe_r1_r2}
    if (params.singleEnd) {
        trimmed_fastq_se_r1.set { star_fastq }
    } else {
        trimmed_fastq_pe_r1_r2.set { star_fastq }
    }
} else {
    if (params.singleEnd) {
        sample_fastq_se_r1.set { star_fastq }
    } else {
        sample_fastq_pe_r1_r2.set { star_fastq }
    }
}

/*
 * STEP 2 - align with STAR
 */
if(!params.bam && params.wasp){
    // WASP mode switch and merging VCF channel
    star_fastq.join(sample_vcf, remainder:true).set{star_fastq}
} else {
    star_fastq.map {
        ch -> [ ch[0],
                ch[1],
                ch[2],
                new nextflow.util.ArrayBag() ]
    }.set { star_fastq }
}
process star {
    tag "$sample_name"
    publishDir "${run_path}/${sample_name}/${params.outPath}/star", mode: 'copy',
               saveAs: {filename -> if (filename.indexOf(".bam") > 0) null
               else filename
               }

    input:
    set val(sample_name), file(r1), file(r2), file(vcf) from star_fastq

    output:
    set val(sample_name), file("*Aligned.out.bam"), val(wasp_status) into unsorted_bam
    set val(sample_name), file("*Aligned.toTranscriptome.out.bam"), val(wasp_status) into transcriptome_bam
    set val(sample_name), file("*.Unmapped.out*.gz") into unmapped_fastq
    file "*.out" into star_log
    file "*SJ.out.tab"

    when:
    !params.bam

    script:
    wasp_status = (params.wasp && (vcf.name =~ /vcf.gz/))
    wasp_option = wasp_status ? "--waspOutputMode SAMtag --varVCFfile <(zcat ${vcf})" : ""
    reads = params.singleEnd ? "${r1.join(',')}" : "${r1.join(',')} ${r2.join(',')}"
    """
    STAR \\
        --genomeDir ${params.star_index} \\
        --sjdbGTFfile ${params.gtf} \\
        --runThreadN ${task.cpus} \\
        --twopassMode ${params.star.twopassMode} ${params.star.twopass1readsN} \\
        --outSAMstrandField intronMotif \\
        --sjdbScore ${params.star.sjdbScore} \\
        --outFilterMismatchNoverReadLmax ${params.star.outFilterMismatchNoverReadLmax} \\
        --outReadsUnmapped Fastx \\
        --outSAMunmapped Within \\
        --outSAMtype BAM Unsorted \\
        --outSAMmode Full \\
        --outSAMattributes ${params.star.SAMattributes}\\
        --outSAMattrRGline ID:${sample_name} LB:${sample_name} PL:${params.star.seq_platform} PU:${params.star.seq_unit} SM:${sample_name} CN:${params.star.seq_center} DS:${params.star.rnaseq} \\
        --outFileNamePrefix ${sample_name}. \\
        --quantMode TranscriptomeSAM \\
        ${wasp_option} \\
        --readFilesCommand zcat \\
        --readFilesIn $reads && \\
    gzip *.Unmapped.out*
    """
}
unsorted_bam.into { fastqc_bam; bam_to_sort; bam_wasp_report }
fastqc_bam.map{ ch -> [ ch[0], ch[1] ] }.set{ fastqc_bam } //fastqc_bam does not need ch[2]:wasp_status

/*
 * STEP 2.1 WASP report
 *   (optional) --wasp
 */
process wasp_stats {
    tag "$sample_name"
    publishDir "${run_path}/${sample_name}/${params.outPath}/bams", mode: 'copy'

    input:
    set val(sample_name), file(bam), val(wasp_status) from bam_wasp_report

    output:
    file "*.wasp.stats"
    file "${sample_name}.wasp.bam"

    when:
    params.wasp && wasp_status

    script:
    tmp_file = "wasp_reads.sam"
    divider = params.singleEnd ? "" : "/2"
    """
    samtools view -h ${bam} | grep 'vW:i:[2-7]' > ${tmp_file}
    samtools view -H ${bam} > header
    grep "vW:i:2" ${tmp_file} | wc -l | awk '{print "vW:i:2\tmulti-mapping read\t" \$1${divider} }' >> ${sample_name}.wasp.stats
    grep "vW:i:3" ${tmp_file} | wc -l | awk '{print "vW:i:3\tvariant base in the read is N (non-ACGT)\t" \$1${divider}}' >> ${sample_name}.wasp.stats
    grep "vW:i:4" ${tmp_file} | wc -l | awk '{print "vW:i:4\tremapped read did not map\t" \$1${divider}}' >> ${sample_name}.wasp.stats
    grep "vW:i:5" ${tmp_file} | wc -l | awk '{print "vW:i:5\tremapped read multi-maps\t" \$1${divider}}' >> ${sample_name}.wasp.stats
    grep "vW:i:6" ${tmp_file} | wc -l | awk '{print "vW:i:6\tremapped read maps to a different locus\t" \$1${divider}}' >> ${sample_name}.wasp.stats
    grep "vW:i:7" ${tmp_file} | wc -l | awk '{print "vW:i:7\tread overlaps too many variants\t" \$1${divider}}' >> ${sample_name}.wasp.stats
    cat header ${tmp_file} | samtools view -b - > ${sample_name}.wasp.bam
    rm ${tmp_file}
    """
}

/*
 * STEP 3 - FASTQC
 *   (optional) --fastqc
 */
process fastqc {
    tag "${sample_name}"
    publishDir "${run_path}/${sample_name}/${params.outPath}/fastqc", mode: 'copy'

    input:
    set val(sample_name), file(reads) from fastqc_bam

    output:
    file "*_fastqc.zip"
    file "*_fastqc.html" into fastqc_report

    when:
    params.fastqc

    script:
    """
    fastqc -q $reads -t 4
    """
}

/*
 * STEP 4 Filter/Sort Bam
 *    (optional) WASP filtering --wasp
 */
process filter_sort {
    tag "$sample_name"

    input:
    set val(sample_name), file(bam), val(wasp_status) from bam_to_sort

    output:
    set val(sample_name), file("*.sorted.bam") into markdup_bam
    set val(sample_name), file("*.{no_filter,filtered}.bam") into wasp_bam

    script:
    filter_and_sort = wasp_status ? "samtools view -h ${bam} | grep -v 'vW:i:[2-7]' | samtools view -bS -@ ${task.cpus} - -o ${sample_name}.filtered.bam && samtools sort -@ ${task.cpus} -m 4G -T temp-samtool-sort -O BAM -o ${sample_name}.sorted.bam ${sample_name}.filtered.bam" : "ln -s -f ${bam} ${sample_name}.no_filter.bam && samtools sort -@ ${task.cpus} -m 4G -T temp-samtool-sort -O BAM -o ${sample_name}.sorted.bam ${bam}"
    """
    ${filter_and_sort}
    """
}
wasp_bam.into{ featureCounts_bam; enhancer_featureCounts_bam; bam_to_fastq }

process transcriptome_bam_filter {
    tag "$sample_name"

    input:
    set val(sample_name), file(bam), val(wasp_status) from transcriptome_bam

    output:
    set val(sample_name), file("*.{no_filter,filtered}.bam") into wasp_transcriptome_bam

    when:
    params.rsem

    script:
    filter_and_sort = wasp_status ? "samtools view -h ${bam} | grep -v 'vW:i:[2-7]' | samtools view -bS -@ ${task.cpus} - -o ${sample_name}.filtered.bam" : "ln -s -f ${bam} ${sample_name}.no_filter.bam"
    """
    ${filter_and_sort}
    """
}
wasp_transcriptome_bam.set { rsem_bam }

process prepare_wasp_fastq {
    tag "$sample_name"

    input:
    set val(sample_name), file(bam) from bam_to_fastq

    output:
    set val(sample_name), file("*.fastq.gz") into wasp_fastq

    when:
    params.kallisto || params.salmon || params.txrevise

    script:
    fastq = params.singleEnd? "-s single.fastq.gz" : "-1 r1.fastq.gz -2 r2.fastq.gz -s /dev/null"
    """
    samtools fastq -c 6 $fastq -0 /dev/null $bam
    """
}
wasp_fastq.into{ kallisto_fastq; salmon_fastq; txrevise_fastq }

/*
 * STEP 5 MarkDuplicates
 */
process markdup {
    tag "$sample_name"
    publishDir "${run_path}/${sample_name}/${params.outPath}", mode: 'copy',
        saveAs: {filename -> if (filename.indexOf(".DupMetrics") > 0) "qc_metrics/${filename}"
        else "bams/${filename}"
        }

    input:
    set val(sample_name), file(bam) from markdup_bam

    output:
    set val(sample_name), file("${sample_name}.bam") into markeddups_bam
    file("${sample_name}.*bai") into markeddups_bai
    file ("*.DupMetrics") optional true into dup_metrics
    file "README"

    script:
    if (params.noMarkDup) {
        markdup_cmd = ""
        ln_cmd = "ln -s ${bam} ${sample_name}.bam"
        readme_marked = ""
    } else {
        markdup_cmd = "java -jar -Xmx24g \${PICARD} MarkDuplicates I=${bam} O=${sample_name}.markdup.bam M=${sample_name}.DupMetrics"
        ln_cmd = "ln -s ${sample_name}.markdup.bam ${sample_name}.bam"
        readme_marked = "/duplicates marked"
    }
    """
    ${markdup_cmd}
    ${ln_cmd}
    samtools index ${sample_name}.bam
    ln -s ${sample_name}.bam.bai ${sample_name}.bai
    echo "[sample_name].bam: mapped and unmapped reads/(optional)wasp marked and filtered/sorted by coordinate${readme_marked} BAM" > README
    echo "[sample_name].wasp.bam (if applicable): removed wasp reads (vW:i:[2-7]) from the main BAM file" >> README
    """
}

// --bam option allow pipeline starting from markeddups_bam
if (params.bam) {
    sample_bam.into{qc_bam; bam_2nd_stage}
    sample_bai.into{qc_bai; bai_2nd_stage}
} else {
    markeddups_bam.into{qc_bam; bam_2nd_stage}
    markeddups_bai.into{qc_bai; bai_2nd_stage}
}

/*
 * STEP 6 QC
 *   (optional) --qc
 */

qc_bam.into {rnaseq_metrics_bam; alignment_metrics_bam; gcbias_metrics_bam; hsmetrics_bam; insert_size_metrics_bam}
qc_bai.into {rnaseq_metrics_bai; alignment_metrics_bai; gcbias_metrics_bai; hsmetrics_bai}

process picard_rnaseq_metrics {
    tag "$sample_name"
    publishDir "${run_path}/${sample_name}/${params.outPath}/qc_metrics", mode: 'copy'

    input:
    set val(sample_name), file(bam) from rnaseq_metrics_bam
    file(bai) from rnaseq_metrics_bai

    output:
    file("*.RNASeqMetrics") into rnaseq_metrics

    when:
    !params.bam && params.qc

    script:
    """
    java -jar \$PICARD CollectRnaSeqMetrics \\
        I=${bam} O=${sample_name}.RNASeqMetrics \\
        REF_FLAT=${params.reflat} \\
        RIBOSOMAL_INTERVALS=${params.rrna} \\
        RRNA_FRAGMENT_PERCENTAGE=0.5 \\
        STRAND=NONE
    """
}

process picard_alignment_metrics {
    tag "$sample_name"
    publishDir "${run_path}/${sample_name}/${params.outPath}/qc_metrics", mode: 'copy'

    input:
    set val(sample_name), file(bam) from alignment_metrics_bam
    file(bai) from alignment_metrics_bai

    output:
    file("*.AlignmentSummaryMetrics") into alignment_metrics

    when:
    !params.bam && params.qc

    script:
    """
    java -jar \$PICARD CollectAlignmentSummaryMetrics \\
        I=${bam} \\
        O=${sample_name}.AlignmentSummaryMetrics \\
        R=${params.fasta}
    """
}

process picard_insert_size_metrics {
    tag "$sample_name"
    publishDir "${run_path}/${sample_name}/${params.outPath}/qc_metrics", mode: 'copy'

    input:
    set val(sample_name), file(bam) from insert_size_metrics_bam

    output:
    file("*.InsertSizeMetrics") into insert_size_metrics
    file("*.InsertSize_histogram.pdf") optional true

    when:
    !params.bam && params.qc && !params.singleEnd

    script:
    """
    java -jar \$PICARD CollectInsertSizeMetrics \\
        I=${bam} \\
        O=${sample_name}.InsertSizeMetrics \\
        H=${sample_name}.InsertSize_histogram.pdf
    """
}

process picard_gcbias_metrics {
    tag "$sample_name"
    publishDir "${run_path}/${sample_name}/${params.outPath}/qc_metrics", mode: 'copy'

    input:
    set val(sample_name), file(bam) from gcbias_metrics_bam
    file(bai) from gcbias_metrics_bai

    output:
    file("*.GcBiasMetrics") into gcbias_metrics
    file("*.GcBiasSummaryMetrics") into gcbias_summary_metrics
    file("*.GcBiasMetrics.pdf") optional true

    when:
    !params.bam && params.qc

    script:
    """
    java -jar \$PICARD CollectGcBiasMetrics \\
        I=${bam} \\
        O=${sample_name}.GcBiasMetrics \\
        CHART=${sample_name}.GcBiasMetrics.pdf \\
        S=${sample_name}.GcBiasSummaryMetrics \\
        R=${params.fasta}
    """
}

process picard_hsmetrics {
    tag "$sample_name"
    publishDir "${run_path}/${sample_name}/${params.outPath}/qc_metrics", mode: 'copy'

    input:
    set val(sample_name), file(bam) from hsmetrics_bam
    file(bai) from hsmetrics_bai

    output:
    file("*.HsMetrics") into hsmetrics

    when:
    !params.bam && params.target && params.qc

    script:
    target = "${params.target}"
    bait = params.bait? "${params.bait}" : "${params.target}"
    """
    java -jar \$PICARD CollectHsMetrics \\
        I=${bam} \\
        O=${sample_name}.HsMetrics \\
        BAIT_INTERVALS=${bait} \\
        TARGET_INTERVALS=${target} \\
        R=${params.fasta}
    """
}

/*
 * STEP 7 Pathogen Analysis
 *   (optional) --pathogen
 */

process pathogen {
    tag "$sample_name"
    publishDir "${run_path}/${sample_name}/${params.outPath}/pathogen", mode: 'copy'

    input:
    set val(sample_name), file(reads) from unmapped_fastq

    output:
    file("*Log{,.final,.progress}.out") into pathogen_star_log
    file("*.pathogen.txt*")
    file("*.pathogen.bw")
    file("README")

    when:
    params.pathogen

    script:
    strand_option = params.stranded == "none" ? "--outSAMstrandField intronMotif": ""
    paired_end = !params.singleEnd ? "-p" : ""
    """
    echo -ne "GeneID\\tChr\\tStart\\tEnd\\tStrand\\n" > annotation.asf
    sed "s/\\(.*\\)|\\(.*\\)|\\(.*\\)|\\(.*\\)\\t\\([0-9]\\+\\)/\\3_\\4\\t\\1|\\2|\\3|\\4\\t1\\t\\5\\t./" ${params.pathogen_star_index}/chrNameLength.txt >> annotation.asf
    cp ${params.pathogen_star_index}/chrNameLength.txt chrNL.txt
    export LC_COLLATE=C

    STAR \\
        --genomeDir ${params.pathogen_star_index} \\
        --runThreadN ${task.cpus} \\
        $strand_option \\
        --sjdbScore ${params.star.sjdbScore} \\
        --outFilterMultimapNmax 700000 \\
        --outFilterMismatchNoverReadLmax ${params.star.outFilterMismatchNoverReadLmax} \\
        --outSAMtype BAM Unsorted \\
        --outSAMunmapped Within \\
        --outSAMattributes ${params.star.SAMattributes}\\
        --outSAMattrRGline ID:${sample_name} LB:${sample_name} PL:${params.star.seq_platform} PU:${params.star.seq_unit} SM:${sample_name} CN:${params.star.seq_center} DS:${params.star.rnaseq} \\
        --outFileNamePrefix ${sample_name}.pathogen. \\
        --readFilesCommand zcat \\
        --readFilesIn $reads && \\
    featureCounts -T ${task.cpus} \\
        -g GeneID \\
        -F SAF \\
        -M -O \\
        ${paired_end} \\
        -a annotation.asf \\
        -o ${sample_name}.pathogen.txt \\
        --donotsort \\
        ${sample_name}.pathogen.Aligned.out.bam && \\
    samtools sort -@ 12 -m 4G -T TMP ${sample_name}.pathogen.Aligned.out.bam | \\
    genomeCoverageBed -bg -ibam stdin | sort -k1,1 -k2,2n > pathogen.bg && \\
    bedGraphToBigWig pathogen.bg \\
                     chrNL.txt \\
                     ${sample_name}.pathogen.bw && \\
    rm ${sample_name}.pathogen.Aligned.out.bam pathogen.bg chrNL.txt && \\
    echo "Unmapped reads aligned to ${params.pathogen_description}" > README && \\
    echo "Counted by featureCounts -M -O" >> README
    """
}

/*
 * STEP 8 Quantification
 */

/*
 * STEP 8.1 Kallisto
 *   (optional) --kallisto
 */

process kallisto {
    tag "$sample_name"
    publishDir "${run_path}/${sample_name}/${params.outPath}/kallisto", mode: 'copy'

    input:
    set val(sample_name), file(reads) from kallisto_fastq

    output:
    file "abundance.{tsv,h5}"

    when:
    params.kallisto

    script:
    single_end_flag = !params.singleEnd ? "" : "--single --fragment-length=${params.kallisto_setting.single_end_fragment_length} --sd=${params.kallisto_setting.single_end_sd}"
    bootstrap = params.kallisto_setting.bootstrap > 0 ? "--bootstrap-samples=${params.kallisto_setting.bootstrap_samples}" : ""
    """
    kallisto quant \\
    -i ${params.kallisto_index} \\
    -o ./ \\
    $bootstrap \\
    $single_end_flag \\
    $reads
    """
}

/*
 * STEP 8.2 RSEM
 *   (optional) --rsem
 */

process rsem {
    tag "$sample_name"
    publishDir "${run_path}/${sample_name}/${params.outPath}/rsem", mode: 'copy'

    input:
    set val(sample_name), file(bam)  from rsem_bam

    output:
    file ("*.results")
    file ("*.quant.pdf")

    when:
    params.rsem

    script:
    stranded = params.stranded == "forward" ? "--forward-prob 1" : params.stranded == "reverse" ? "--forward-prob 0" : ""
    paired_end = !params.singleEnd ? "--paired-end" : ""
    """
    rsem-calculate-expression \\
        --num-threads ${task.cpus} \\
        --estimate-rspd \\
        --calc-ci \\
        --ci-memory 4096 \\
        --no-bam-output \\
        --seed 12345 \\
        --bam \\
        ${stranded} \\
        ${paired_end} \\
        ${bam} \\
        ${params.rsem_index} \\
        ${sample_name} \\
        > ${sample_name}.Log.rsem && \\
rsem-plot-model \\
        ${sample_name} \\
        ${sample_name}.quant.pdf
    """
}

/*
 * STEP 8.3.1 FeauterCounts
 *   (optional) --featureCounts
 */

process featurecounts {
    tag "$sample_name"
    publishDir "${run_path}/${sample_name}/${params.outPath}/featureCounts", mode: 'copy'

    input:
    set val(sample_name), file(bam) from featureCounts_bam

    output:
    file ("*.txt")
    file ("*.exon.geneID.txt.summary") into featureCounts_log
    file ("*.exon.transcriptID.txt.summary")
    file ("*.exon.txt.summary")
    file ("*.primary.txt.summary")

    when:
    params.featureCounts

    script:
    strand = params.stranded == "forward" ? "-s 1" : params.stranded == "reverse" ? "-s 2" : ""
    paired_end = !params.singleEnd ? "-p" : ""
    """
    featureCounts -T ${task.cpus} \\
        -t exon \\
        -g transcript_id \\
        ${paired_end} \\
        ${strand} \\
        -a ${params.gtf} \\
        -o ${sample_name}.exon.transcriptID.txt \\
        --donotsort \\
        $bam && \\
    featureCounts -T ${task.cpus} \\
        -t exon \\
        -g gene_id \\
        ${paired_end} \\
        ${strand} \\
        -a ${params.gtf} \\
        -o ${sample_name}.exon.geneID.txt \\
        --donotsort \\
        $bam && \\
    featureCounts -T ${task.cpus} \\
        -t exon \\
        -f -O \\
        ${paired_end} \\
        ${strand} \\
        -a ${params.gtf} \\
        -o ${sample_name}.exon.txt \\
        --donotsort \\
        $bam && \\
    featureCounts -T ${task.cpus} \\
        -t exon \\
        -g gene_id \\
        --primary -O \\
        ${paired_end} \\
        ${strand} \\
        -a ${params.gtf} \\
        -o ${sample_name}.primary.txt \\
        --donotsort \\
        $bam && \\
    sed -i "s/[.]no_filter[.]bam//" *.summary && \\
    sed -i "s/[.]filtered[.]bam//" *.summary
    """
}

/*
 * STEP 8.3.2 FeauterCounts for Enhancer RNA
 *   (optional) --enhancer
 */

process enhancer_featurecounts {
    tag "$sample_name"
    publishDir "${run_path}/${sample_name}/${params.outPath}/featureCounts", mode: 'copy'

    input:
    each location from "postive", "negative", "both"
    set val(sample_name), file(bam) from enhancer_featureCounts_bam

    output:
    file ("*.primary.txt")
    file ("*.primary.txt.summary")

    when:
    params.enhancer && (params.stranded == "none" && location=="both" || params.stranded != "none" && location!="both")

    script:
    strand = params.stranded == "forward" ? "-s 1" : params.stranded == "reverse" ? "-s 2" : ""
    location_strand = location == "negative" ? "strand='-'" : "strand='+'"
    paired_end = !params.singleEnd ? "-p" : ""
    annotation = location == "negative" ? "negative.saf" : "positive.saf"
    """
    echo -ne "GeneID\\tChr\\tStart\\tEnd\\tStrand\\n" > ${annotation}
    cat ${params.enhancer_bed} | awk -v OFS='\\t' -v ${location_strand} '{print \$4, \$1, \$2, \$3, strand}' >> ${annotation}

    featureCounts -T ${task.cpus} \\
        -g GeneID \\
        -F SAF \\
        --primary -O \\
        ${paired_end} \\
        ${strand} \\
        -a ${annotation} \\
        -o ${sample_name}.GeneHancerID.${location}.primary.txt \\
        --donotsort \\
        $bam
    """
}


/*
 * STEP 8.4.1 Salmon
 *   (optional) --salmon
 */

process salmon {
    tag "$sample_name"
    publishDir "${run_path}/${sample_name}/${params.outPath}/salmon", mode: 'copy'

    input:
    set val(sample_name), file(fastq) from salmon_fastq

    output:
    file ("*quant.sf")

    when:
    params.salmon

    script:
    reads = params.singleEnd? "-r ${fastq}" : "-1 ${fastq[0]} -2 ${fastq[1]}"
    orientation = params.singleEnd?"":"I" //usually inward for Illumina data
    strandness = params.stranded == "forward" ? "SF" : params.stranded == "reverse" ? "SR" : "U"
    """
    salmon quant --seqBias --useVBOpt --gcBias \\
        --libType ${orientation}${strandness} \\
        --index ${params.salmon_index} \\
        ${reads} \\
        -p ${task.cpus} \\
        -o . && \\
    mv quant.sf ${sample_name}.${params.stranded}.quant.sf
    """
}

/*
 * STEP 8.4.2 txrevise Salmon
 *   (optional) --txrevise
 */

process txrevise_salmon {
    tag "${sample_name}"
    publishDir "${run_path}/${sample_name}/${params.outPath}/salmon", mode: 'copy'

    input:
    set val(sample_name), file(fastq) from txrevise_fastq
    each index from txrevise_salmon_indices

    output:
    file ("*quant.sf")

    when:
    params.txrevise

    script:
    reads = params.singleEnd? "-r ${fastq}" : "-1 ${fastq[0]} -2 ${fastq[1]}"
    orientation = params.singleEnd?"":"I" //usually inward for Illumina data
    strandness = params.stranded == "forward" ? "SF" : params.stranded == "reverse" ? "SR" : "U"
    """
    salmon quant --seqBias --useVBOpt --gcBias \\
        --libType ${orientation}${strandness} \\
        --index ${index} \\
        ${reads} \\
        -p ${task.cpus} \\
        -o . && \\
    mv quant.sf ${sample_name}.${index.baseName}.${params.stranded}.quant.sf
    """
}

/*
 * STEP 9 Summarizing Stats
 */
star_log.mix(
    featureCounts_log,
    fastqc_report,
    rnaseq_metrics,
    alignment_metrics,
    insert_size_metrics,
    gcbias_metrics,
    gcbias_summary_metrics,
    hsmetrics,
    dup_metrics,
    trimmomatic_log,
    pathogen_star_log
    ).set{ all_metrics }

process multiqc {
    publishDir "${run_path}", mode: 'copy'

    input:
    file(star) from all_metrics.collect()

    output:
    file("multiqc_report.html")
    file("multiqc_general_stats.txt")

    when:
    !params.bam

    script:
    """
    multiqc -m star -m featureCounts -m picard -m trimmomatic . && \\
    ln -s multiqc_data/multiqc_general_stats.txt
    """
}


/*
 * STEP 10 Additional Tools
 */

bam_2nd_stage.into {leafcutter_bam; bam_file_to_bigwig; verifybamid_bam; qtltools_mbv_bam; bam_coverage_bam; ase_read_counter_bam}
bai_2nd_stage.into {leafcutter_bai; bai_file_to_bigwig; verifybamid_bai; qtltools_mbv_bai; bam_coverage_bai; ase_read_counter_bai}
 
/* Leafcutter
 *   (optional) --leafcutter
 */

process leafcutter {
    tag "$sample_name"
    publishDir "${run_path}/${sample_name}/${params.outPath}/leafcutter", mode: 'copy'

    input:
    set val(sample_name), file(bam) from leafcutter_bam
    file(bai) from leafcutter_bai

    output:
    file ("*.junc")

    when:
    params.leafcutter
    stranded = params.stranded == "none" ? 0 : params.stranded == "forward" ? 1 : params.stranded == "reverse" ? 2 : "ERROR"

    script:
    """
    samtools view -b -q ${params.leafcutter_setting.quality_filter_threshold} \\
    $bam > ${sample_name}.Aligned.Quality.Sorted.bam && \\
    regtools junctions extract -a 8 -m 50 -M 500000 -s ${stranded} -o ${sample_name}.${params.stranded}.output.junc $bam && \\
    echo ${sample_name}.${params.stranded}.output.junc >> junction_files.fofn
    #python \$LEAFCUTTER_CLUSTER_REGTOOLS -j junction_files.fofn -m 50 -o ${sample_name}.${params.stranded}.cluster -l 500000
    """
}

/*
 * BigWig generation
 *   (optional) --bigwig
 *   (optional) --bamCoverage
 */

process bam_to_bigwig {
    tag "$sample_name"
    publishDir "${run_path}/${sample_name}/${params.outPath}/bigwig", mode: 'copy'

    input:
    each strand from "forward", "reverse", "none"
    set val(sample_name), file(bam) from bam_file_to_bigwig

    output:
    file ("*.bw")

    when:
    params.bigwig && ((params.stranded == "none" && strand == "none") || (params.stranded != "none" && strand != "none"))

    script:
    bed_graph_strand = strand == "forward" ? "-strand +" : strand == "reverse" ? "-strand -" : ""
    bg = strand == "forward" ? "${sample_name}.forward.bg" : strand == "reverse" ? "${sample_name}.reverse.bg" : "${sample_name}.unstranded.bg"
    bw = strand == "forward" ? "${sample_name}.forward.bw" : strand == "reverse" ? "${sample_name}.reverse.bw" : "${sample_name}.unstranded.bw"
    """
    export LC_COLLATE=C
    genomeCoverageBed -bg ${bed_graph_strand} -ibam ${bam} | sort -k1,1 -k2,2n > ${bg}
    cp ${params.star_index}/chrNameLength.txt chrNL.txt
    bedGraphToBigWig ${bg} \\
                     chrNL.txt \\
                     ${bw}
    """
}

process bam_coverage {
    //TODO
    tag "$sample_name"
    publishDir "${run_path}/${sample_name}/${params.outPath}/bam_coverage", mode: 'copy'

    input:
    each strand from "forward", "reverse", "none"
    set val(sample_name), file(bam) from bam_coverage_bam
    file (bai) from bam_coverage_bai

    output:
    file ("*.bw")
    
    when:
    params.bamCoverage && ((params.stranded == "none" && strand == "none") || (params.stranded != "none" && strand != "none"))
    
    script:
    filterRNAstrand = params.stranded != "none" ? "--filterRNAstrand ${params.stranded}" : "" 
    """
    bamCoverage \\
        --bam ${bam} \\
        --numberOfProcessors ${task.cpus} \\
        --outFileFormat bigwig \\
        ${filterRNAstrand} \\
        --outFileName ${sample_name}.${strand}.bw 
    """
}

process verifybamid {
    tag "$sample_name"
    publishDir "${run_path}/${sample_name}/${params.outPath}/verifyBamID", mode: 'copy'

    input:
    set val(sample_name), file(bam) from verifybamid_bam
    file(bai) from verifybamid_bai

    output:
    file("${sample_name}.*SM")
    file("${sample_name}.*RG")

    when:
    params.globalVCF && params.verifyBamID
    
    script:
    """
    verifyBamID --bam $bam --vcf ${params.globalVCF} --out ${sample_name} --ignoreRG
    """
}

process qtltools_mbv {
    tag "$sample_name"
    publishDir "${run_path}/${sample_name}/${params.outPath}/qtltools_mbv", mode: 'copy'

    input:
    set val(sample_name), file(bam) from qtltools_mbv_bam
    file(bai) from qtltools_mbv_bai

    output:
    file("${sample_name}.match.txt")

    when:
    params.globalVCF && params.QTLtools_mbv
    
    script:
    """
    QTLtools mbv --bam $bam --vcf ${params.globalVCF} --filter-mapping-quality 150 --out ${sample_name}.match.txt
    """
}

/*
process ase_read_counter {
    tag "$sample_name"
    publishDir "${run_path}/${sample_name}/${params.outPath}/ase_read_counter", mode: 'copy'

    input:
    set val(sample_name), file(bam) from ase_read_counter_bam
    file(bai) from ase_read_counter_bai

    output:
    file ("${sample_name}.csv")

    when:
    params.ASEReadCounterSiteVCF

    script:
    """
    java -jar \$GATK_JAR ASEReadCounter \\
    -R ${params.fasta} \\
    -O ${sample_name}.ASEReadCount.tsv \\
    -I ${bam} \\
    -V ${params.ASEReadCounterSiteVCF} \\
    -DF NotDuplicateReadFilter
    """
}
*/
