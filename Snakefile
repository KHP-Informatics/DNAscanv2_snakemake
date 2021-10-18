configfile: "config.yaml"

import os
import os.path
import sys
import glob
import gzip
import bgzip

#INPUT OPTIONS
##bam, cram, fastq and vcf options
format = config["INPUT_FORMAT"]
sample_name = config["SAMPLE_NAMES"]
alsgenescanner = config["ALSGENESCANNER"]
reference = config["REFERENCE_VERSION"]
paired = config["READ_TYPE"]
input_dir = config["INPUT_FILE_DIR"]
vcf = config["ADDITIONAL_VCF_FILE_DIR"]
filter_string = config["FILTER_VARIANTS"]
debug = config["DEBUG"]
exome = config["EXOME"]
BED = config["USE_BED"]
use_gene_list = config["USE_GENE_LIST_TO_GENERATE_BED"]
RG = config["ADD_READ_GROUP"]
rm_dup = config["REMOVE_DUPLICATES"]
memory = config["MEM_GB"] * 1000

#MAIN OUTPUT DIRECTORIES
log_dir = config["OUT_DIR"] + "logs/"
reports_dir = config["OUT_DIR"] + "reports/"
results_dir = config["OUT_DIR"] + "results/"

#OPTIONS
alignment = config["ALIGNMENT"]
variantcalling = config["SNP_INDEL"]
SV = config["STRUCTURAL_VARIANT"]
expansion = config["EXPANSION"]
STR = config["SHORT_TANDEM_REPEAT"]
genotypeSTR = config["genotype_SHORT_TANDEM_REPEAT"]
MEI = config["MOBILE_ELEMENT_INSERTION"]
annotation = config["ANNOTATION"]
annovar_protocols = config["ANNOVAR_PROTOCOLS"]
annovar_operations = config["ANNOVAR_OPERATIONS"]
sequencing_report = config["SEQUENCING_REPORT"]
alignment_report = config["ALIGNMENT_REPORT"]
calls_report = config["SNP_INDEL_CALLS_REPORT"]
results_report = config["ANNOTATION_RESULTS_REPORT"]
virus = config["VIRUS"]
bacteria = config["BACTERIA"]
microbes = config["CUSTOM_MICROBES"]
RG_ID = config["ID"]
RG_LB = config["LIBRARY"]
RG_PL = config["PLATFORM"]
RG_PU = config["PLATFORM_UNIT"]
RG_SM = config["SAMPLE"]

#PATHS
path_bed = config["BED_FILE"]
path_gene_list = config["GENE_LIST"]
path_reference = config["REFERENCE_FILE"]
path_melt = config["MELT_DIR"]

#FILES FOR REFERENCE VERSIONS BASED ON ANALYSIS OPTIONS
if reference == "hg19" or reference == "grch37":
    path_expansionHunter_catalog = "resources/repeats/hg19_variant_catalog.json"
    path_delly_exclude_regions = "resources/delly_hg19.excl.tsv"
    melt_zipped_files = path_melt + "me_refs/1KGP_Hg19/*zip"
    melt_bed = path_melt + "add_bed_files/1KGP_Hg19/hg19.genes.bed"
    annovar_ref_version = "hg19"
    annotsv_ref_version = "GRCh37"
    if alsgenescanner == "true":
        alsgenescanner_bed = "resources/alsgenescanner/als_gene_scanner_hg19.bed"
    if reference == "grch37":
        os.system("zcat resources/exome_hg19.bed.gz | sed 's/chr//g' | bgzip -c > resources/exome_grch37.bed.gz")
        os.system("zcat resources/hg19_gene_db.txt.gz | sed 's/chr//g' | bgzip -c > resources/grch37_gene_db.txt.gz")
        os.system("cp resources/hg19_gene_names.txt.gz resources/grch37_gene_names.txt.gz")


if reference == "hg38" or reference == "grch38":
    path_expansionHunter_catalog = "resources/repeats/hg38_variant_catalog.json"
    path_delly_exclude_regions = "resources/delly_hg38.excl.tsv"
    melt_zipped_files = path_melt + "me_refs/Hg38/*zip"
    melt_bed = path_melt + "add_bed_files/Hg38/Hg38.genes.bed"
    annovar_ref_version = "hg38"
    annotsv_ref_version = "GRCh38"
    if alsgenescanner == "true":
        alsgenescanner_bed = "resources/alsgenescanner/als_gene_scanner_hg38.bed"
    if reference == "grch38":
        os.system("zcat resources/exome_hg38.bed.gz | sed 's/chr//g' | bgzip -c > resources/exome_grch38.bed.gz")
        os.system("zcat resources/hg38_gene_db.txt.gz | sed 's/chr//g' | bgzip -c > resources/grch38_gene_db.txt.gz")
        os.system("cp resources/hg38_gene_names.txt.gz resources/grch38_gene_names.txt.gz")

#ADAPTING INPUT FILE FORMATS
if format == "fastq" and paired == "single":
    expand(input_dir + "{sample}.1.fq.gz", sample=sample_name)
if format == "fastq" and paired == "paired":
    expand(input_dir + "{sample}.1.fq.gz", sample=sample_name) + expand(input_dir + "{sample}.2.fq.gz", sample=sample_name)
if format == "sam":
    expand(input_dir + "{sample}.sam", sample=sample_name)
if format == "bam":
    expand(input_dir + "{sample}.bam", sample=sample_name)
if format == "cram":
    expand(input_dir + "{sample}.cram", sample=sample_name)

if alsgenescanner == "true":
    alsgene_annovar_protocols = "refGene,dbnsfp33a,clinvar_20210501,intervar_20180118"
    alsgene_annovar_operations = "g,f,f,f"
    path_gene_list = ""
    filter_string = "false"
    BED = "true"
    annotation = "true"
    variantcalling = "true"
    SV = "true"
    expansion = "true"

if rm_dup == "true":
    if exome == "true":
        samblaster_cmq = "samblaster --ignoreUnmated |"
    else:
        samblaster_cmq = "samblaster |"
        samblaster_bwa = "samblaster --ignoreUnmated |"

if config["USE_OWN_TEMP_DIR"] == "true":
    tmp_dir = config["TEMPORARY_DIR"]
else:
    tmp_dir = results_dir + "/tmp"

if RG:
    rg_option_hisat2 = " --rg-id %s --rg LB:%s --rg PL:%s --rg PU:%s --rg SM:%s" % (RG_ID, RG_LB, RG_PL, RG_PU, RG_SM)
    rg_option_bwa = " -R '@RG\\tID:%s\\tLB:%s\\tPL:%s\\tPU:%s\\tSM:%s' " % (RG_ID, RG_LB, RG_PL, RG_PU, RG_SM)
else:
    rg_option_hisat2 = ""
    rg_option_bwa = ""

if config["HISAT_CUSTOM_OPTIONS"] == "None":
    config["HISAT_CUSTOM_OPTIONS"] = ""

if config["BWA_CUSTOM_OPTIONS"] == "None":
    config["BWA_CUSTOM_OPTIONS"] = ""

if config["MELT_CUSTOM_OPTIONS"] == "None":
    config["MELT_CUSTOM_OPTIONS"] = ""

if config["ANNOTSV_CUSTOM_OPTIONS"] == "None":
    config["ANNOTSV_CUSTOM_OPTIONS"] = ""

rule all:
    input:
        expand(results_dir + "{sample}/custom.bed", sample=sample_name) if use_gene_list == "true" and path_gene_list else [] +
        expand(results_dir + "{sample}/{sample}_sorted_aligned.bam", sample=sample_name) if format == "fastq" and alignment == "true" or alsgenescanner == "true" else [] +
        expand(results_dir + "{sample}/{sample}_sorted_aligned.bam.bai", sample=sample_name) if format == "fastq" and alignment == "true" or alsgenescanner == "true" else [] +
        expand(results_dir + "{sample}/{sample}_sorted.bam", sample=sample_name) if format == "sam" else [] +
        expand(results_dir + "{sample}/{sample}_sorted.bam.bai", sample=sample_name) if format == "sam" else [] +
        expand(results_dir + "{sample}/{sample}_delly.bam", sample=sample_name) if format == "cram" and SV == "true" else [] +
        expand(results_dir + "{sample}/{sample}_delly.bam.bai", sample=sample_name) if format == "cram" and SV == "true" else [] +
        expand(results_dir + "{sample}/{sample}_sorted.vcf.gz", sample=sample_name) if variantcalling == "true" and filter_string != "true" or alsgenescanner == "true" else [] +
        expand(results_dir + "{sample}/{sample}_sorted_filtered.vcf.gz", sample=sample_name) if variantcalling == "true" and filter_string == "true" else [] +
        expand(results_dir + "{sample}/{sample}_expansions.vcf.gz", sample=sample_name) if expansion == "true" or alsgenescanner == "true" else [] +
        expand(results_dir + "{sample}/{sample}_expansions.vcf.gz.tbi", sample=sample_name) if expansion == "true" or alsgenescanner == "true" else [] +
        expand(results_dir + "{sample}/{sample}_expansiondenovo.str_profile.json", sample=sample_name) if STR == "true" else [] +
        expand(results_dir + "{sample}/{sample}_genotypeSTRinput.txt", sample=sample_name) if STR == "true" and genotypeSTR == "true" else [] +
        expand(results_dir + "{sample}/{sample}_EHDN_variant_catalog.json", sample=sample_name) if STR == "true" and genotypeSTR == "true" else [] +
        expand(results_dir + "{sample}/{sample}_EHDNexpansions.vcf.gz", sample=sample_name) if STR == "true" and genotypeSTR == "true" else [] +
        expand(results_dir + "{sample}/{sample}_EHDNexpansions.vcf.gz.tbi", sample=sample_name) if STR == "true" and genotypeSTR == "true" else [] +
        expand(results_dir + "{sample}/{sample}_merged_SV.vcf.gz", sample=sample_name) if SV == "true" or alsgenescanner == "true" else [] +
        expand(results_dir + "{sample}/{sample}_merged_SV.vcf.gz.tbi", sample=sample_name) if SV == "true" or alsgenescanner == "true" else [] +
        expand(results_dir + "{sample}/{sample}_MEI.vcf.gz", sample=sample_name) if MEI == "true" else [] +
        expand(results_dir + "{sample}/{sample}_MEI.vcf.gz.tbi", sample=sample_name) if MEI == "true" else [] +
        expand(results_dir + "{sample}/{sample}_SV_MEI.merged.vcf.gz", sample=sample_name) if SV == "true" and MEI == "true" else [] +
        expand(results_dir + "{sample}/{sample}_SV_MEI.merged.vcf.gz.tbi", sample=sample_name) if SV == "true" and MEI == "true" else [] +
        expand(results_dir + "{sample}/unaligned_reads.fastq.gz", sample=sample_name) if (virus or bacteria or microbes) == "true" else [] +
        expand(results_dir + "{sample}/{sample}_virus_stats.txt", sample=sample_name) if virus == "true" else [] +
        expand(reports_dir + "{sample}/{sample}_virus_report.txt", sample=sample_name) if virus == "true" else [] +
        expand(results_dir + "{sample}/{sample}_bacteria_stats.txt", sample=sample_name) if bacteria == "true" else [] +
        expand(reports_dir + "{sample}/{sample}_bacteria_report.txt", sample=sample_name) if bacteria == "true" else [] +
        expand(results_dir + "{sample}/{sample}_microbes_stats.txt", sample=sample_name) if microbes == "true" else [] +
        expand(reports_dir + "{sample}/{sample}_microbes_report.txt", sample=sample_name) if microbes == "true" else [] +
        expand(results_dir + "{sample}/{sample}_SNPindel_annotated.vcf.gz", sample=sample_name) if variantcalling == "true" and annotation == "true" or alsgenescanner == "true" else [] +
        expand(results_dir + "{sample}/{sample}_SNPindel_annotated.vcf.gz.tbi", sample=sample_name) if variantcalling == "true" and annotation == "true" or alsgenescanner == "true" else [] +
        expand(results_dir + "{sample}/{sample}_SNPindel_annotated.txt", sample=sample_name) if (variantcalling and annotation and alsgenescanner) == "true" else [] +
        expand(results_dir + "{sample}/{sample}_expansions_annotated.vcf.gz", sample=sample_name) if expansion == "true" and annotation == "true" or alsgenescanner == "true" else [] +
        expand(results_dir + "{sample}/{sample}_expansions_annotated.vcf.gz.tbi", sample=sample_name) if expansion == "true" and annotation == "true" or alsgenescanner == "true" else [] +
        expand(results_dir + "{sample}/{sample}_STR_annotated.vcf.gz", sample=sample_name) if (STR and genotypeSTR and annotation) == "true" else [] +
        expand(results_dir + "{sample}/{sample}_STR_annotated.vcf.gz.tbi", sample=sample_name) if (STR and genotypeSTR and annotation) == "true" else [] +
        expand(results_dir + "{sample}/{sample}_SV_annotated.tsv", sample=sample_name) if SV == "true" and annotation == "true" or alsgenescanner == "true" else [] +
        expand(results_dir + "{sample}/{sample}_MEI_annotated.tsv", sample=sample_name) if MEI == "true" and annotation == "true" else [] +
        expand(results_dir + "{sample}/{sample}_SV_MEI_annotated.tsv", sample=sample_name) if (SV and MEI and annotation) == "true" else [] +
        expand(reports_dir + "{sample}/{sample}_alignment_flagstat.txt", sample=sample_name) if alignment == "true" and alignment_report == "true" else [] +
        expand(reports_dir + "{sample}/{sample}_alignment_stats.txt", sample=sample_name) if alignment == "true" and alignment_report == "true" else [] +
        expand(reports_dir + "{sample}/{sample}_sequencing_report.txt", sample=sample_name) if format == "fastq" and sequencing_report == "true" else []+
        expand(reports_dir + "{sample}/{sample}_calls_vcfstats.txt", sample=sample_name) if variantcalling == "true" and calls_report == "true" else [] +
        expand(reports_dir + "multiqc_report.html", sample=sample_name) if (alignment_report and sequencing_report and calls_report) == "true" else [] +
        expand(reports_dir + "{sample}/{sample}_annovar_SNPindel.txt", sample=sample_name) if (variantcalling and annotation and results_report) == "true" or alsgenescanner == "true" else [] +
        expand(reports_dir + "{sample}/{sample}_annovar_expansions.txt", sample=sample_name) if (expansion and annotation and results_report) == "true" or alsgenescanner == "true" else [] +
        expand(reports_dir + "{sample}/{sample}_annovar_STR.txt", sample=sample_name) if (STR and genotypeSTR and annotation and results_report) == "true" else [] +
        expand(reports_dir + "{sample}/{sample}_SV_MEI_annotated.html", sample=sample_name) if (SV and MEI and annotation and results_report) == "true" else [] +
        expand(reports_dir + "{sample}/{sample}_SV_annotated.html", sample=sample_name) if (SV and annotation and results_report) == "true" or alsgenescanner == "true" else [] +
        expand(reports_dir + "{sample}/{sample}_MEI_annotated.html", sample=sample_name) if (MEI and annotation and results_report) == "true" else [] +
        expand(reports_dir + "{sample}/{sample}_all_variants.tsv", sample=sample_name) if (results_report and variantcalling and SV and (MEI or expansion or genotypeSTR)) == "true" else [] +
        expand(reports_dir + "{sample}/alsgenescanner/{sample}_alsgenescanner_all.txt", sample=sample_name) if alsgenescanner == "true" else [] +
        expand(reports_dir + "{sample}/alsgenescanner/{sample}_alsgenescanner_alsod.txt", sample=sample_name) if alsgenescanner == "true" else [] +
        expand(reports_dir + "{sample}/alsgenescanner/{sample}_alsgenescanner_clinvar.txt", sample=sample_name) if alsgenescanner == "true" else [] +
        expand(reports_dir + "{sample}/alsgenescanner/{sample}_alsgenescanner_manual_review.txt", sample=sample_name) if alsgenescanner == "true" else [] +
        expand(reports_dir + "{sample}/alsgenescanner/{sample}_alsgenescanner_all_ranked.txt", sample=sample_name) if alsgenescanner == "true" else []

rule custombed:
    input:
        path_gene_list
    output:
        matched_genes = results_dir + "{sample}/matched_genes.txt",
        unmatched_genes = results_dir + "{sample}/unmatched_genes.txt",
        matched_genes_codes = results_dir + "{sample}/matched_genes_codes.txt",
        custom_temp = results_dir + "{sample}/custom_tmp.bed",
        custom_sorted = results_dir + "{sample}/custom_sorted.bed",
        custom_bed = results_dir + "{sample}/custom.bed",
    params:
        use_gene_list
    conda:
        "envs/bedtools.yaml"
    log:
        log_dir + "{sample}/custombed.log"
    threads: config["NUMBER_CPU"]
    resources:
        mem_mb = memory
    shell:
        """
        gene_list="{input.path_gene_list}"
        use_list="{params[0]}"
        if [[ "$gene_list" ]] && [[ "$use_list" == "True" ]]; then
            zgrep -iwf {input[0]} resources/{config[REFERENCE_VERSION]}_gene_names.txt.gz | awk '{{print $2}}' > {output.matched_genes}
            zgrep -viwf {output.matched_genes} {input[0]} > {output.unmatched_genes}
            zgrep -iwf {input[0]} resources/{config[REFERENCE_VERSION]}_gene_names.txt.gz | awk '{{print $1}}' > {output.matched_genes_codes}
            zgrep -wf {output.matched_genes_codes} resources/{config[REFERENCE_VERSION]}_gene_db.txt | awk '{{i=1; while (i<= int($8)) {n=split($9,a,/,/);n=split($10,b,/,/); print $2\"\t\"a[i]\"\t\"b[i]; i+=1}}}' > {output.custom_temp}
            bedtools sort -i {output.custom_temp} > {output.custom_sorted}
            bedtools merge -i {output.custom_sorted} > {output.custom_bed}
            rm {output.custom_sorted} {output.custom_temp}
        fi
        """

rule alignment:
    input:
        fastq1 = input_dir + "{sample}.1.fq.gz" if format == "fastq" else [],
        fastq2 = input_dir + "{sample}.2.fq.gz" if format == "fastq" and paired == "paired" else []
    output:
        bam_file = results_dir + "{sample}/{sample}_sorted_aligned.bam",
        bam_file_index = results_dir + "{sample}/{sample}_sorted_aligned.bam.bai"
    params:
        format,
        alignment,
        paired,
        variantcalling,
        SV,
        MEI,
        STR,
        genotypeSTR,
        expansion,
        hisat2_bam = results_dir + "{sample}/{sample}_hisat2.bam",
        unaligned_reads = results_dir + "{sample}/{sample}_unaligned_reads.fq",
        bwa_bam = results_dir + "{sample}/{sample}_bwa.bam",
        header = results_dir + "{sample}/header.txt",
        sample = "{sample}",
        samblaster = "samblaster --ignoreUnmated" if (rm_dup and exome) == "true" else ("samblaster |" if rm_dup == "true" else []),
        out_dir = results_dir,
        rg_hisat2 = rg_option_hisat2,
        rg_bwa = rg_option_bwa,
        tmp = tmp_dir,
        bwa_samblaster = "samblaster --ignoreUnmated |"
    conda:
        "envs/alignmentfast.yaml"
    resources:
        mem_mb = memory
    log:
        log_dir + "{sample}/alignmentSNPindel.log"
    shell:
        """
        format="{params[0]}"
        paired="{params[2]}"
        alignment="{params[1]}"
        variantcalling="{params[3]}"
        SV="{params[4]}"
        MEI="{params[5]}"
        STR="{params[6]}"
        genotype="{params[7]}"
        expansion="{params[8]}"
        if [[ "$format" == "fastq" ]] && [[ "$alignment" == "True" ]]; then
            if [[ "$variantcalling" == "True" ]]; then
                if [[ "$SV" == "True" || "$MEI" == "True" || "$STR" == "True" || "$genotypeSTR" == "True" || "$expansion" == "True" ]]; then
                    if [[ "$paired" == "paired" ]]; then
                        hisat2 {config[HISAT_CUSTOM_OPTIONS]} --no-spliced-alignment -p {config[NUMBER_CPU]} -x {config[HISAT2_INDEX]} -1 {input.fastq1} -2 {input.fastq2} | {params.samblaster} samtools view -@ {config[NUMBER_CPU]} -Sb - | sambamba sort -t {config[NUMBER_CPU]} --tmpdir={params.tmp} -o {output.bam_file} /dev/stdin
                        samtools index -@ {config[NUMBER_CPU]} {output.bam_file}
                    else
                        hisat2 {config[HISAT_CUSTOM_OPTIONS]} --no-spliced-alignment -p {config[NUMBER_CPU]} -x {config[HISAT2_INDEX]} -U {input[0]} | {params.samblaster} samtools view -@ {config[NUMBER_CPU]} -Sb - | sambamba sort -t {config[NUMBER_CPU]} --tmpdir={params.tmp} -o {output.bam_file} /dev/stdin
                        samtools index -@ {config[NUMBER_CPU]} {output.bam_file}
                    fi
                else
                    if [[ "$paired" == "paired" ]]; then
                        hisat2 {config[HISAT_CUSTOM_OPTIONS]} {params.rg_hisat2} --no-softclip --no-spliced-alignment -p {config[NUMBER_CPU]} -x {config[HISAT2_INDEX]} -1 {input.fastq1} -2 {input.fastq2} | {params.samblaster} samtools view -@ {config[NUMBER_CPU]} -Sb - | sambamba sort -t {config[NUMBER_CPU]} --tmpdir={params.tmp} -o {params.hisat2_bam} /dev/stdin
                        samtools index -@ {config[NUMBER_CPU]} {params.hisat2_bam}
                        samtools view -@ {config[NUMBER_CPU]} -bhf 4 {params.hisat2_bam} | samtools bam2fq - > {params.unaligned_reads}
                        bwa mem {config[BWA_CUSTOM_OPTIONS]} {params.rg_bwa} -t {config[NUMBER_CPU]} {config[BWA_INDEX]} {params.unaligned_reads} | {params.bwa_samblaster} samtools view -@ {config[NUMBER_CPU]} -Sb - | sambamba sort -t {config[NUMBER_CPU]} --tmpdir={params.tmp} -o {params.bwa_bam} /dev/stdin
                        samtools index -@ {config[NUMBER_CPU]} {params.bwa_bam}
                        samtools view -H {params.hisat2_bam} > {params.header}
                        samtools merge -c -@ {config[NUMBER_CPU]} -f -h {params.header} {output.bam_file} {params.hisat2_bam} {params.bwa_bam}
                        rm {params.unaligned_reads} {params.header} {params.hisat2_bam} {params.bwa_bam} {params.out_dir}{params.sample}/{params.sample}_hisat2.bam.bai {params.out_dir}{params.sample}/{params.sample}_bwa.bam.bai
                        samtools index -@ {config[NUMBER_CPU]} {output.bam_file}
                    else
                        hisat2 {config[HISAT_CUSTOM_OPTIONS]} {params.rg_hisat2} --no-softclip --no-spliced-alignment -p {config[NUMBER_CPU]} -x {config[HISAT2_INDEX]} -1 {input.fastq1} -2 {input.fastq2} | {params.samblaster} samtools view -@ {config[NUMBER_CPU]} -Sb - | sambamba sort -t {config[NUMBER_CPU]} --tmpdir={params.tmp} -o {params.hisat2_bam} /dev/stdin
                        samtools index -@ {config[NUMBER_CPU]} {params.hisat2_bam}
                        samtools view -@ {config[NUMBER_CPU]} -bhf 4 {params.hisat2_bam} | samtools bam2fq - > {params.unaligned_reads}
                        bwa mem {config[BWA_CUSTOM_OPTIONS]} {params.rg_bwa} -t {config[NUMBER_CPU]} {config[BWA_INDEX]} {params.unaligned_reads} | {params.bwa_samblaster} samtools view -@ {config[NUMBER_CPU]} -Sb - | sambamba sort -t {config[NUMBER_CPU]} --tmpdir={params.tmp} -o {params.bwa_bam} /dev/stdin
                        samtools index -@ {config[NUMBER_CPU]} {params.bwa_bam}
                        samtools view -H {params.hisat2_bam} > {params.header}
                        samtools merge -c -@ {config[NUMBER_CPU]} -f -h {params.header} {output.bam_file} {params.hisat2_bam} {params.bwa_bam}
                        rm {params.unaligned_reads} {params.header} {params.hisat2_bam} {params.bwa_bam} {params.out_dir}{params.sample}/{params.sample}_hisat2.bam.bai {params.out_dir}{params.sample}/{params.sample}_bwa.bam.bai
                        samtools index -@ {config[NUMBER_CPU]} {output.bam_file}
                    fi
                fi
            fi
        fi
        """

rule sam2bam:
    input:
        input_file = input_dir + "{sample}.sam" if format == "sam" else []
    output:
        bam_file = results_dir + "{sample}/{sample}_sorted.bam",
        bam_file_index = results_dir + "{sample}/{sample}_sorted.bam.bai"
    params:
        format
    conda:
        "envs/samtools.yaml"
    log:
        log_dir + "{sample}/samtobam.log"
    resources:
        mem_mb = memory
    shell:
        """
        format="{params[0]}"
        if [[ "$format" == "sam" ]]; then
            samtools view -Sb {input.input_file} > {output.bam_file}
            samtools index -@ {config[NUMBER_CPU]} {output.bam_file}
        fi
        """

rule variantcallingfiltered:
    input:
        path_reference,
        bed = rules.custombed.output.custom_bed if use_gene_list == "true" else (path_bed if BED == "true" and not use_gene_list else (alsgenescanner_bed if alsgenescanner == "true" else [])),
        bam_file = input_dir + "{sample}.bam" if config["INPUT_FORMAT"] == "bam" else (input_dir + "{sample}.cram" if config["INPUT_FORMAT"] == "cram" else (rules.sam2bam.output.bam_file if config["INPUT_FORMAT"] == "sam" else (rules.alignment.output.bam_file if config["INPUT_FORMAT"] == "fastq" else [])))
    output:
        variant_results_file_filtered = results_dir + "{sample}/{sample}_sorted_filtered.vcf.gz",
        variant_results_file_filtered_index = results_dir + "{sample}/{sample}_sorted_filtered.vcf.gz.tbi"
    log:
        log_dir + "{sample}/variantcallingfiltered.log"
    conda:
        "envs/variantcalling.yaml"
    resources:
        mem_mb = memory
    params:
        variantcalling,
        exome,
        paired,
        BED,
        filter_string,
        out_dir = results_dir,
        sample = "{sample}",
        temp_bed = results_dir + "{sample}/temp.bed.gz",
        sorted_bed = results_dir + "{sample}/sorted.bed.gz",
        temp_results_filtered = results_dir + "{sample}/{sample}_sorted_unfiltered.vcf.gz",
        temp_results_filtered_index = results_dir + "{sample}/{sample}_sorted_unfiltered.vcf.gz.tbi"
    shell:
        """
        variantcalling="{params[0]}"
        exome="{params[1]}"
        paired="{params[2]}"
        bed="{params[3]}"
        filter_string="{params[4]}"
        if [[ "$variantcalling" == "True" ]] && [[ "$paired" == "paired" ]]; then
            if [[ "$filter_string" == "True" ]]; then
                {config[STRELKA_DIR]}bin/configureStrelkaGermlineWorkflow.py --bam {input.bam_file} --referenceFasta {input[0]} --runDir {params.out_dir}/{params.sample}/strelka
                {params.out_dir}{params.sample}/strelka/runWorkflow.py -j {config[NUMBER_CPU]} -m local
                mv {params.out_dir}{params.sample}/strelka/results/variants/genome.S1.vcf.gz {params.temp_results_filtered}
                mv {params.out_dir}{params.sample}/strelka/results/variants/genome.S1.vcf.gz.tbi {params.temp_results_filtered_index}
                rm -r {params.out_dir}{params.sample}/strelka
                bcftools filter -i '{config[VARIANT_FILTER_STRING]}' {params.temp_results_filtered} | bgzip -c > {output.variant_results_file_filtered} ; tabix -fp vcf {output.variant_results_file_filtered}

                if [[ "$exome" == "True" ]]; then
                    {config[STRELKA_DIR]}bin/configureStrelkaGermlineWorkflow.py --bam {input.bam_file} --referenceFasta {input[0]} --runDir {params.out_dir}/{params.sample}/strelka --exome
                    {params.out_dir}{params.sample}/strelka/runWorkflow.py -j {config[NUMBER_CPU]} -m local
                    mv {params.out_dir}{params.sample}/strelka/results/variants/genome.S1.vcf.gz {params.temp_results_filtered}
                    mv {params.out_dir}{params.sample}/strelka/results/variants/genome.S1.vcf.gz.tbi {params.temp_results_filtered_index}
                    bcftools filter -i '{config[VARIANT_FILTER_STRING]}' {params.temp_results_filtered} | bgzip -c > {output.variant_results_file_filtered} ; tabix -fp vcf {output.variant_results_file_filtered}

                    if [[ "$bed" == "True" ]]; then
                        bgzip -c {input.bed} > {params.temp_bed}
                        sortBed -i {params.temp_bed} | bgzip -c > {params.sorted_bed}
                        tabix -p bed {params.sorted_bed}
                        {config[STRELKA_DIR]}bin/configureStrelkaGermlineWorkflow.py --bam {input.bam_file} --referenceFasta {input[0]} --runDir {params.out_dir}/{params.sample}/strelka --callRegions {params.sorted_bed}
                        {params.out_dir}{params.sample}/strelka/runWorkflow.py -j {config[NUMBER_CPU]} -m local
                        mv {params.out_dir}{params.sample}/strelka/results/variants/genome.S1.vcf.gz {params.temp_results_filtered}
                        mv {params.out_dir}{params.sample}/strelka/results/variants/genome.S1.vcf.gz.tbi {params.temp_results_filtered_index}
                        bcftools filter -i '{config[VARIANT_FILTER_STRING]}' {params.temp_results_filtered} | bgzip -c > {output.variant_results_file_filtered} ; tabix -fp vcf {output.variant_results_file_filtered}
                    fi
                fi
            fi
        fi
        """

rule variantcalling:
    input:
        path_reference,
        bed = rules.custombed.output.custom_bed if use_gene_list == "true" else (path_bed if BED == "true" and not use_gene_list else (alsgenescanner_bed if alsgenescanner == "true" else [])),
        bam_file = rules.alignment.output.bam_file if config["INPUT_FORMAT"] == "fastq" else (input_dir + "{sample}.bam" if config["INPUT_FORMAT"] == "bam" else (input_dir + "{sample}.cram" if config["INPUT_FORMAT"] == "cram" else (rules.sam2bam.output.bam_file if config["INPUT_FORMAT"] == "sam" else [])))
    output:
        variant_results_file = results_dir + "{sample}/{sample}_sorted.vcf.gz",
        variant_results_file_index = results_dir + "{sample}/{sample}_sorted.vcf.gz.tbi"
    log:
        log_dir + "{sample}/variantcalling.log"
    conda:
        "envs/variantcalling.yaml"
    resources:
        mem_mb = memory
    params:
        variantcalling,
        exome,
        paired,
        BED,
        filter_string,
        out_dir = results_dir,
        sample = "{sample}",
        temp_bed = results_dir + "{sample}/temp.bed.gz",
        sorted_bed = results_dir + "{sample}/sorted.bed.gz"
    shell:
        """
        variantcalling="{params[0]}"
        exome="{params[1]}"
        paired="{params[2]}"
        bed="{params[3]}"
        filter_string="{params[4]}"
        if [[ "$variantcalling" == "True" ]] && [[ "$filter_string" != "True" ]]; then
            {config[STRELKA_DIR]}bin/configureStrelkaGermlineWorkflow.py --bam {input.bam_file} --referenceFasta {input[0]} --runDir {params.out_dir}{params.sample}/strelka
            {params.out_dir}{params.sample}/strelka/runWorkflow.py -j {config[NUMBER_CPU]} -m local
            mv {params.out_dir}{params.sample}/strelka/results/variants/genome.S1.vcf.gz {output.variant_results_file}
            mv {params.out_dir}{params.sample}/strelka/results/variants/genome.S1.vcf.gz.tbi {output.variant_results_file_index}

            if [[ "$exome" == "True" ]]; then
                {config[STRELKA_DIR]}bin/configureStrelkaGermlineWorkflow.py --bam {input.bam_file} --referenceFasta {input[0]} --runDir {params.out_dir}/{params.sample}/strelka --exome
                {params.out_dir}{params.sample}/strelka/runWorkflow.py -j {config[NUMBER_CPU]} -m local
                mv {params.out_dir}{params.sample}/strelka/results/variants/genome.S1.vcf.gz {output.variant_results_file}
                mv {params.out_dir}{params.sample}/strelka/results/variants/genome.S1.vcf.gz.tbi {output.variant_results_file_index}

                if [[ "$bed" == "True" ]]; then
                    bgzip -c {input.bed} > {params.temp_bed}
                    sortBed -i {params.temp_bed} | bgzip -c > {params.sorted_bed}
                    tabix -p bed {params.sorted_bed}
                    {config[STRELKA_DIR]}bin/configureStrelkaGermlineWorkflow.py --bam {input.bam_file} --referenceFasta {input[0]} --runDir {params.out_dir}/{params.sample}/strelka --callRegions {params.sorted_bed}
                    {params.out_dir}{params.sample}/strelka/runWorkflow.py -j {config[NUMBER_CPU]} -m local
                    mv {params.out_dir}{params.sample}/strelka/results/variants/genome.S1.vcf.gz {output.variant_results_file}
                    mv {params.out_dir}{params.sample}/strelka/results/variants/genome.S1.vcf.gz.tbi {output.variant_results_file_index}
                fi
            fi
        fi
        """

rule variantannotation:
    input:
        variant_file = rules.variantcalling.output.variant_results_file if alsgenescanner == "true" or filter_string != "true" else (rules.variantcallingfiltered.output.variant_results_file_filtered if filter_string == "true" else []),
    output:
        annotated_variant_results_file = results_dir + "{sample}/{sample}_SNPindel_annotated.vcf.gz",
        annotated_variant_results_file_index = results_dir + "{sample}/{sample}_SNPindel_annotated.vcf.gz.tbi",
        annotated_variant_results_text = results_dir + "{sample}/{sample}_SNPindel_annotated.txt"
    params:
        annotation,
        variantcalling,
        ref_version = annovar_ref_version,
        out_dir = results_dir,
        sample = "{sample}",
        operations = alsgene_annovar_operations if alsgenescanner == "true" else annovar_operations,
        protocols = alsgene_annovar_protocols if alsgenescanner == "true" else annovar_protocols
    log:
        log_dir + "{sample}/variantannotation.log"
    resources:
        mem_mb = memory
    shell:
        """
        variantcalling="{params[1]}"
        annotation="{params[0]}"
        if [[ "$variantcalling" == "True" ]] && and [[ "$annotation" == "True" ]]; then
            perl {config[ANNOVAR_DIR]}table_annovar.pl --thread {config[NUMBER_CPU]} --vcfinput {input.variant_file} {config[ANNOVAR_DB]} -buildver {params.ref_version} -remove -protocol {params.protocols} -operation {params.operations} -nastring . --outfile {params.out_dir}{params.sample}/{params.sample}_annovar_SNPindel.vcf
            mv {params.out_dir}{params.sample}/{params.sample}_annovar_SNPindel.vcf.{params.ref_version}_multianno.vcf {params.out_dir}{params.sample}/{params.sample}_SNPindel_annotated.vcf
            mv {params.out_dir}{params.sample}/annovar_SNPindel.vcf.{params.ref_version}_multianno.txt {output.annotated_variant_results_text}
            bgzip -f {params.out_dir}{params.sample}/{params.sample}_SNPindel_annotated.vcf ; tabix -fp vcf {params.out_dir}{params.sample}/{params.sample}_SNPindel_annotated.vcf.gz
        fi
        """

rule expansion:
    input:
        path_reference,
        bam_file = input_dir + "{sample}.bam" if config["INPUT_FORMAT"] == "bam" else (input_dir + "{sample}.cram" if config["INPUT_FORMAT"] == "cram" else (rules.sam2bam.output.bam_file if config["INPUT_FORMAT"] == "sam" else (rules.alignment.output.bam_file if config["INPUT_FORMAT"] == "fastq" else []))),
    output:
        expansion_file = results_dir + "{sample}/{sample}_expansions.vcf.gz",
        expansion_file_index = results_dir + "{sample}/{sample}_expansions.vcf.gz.tbi"
    conda:
        "envs/variantcalling.yaml"
    params:
        expansion,
        debug,
        out_dir = results_dir,
        sample = "{sample}",
        variant_catalog = path_expansionHunter_catalog
    log:
        log_dir + "{sample}/expansion.log"
    resources:
        mem_mb = memory
    shell:
        """
        expansion="{params[0]}"
        debug="{params[1]}"
        if [[ "$expansion" == "True" ]]; then
            ExpansionHunter --reads {input.bam_file} --reference {input[0]} --variant-catalog {params.variant_catalog} --output-prefix {params.out_dir}{params.sample}/{params.sample}_expansions
            bgzip {params.out_dir}{params.sample}/{params.sample}_expansions.vcf ; tabix -p vcf {output.expansion_file}
            if [[ "$debug" != "True" ]]; then
                rm {params.out_dir}{params.sample}/{params.sample}_expansions_realigned.bam {params.out_dir}{params.sample}/{params.sample}_expansions.json
            fi
        fi
        """

rule expansionannotation:
    input:
        rules.expansion.output.expansion_file
    output:
        annotated_expansion_file = results_dir + "{sample}/{sample}_expansions_annotated.vcf.gz",
        annotated_expansion_file_index = results_dir + "{sample}/{sample}_expansions_annotated.vcf.gz.tbi"
    conda:
        "envs/annotation.yaml"
    params:
        expansion,
        annotation,
        out_dir = results_dir,
        sample = "{sample}",
        ref_version = annovar_ref_version,
        operations = alsgene_annovar_operations if alsgenescanner == "true" else annovar_operations,
        protocols = alsgene_annovar_protocols if alsgenescanner == "true" else annovar_protocols
    log:
        log_dir + "{sample}/expansionannotation.log"
    resources:
        mem_mb = memory
    shell:
        """
        expansion="{params[0]}"
        annotation="{params[1]}"
        if [[ "$expansion" == "True" ]] && [[ "$annotation" == "True" ]]; then
            perl {config[ANNOVAR_DIR]}table_annovar.pl --thread {config[NUMBER_CPU]} --vcfinput {input[0]} {config[ANNOVAR_DB]} -buildver {params.ref_version} -remove -protocol {params.protocols} -operation {params.operations} -nastring . --outfile {params.out_dir}{params.sample}/annovar_expansions.vcf
            mv {params.out_dir}{params.sample}/annovar_expansions.vcf.{params.ref_version}_multianno.vcf {params.out_dir}{params.sample}/{params.sample}_expansions_annotated.vcf
            bgzip -f {params.out_dir}{params.sample}/{params.sample}_expansions_annotated.vcf ; tabix -fp vcf {params.out_dir}{params.sample}/{params.sample}_expansions_annotated.vcf.gz
        fi
        """

rule STRprofile:
    input:
        path_reference,
        bam_file = input_dir + "{sample}.bam" if config["INPUT_FORMAT"] == "bam" else (input_dir + "{sample}.cram" if config["INPUT_FORMAT"] == "cram" else (rules.sam2bam.output.bam_file if config["INPUT_FORMAT"] == "sam" else (rules.alignment.output.bam_file if config["INPUT_FORMAT"] == "fastq" else []))),
    output:
        STR_profile = results_dir + "{sample}/{sample}_expansiondenovo.str_profile.json",
        genotypeSTR_input = results_dir + "{sample}/{sample}_genotypeSTRinput.txt"
    conda:
        "envs/variantcalling.yaml"
    params:
        STR,
        out_dir = results_dir,
        sample = "{sample}",
    log:
        log_dir + "{sample}/STRprofile.log"
    resources:
        mem_mb = memory
    shell:
        """
        STR="{params[0]}"
        if [[ "$STR" == "True" ]]; then
            {config[EXPANSIONHUNTERDENOVO_DIR]}bin/ExpansionHunterDenovo profile --reads {input.bam_file} --reference {input[0]} --output-prefix {params.out_dir}{params.sample}/{params.sample}_expansiondenovo --min-anchor-mapq 50 --max-irr-mapq 40 --log-reads
            cat {params.out_dir}{params.sample}/{params.sample}_expansiondenovo.locus.tsv | sed 's/contig/chr/g' | cut -f1-4 > {output.genotypeSTR_input}
        fi
        """

rule genotypeSTR:
    input:
        rules.STRprofile.output.genotypeSTR_input,
        path_reference,
        bam_file = input_dir + "{sample}.bam" if config["INPUT_FORMAT"] == "bam" else (input_dir + "{sample}.cram" if config["INPUT_FORMAT"] == "cram" else (rules.sam2bam.output.bam_file if config["INPUT_FORMAT"] == "sam" else (rules.alignment.output.bam_file if config["INPUT_FORMAT"] == "fastq" else []))),
    output:
        EHDN_variant_catalog = results_dir + "{sample}/{sample}_EHDN_variant_catalog.json",
        EHDN_expansion_file = results_dir + "{sample}/{sample}_EHDNexpansions.vcf.gz",
        EHDN_expansion_file_index = results_dir + "{sample}/{sample}_EHDNexpansions.vcf.gz.tbi"
    conda:
        "envs/variantcalling.yaml"
    params:
        STR,
        genotypeSTR,
        out_dir = results_dir,
        sample = "{sample}",
        EHDN_excluded = results_dir + "{sample}/{sample}_EHDN_excluded.csv",
        EHDN_unmatched = results_dir + "{sample}/{sample}_EHDN_unmatched.csv"
    log:
        log_dir + "{sample}/genotypeSTR.log"
    resources:
        mem_mb = memory
    shell:
        """
        STR="{params[0]}"
        genotype="{params[1]}"
        if [[ "$STR" == "True" ]] && [[ "$genotype" == "True" ]]; then
            python scripts/conversion_EHDN_catalog.py {input[0]} {input[1]} {output.EHDN_variant_catalog} {params.EHDN_unmatched} {params.EHDN_excluded}
            ExpansionHunter --reads {input[2]} --reference {input[1]} --variant-catalog {output.EHDN_variant_catalog} --output-prefix {params.out_dir}{params.sample}_EHDNexpansions
            bgzip {params.out_dir}{params.sample}_EHDNexpansions.vcf
            tabix -p vcf {output.EHDN_expansion_file}
        fi
        """

rule STRannotation:
    input:
        rules.genotypeSTR.output.EHDN_expansion_file
    output:
        annotated_EHDN_expansion_file = results_dir + "{sample}/{sample}_STR_annotated.vcf.gz",
        annotated_EHDN_expansion_file_index = results_dir + "{sample}/{sample}_STR_annotated.vcf.gz.tbi"
    conda:
        "envs/annotation.yaml"
    params:
        STR,
        genotypeSTR,
        annotation,
        out_dir = results_dir,
        sample = "{sample}",
        ref_version = annovar_ref_version,
        protocols = alsgene_annovar_protocols if alsgenescanner == "true" else annovar_protocols,
        operations = alsgene_annovar_operations if alsgenescanner == "true" else annovar_operations
    log:
        log_dir + "{sample}/STRannotation.log"
    resources:
        mem_mb = memory
    shell:
        """
        STR="{params[0]}"
        genotype="{params[1]}"
        annotation="{params[2]}"
        if [[ "$STR" == "True" ]] && [[ "$genotype" == "True" ]] && [[ "$annotation" == "True" ]]; then
            perl {config[ANNOVAR_DIR]}table_annovar.pl --thread {config[NUMBER_CPU]} --vcfinput {input[0]} {config[ANNOVAR_DB]} -buildver {params.ref_version} -remove -protocol {params.protocols} -operation {params.operations} -nastring . --outfile {params.out_dir}{params.sample}/annovar_EHDNexpansions.vcf
            mv {params.out_dir}{params.sample}/annovar_EHDNexpansions.vcf.{params.ref_version}_multianno.vcf {params.out_dir}{params.sample}/{params.sample}_STR_annotated.vcf
            bgzip -f {params.out_dir}{params.sample}/{params.sample}_STR_annotated.vcf ; tabix -fp vcf {output.annotated_EHDN_expansion_file}
        fi
        """

rule CramToBam:
    input:
        path_reference,
        input_file = input_dir + "{sample}.cram"
    output:
        delly_bam = results_dir + "{sample}/{sample}_delly.bam",
        delly_bam_index = results_dir + "{sample}/{sample}_delly.bam.bai"
    params:
        SV,
        format,
        out_dir = results_dir,
        sample = "{sample}"
    conda:
        "envs/samtools.yaml"
    log:
        log_dir + "{sample}/cramtobam.log"
    resources:
        mem_mb = memory
    shell:
        """
        SV="{params[0]}"
        format="{params[1]}"
        if [[ "$SV" == "True" ]] && [[ "$format" == "cram" ]]; then
            samtools view -b -h -@ {config[NUMBER_CPU]} -T {input[0]} -o {output.delly_bam} {input.input_file}
            samtools index -@ {config[NUMBER_CPU]} {output.delly_bam}
        fi
        """

rule SV:
    input:
        path_reference,
        bed = rules.custombed.output.custom_bed if use_gene_list == "true" else (path_bed if BED == "true" and not use_gene_list else (alsgenescanner_bed if alsgenescanner == "true" else [])),
        bam_file = input_dir + "{sample}.bam" if config["INPUT_FORMAT"] == "bam" else (rules.CramToBam.output.delly_bam if config["INPUT_FORMAT"] == "cram" else (rules.sam2bam.output.bam_file if config["INPUT_FORMAT"] == "sam" else (rules.alignment.output.bam_file if config["INPUT_FORMAT"] == "fastq" else [])))
    output:
        manta_SV = results_dir + "{sample}/{sample}_SV_manta.vcf",
        delly_SV = results_dir + "{sample}/{sample}_delly_SV.vcf",
        merged_SV = results_dir + "{sample}/{sample}_merged_SV.vcf.gz",
        merged_SV_index = results_dir + "{sample}/{sample}_merged_SV.vcf.gz.tbi"
    conda:
        "envs/variantcalling.yaml"
    params:
        paired,
        BED,
        SV,
        temp_bed = results_dir + "{sample}/temp.bed.gz",
        sorted_bed = results_dir + "{sample}/sorted.bed.gz",
        out_dir = results_dir,
        sample = "{sample}",
        delly_exclude_regions = path_delly_exclude_regions,
    log:
        log_dir + "{sample}/SV.log"
    resources:
        mem_mb = memory
    shell:
        """
        paired="{params[0]}"
        use_bed="{params[1]}"
        SV="{params[2]}"
        if [[ "$SV" == "True" ]] && [[ "$paired" == "paired" ]]; then
            if [[ "$use_bed" == "True" ]]; then
                bgzip -c {input.bed} > {params.temp_bed}
                sortBed -i {params.temp_bed} | bgzip -c > {params.sorted_bed}
                tabix -p bed {params.sorted_bed}
                {config[MANTA_DIR]}bin/configManta.py --bam {input.bam_file} --referenceFasta {input[0]} --runDir {params.out_dir}{params.sample}/SV_manta --callRegions {params.sorted_bed}
                {params.out_dir}{params.sample}/SV_manta/runWorkflow.py -j {config[NUMBER_CPU]} -m local
                mv {params.out_dir}{params.sample}/SV_manta/results/variants/diploidSV.vcf.gz {params.out_dir}{params.sample}/{params.sample}_SV_manta.vcf.gz
                gzip -d {params.out_dir}{params.sample}/{params.sample}_SV_manta.vcf.gz
                rm -r {params.out_dir}{params.sample}/SV_manta
                delly call -g {input[0]} -o {params.out_dir}{params.sample}/{params.sample}_delly_SV.bcf -x {params.delly_exclude_regions} {input.bam_file}
                bcftools view {params.out_dir}{params.sample}/{params.sample}_delly_SV.bcf > {output.delly_SV}
                ls {params.out_dir}{params.sample}/*SV.vcf > {params.out_dir}{params.sample}/survivor_sample_files
                SURVIVOR merge {params.out_dir}{params.sample}/survivor_sample_files 1000 1 1 1 0 30 {params.out_dir}{params.sample}/{params.sample}_merged_SV.vcf
                perl scripts/vcf-sort.pl {params.out_dir}{params.sample}/{params.sample}_merged_SV.vcf | bgzip -c > {output.merged_SV}
                tabix -p vcf {output.merged_SV}
                rm {params.out_dir}{params.sample}/survivor_sample_files {params.out_dir}{params.sample}/{params.sample}_delly_SV.bcf {params.out_dir}{params.sample}/{params.sample}_delly_SV.bcf.csi {params.out_dir}{params.sample}/{params.sample}_merged_SV.vcf

            else
                {config[MANTA_DIR]}bin/configManta.py --bam {input.bam_file} --referenceFasta {input[0]} --runDir {params.out_dir}{params.sample}/SV_manta
                {params.out_dir}{params.sample}/SV_manta/runWorkflow.py -j {config[NUMBER_CPU]} -m local
                mv {params.out_dir}{params.sample}/SV_manta/results/variants/diploidSV.vcf.gz {params.out_dir}{params.sample}/{params.sample}_SV_manta.vcf.gz
                gzip -d {params.out_dir}{params.sample}/{params.sample}_SV_manta.vcf.gz
                rm -r {params.out_dir}{params.sample}/SV_manta
                delly call -g {input[0]} -o {params.out_dir}{params.sample}/{params.sample}_delly_SV.bcf -x {params.delly_exclude_regions} {input.bam_file}
                bcftools view {params.out_dir}{params.sample}/{params.sample}_delly_SV.bcf > {output.delly_SV}
                ls {params.out_dir}{params.sample}/*SV.vcf > {params.out_dir}{params.sample}/survivor_sample_files
                SURVIVOR merge {params.out_dir}{params.sample}/survivor_sample_files 1000 1 1 1 0 30 {params.out_dir}{params.sample}/{params.sample}_merged_SV.vcf
                perl scripts/vcf-sort.pl {params.out_dir}{params.sample}/{params.sample}_merged_SV.vcf | bgzip -c > {output.merged_SV}
                tabix -p vcf {output.merged_SV}
                rm {params.out_dir}{params.sample}/survivor_sample_files {params.out_dir}{params.sample}/{params.sample}_delly_SV.bcf {params.out_dir}{params.sample}/{params.sample}_delly_SV.bcf.csi {params.out_dir}{params.sample}/{params.sample}_merged_SV.vcf
            fi
        fi
        """

rule SVannotation:
    input:
        rules.SV.output.merged_SV,
    output:
        SV_annotation_file = results_dir + "{sample}/{sample}_SV_annotated.tsv"
    conda:
        "envs/SVannotation.yaml"
    params:
        alsgenescanner,
        annotation,
        SV,
        out_dir = results_dir,
        sample = "{sample}",
        genomebuild = annotsv_ref_version,
    log:
        log_dir + "{sample}/SVannotation.log"
    resources:
        mem_mb = memory
    shell:
        """
        alsgenescanner="{params[0]}"
        annotation="{params[1]}"
        SV="{params[2]}"
        if [[ "$SV" == "True" ]] && [[ "$annotation" == "True" ]]; then
            if [[ "$alsgenescanner" == "True" ]]; then
                cpan YAML::XS
                cpan Sort::Key::Natural
                export ANNOTSV={config[ANNOTSV_DIR]}
                {config[ANNOTSV_DIR]}bin/AnnotSV -annotationsDir {config[ANNOTSV_DIR]}/share/AnnotSV -SVinputFile {input[0]} -genomeBuild {params.genomebuild} -candidateGenesFiltering yes -candidateGenesFile resources/alsgenescanner/list_genes_all.txt outputDir {params.out_dir}{params.sample} -outputFile {params.sample}_SV_annotated -SVminSize 30 {config[ANNOTSV_CUSTOM_OPTIONS]}
            else
                cpan YAML::XS
                cpan Sort::Key::Natural
                export ANNOTSV={config[ANNOTSV_DIR]}
                {config[ANNOTSV_DIR]}bin/AnnotSV -annotationsDir {config[ANNOTSV_DIR]}/share/AnnotSV -SVinputFile {input[0]} -genomeBuild {params.genomebuild} outputDir {params.out_dir}{params.sample} -outputFile {params.sample}_SV_annotated -SVminSize 30 {config[ANNOTSV_CUSTOM_OPTIONS]}
            fi
        fi
        """

rule MEI:
    input:
        path_reference,
        melt_bed,
        bam_file = input_dir + "{sample}.bam" if config["INPUT_FORMAT"] == "bam" else (input_dir + "{sample}.cram" if config["INPUT_FORMAT"] == "cram" else (rules.sam2bam.output.bam_file if config["INPUT_FORMAT"] == "sam" else (rules.alignment.output.bam_file if config["INPUT_FORMAT"] == "fastq" else [])))
    output:
        MEI_file = results_dir + "{sample}/{sample}_MEI.vcf.gz",
        MEI_file_index = results_dir + "{sample}/{sample}_MEI.vcf.gz.tbi"
    conda:
        "envs/MEI.yaml"
    params:
        exome,
        MEI,
        out_dir = results_dir,
        sample="{sample}",
        zipped = melt_zipped_files,
        transposon_list = path_melt + "transposon_reference.list",
        removal_dir=input_dir + "{sample}.bam" if config["INPUT_FORMAT"] == "bam" or config["INPUT_FORMAT"] == "cram" else (rules.sam2bam.output.bam_file if config["INPUT_FORMAT"] == "sam" else (rules.alignment.output.bam_file if config["INPUT_FORMAT"] == "fastq" else []))
    log:
        log_dir + "{sample}/MEI.log"
    resources:
        mem_mb = memory
    shell:
        """
        exome="{params[0]}"
        MEI="{params[1]}"
        if [[ "$MEI" == "True" ]]; then
            if [[ "$exome" == "True" ]]; then
                mkdir {params.out_dir}/{params.sample}/melt
                ls {params[1]} | sed 's/\*//g' > {params.transposon_list}
                java -Xmx{config[MEM_GB]}G -jar {config[MELT_DIR]}MELT.jar Single -bamfile {input.bam_file} -h {input[0]} -t {params.transposon_list} -n {input[1]} -w {params.out_dir}/{params.sample}/melt -exome {config[MELT_CUSTOM_OPTIONS]}
                cat {params.out_dir}/{params.sample}/melt/SVA.final_comp.vcf | grep '^#' > {params.out_dir}/{params.sample}/melt/{params.sample}.header.txt
                cat {params.out_dir}/{params.sample}/melt/SVA.final_comp.vcf | grep -v '^#' > {params.out_dir}/{params.sample}/melt/{params.sample}.sva.vcf
                cat {params.out_dir}/{params.sample}/melt/LINE1.final_comp.vcf | grep -v '^#' > {params.out_dir}/{params.sample}/melt/{params.sample}.line1.vcf
                cat {params.out_dir}/{params.sample}/melt/ALU.final_comp.vcf | grep -v '^#' > {params.out_dir}/{params.sample}/melt/{params.sample}.alu.vcf
                cat {params.out_dir}/{params.sample}/melt/{params.sample}.header.txt {params.out_dir}/{params.sample}/melt/{params.sample}.sva.vcf {params.out_dir}/{params.sample}/melt/{params.sample}.line1.vcf {params.out_dir}/{params.sample}/melt/{params.sample}.alu.vcf | perl scripts/vcf-sort.pl -c | bgzip -c > {output.MEI_file}
                tabix -p vcf {output.MEI_file}
                rm -r {params.out_dir}/{params.sample}/melt
                rm {params.removal_dir}.disc {params.removal_dir}.disc.bai {params.removal_dir}.fq
            else
                ls {params[1]} | sed 's/\*//g' > {params.transposon_list}
                java -Xmx{config[MEM_GB]}G -jar {config[MELT_DIR]}MELT.jar Single -bamfile {input.bam_file} -h {input[0]} -t {params.transposon_list} -n {input[1]} -w {params.out_dir}/{params.sample}/melt {config[MELT_CUSTOM_OPTIONS]}
                cat {params.out_dir}/{params.sample}/melt/SVA.final_comp.vcf | grep '^#' > {params.out_dir}/{params.sample}/melt/{params.sample}.header.txt
                cat {params.out_dir}/{params.sample}/melt/SVA.final_comp.vcf | grep -v '^#' > {params.out_dir}/{params.sample}/melt/{params.sample}.sva.vcf
                cat {params.out_dir}/{params.sample}/melt/LINE1.final_comp.vcf | grep -v '^#' > {params.out_dir}/{params.sample}/melt/{params.sample}.line1.vcf
                cat {params.out_dir}/{params.sample}/melt/ALU.final_comp.vcf | grep -v '^#' > {params.out_dir}/{params.sample}/melt/{params.sample}.alu.vcf
                cat {params.out_dir}/{params.sample}/melt/{params.sample}.header.txt {params.out_dir}/{params.sample}/melt/{params.sample}.sva.vcf {params.out_dir}/{params.sample}/melt/{params.sample}.line1.vcf {params.out_dir}/{params.sample}/melt/{params.sample}.alu.vcf | perl scripts/vcf-sort.pl -c | bgzip -c > {output.MEI_file}
                tabix -p vcf {output.MEI_file}
                rm -r {params.out_dir}/{params.sample}/melt
                rm {params.removal_dir}.disc {params.removal_dir}.disc.bai {params.removal_dir}.fq
            fi
        fi
        """

rule MEIannotation:
    input:
        rules.MEI.output.MEI_file,
    output:
        MEI_annotation_file = results_dir + "{sample}/{sample}_MEI_annotated.tsv"
    conda:
        "envs/SVannotation.yaml"
    params:
        alsgenescanner,
        MEI,
        annotation,
        out_dir = results_dir,
        sample = "{sample}",
        genomebuild = annotsv_ref_version,
    log:
        log_dir + "{sample}/MEIannotation.log"
    resources:
        mem_mb = memory
    shell:
        """
        alsgenescanner="{params[0]}"
        MEI="{params[1]}"
        annotation="{params[2]}"
        if [[ "$MEI" == "True" ]] && [[ "$annotation" == "True" ]]; then
            if [[ "$alsgenescanner" == "True" ]]; then
                cpan YAML::XS
                cpan Sort::Key::Natural
                export ANNOTSV={config[ANNOTSV_DIR]}
                {config[ANNOTSV_DIR]}bin/AnnotSV -annotationsDir {config[ANNOTSV_DIR]}/share/AnnotSV -SVinputFile {input[0]} -genomeBuild {params.genomebuild} -candidateGenesFiltering yes -candidateGenesFile resources/alsgenescanner/list_genes_all.txt outputDir {params.out_dir}{params.sample} -outputFile {params.sample}_MEI_annotated -SVminSize 30 {config[ANNOTSV_CUSTOM_OPTIONS]}
            else
                cpan YAML::XS
                cpan Sort::Key::Natural
                export ANNOTSV={config[ANNOTSV_DIR]}
                {config[ANNOTSV_DIR]}bin/AnnotSV -annotationsDir {config[ANNOTSV_DIR]}/share/AnnotSV -SVinputFile {input[0]} -genomeBuild {params.genomebuild} outputDir {params.out_dir}{params.sample} -outputFile {params.sample}_MEI_annotated -SVminSize 30 {config[ANNOTSV_CUSTOM_OPTIONS]}
            fi
        fi
        """

rule SVandMEImergingandannotation:
    input:
        rules.SV.output.merged_SV,
        rules.MEI.output.MEI_file
    output:
        merged_SV_MEI = results_dir + "{sample}/{sample}_SV_MEI.merged.vcf.gz",
        merged_SV_MEI_index = results_dir + "{sample}/{sample}_SV_MEI.merged.vcf.gz.tbi",
        SV_MEI_annotation_file = results_dir + "{sample}/{sample}_SV_MEI_annotated.tsv"
    conda:
        "envs/variantcalling.yaml"
    params:
        SV,
        MEI,
        annotation,
        alsgenescanner,
        merged_dir = results_dir + "{sample}/merging",
        out_dir = results_dir,
        sample = "{sample}",
        SV_vcf = results_dir + "{sample}/merging/{sample}_merged_SV.vcf",
        MEI_vcf = results_dir + "{sample}/merging/{sample}_MEI.vcf",
        genomebuild = annotsv_ref_version
    log:
        log_dir + "{sample}/SVandMEImerging.log"
    resources:
        mem_mb = memory
    shell:
        """
        SV="{params[0]}"
        MEI="{params[1]}"
        annotation="{params[2]}"
        alsgenescanner="{params[3]}"
        if [[ "$SV" == "True" ]] && [[ "$MEI" == "True" ]]; then
            mkdir {params.merged_dir}
            bgzip -d {input[0]} > {params.SV_vcf}
            bgzip -d {input[1]} > {params.MEI_vcf}
            ls {params.merged_dir}*.vcf > {params.out_dir}{params.sample}survivor_sample_files
            SURVIVOR merge {params.out_dir}{params.sample}survivor_sample_files 1000 1 1 1 0 30 {params.out_dir}{params.sample}/{params.sample}_SV_MEI.merged.vcf
            perl scripts/vcf-sort.pl {params.out_dir}{params.sample}/{params.sample}_SV_MEI.merged.vcf | bgzip -c > {output.merged_SV_MEI}
            tabix -p vcf {output.merged_SV_MEI}
            rm {params.out_dir}{params.sample}/survivor_sample_files

            if [[ "$annotation" == "True" ]]; then
                if [[ "$alsgenescanner" == "True" ]]; then
                    cpan YAML::XS
                    cpan Sort::Key::Natural
                    export ANNOTSV={config[ANNOTSV_DIR]}
                    {config[ANNOTSV_DIR]}bin/AnnotSV -annotationsDir {config[ANNOTSV_DIR]}/share/AnnotSV -SVinputFile {input[0]} -genomeBuild {params.genomebuild} -candidateGenesFiltering yes -candidateGenesFile resources/alsgenescanner/list_genes_all.txt outputDir {params.out_dir}{params.sample} -outputFile {params.sample}_SV_MEI_annotated -SVminSize 30 {custom[ANNOTSV_CUSTOM_OPTIONS]}

                else
                    cpan YAML::XS
                    cpan Sort::Key::Natural
                    export ANNOTSV={config[ANNOTSV_DIR]}
                    {config[ANNOTSV_DIR]}bin/AnnotSV -annotationsDir {config[ANNOTSV_DIR]}/share/AnnotSV -SVinputFile {input[0]} -genomeBuild {params.genomebuild} outputDir {params.out_dir}{params.sample} -outputFile {params.sample}_SV_MEI_annotated -SVminSize 30 {config[ANNOTSV_CUSTOM_OPTIONS]}
                fi
            fi
        fi
        """

rule extractnonhumanreads:
    input:
        bam_file = input_dir + "{sample}.bam" if config["INPUT_FORMAT"] == "bam" else (input_dir + "{sample}.cram" if config["INPUT_FORMAT"] == "cram" else (rules.sam2bam.output.bam_file if config["INPUT_FORMAT"] == "sam" else (rules.alignment.output.bam_file if config["INPUT_FORMAT"] == "fastq" else []))),
    output:
        unaligned_fastq = results_dir + "{sample}/unaligned_reads.fastq.gz"
    conda:
        "envs/samtools.yaml"
    params:
        virus,
        bacteria,
        microbes,
        out_dir = results_dir,
        sample = "{sample}",
    log:
        log_dir + "{sample}/extractnonhumanreads.log"
    resources:
        mem_mb = memory
    shell:
        """
        virus="{params[0]}"
        bacteria="{params[1]}"
        microbes="{params[2]}"
        if [[ "$virus" == "True" || "$bacteria" == "True" || "$microbes" == "True" ]]; then
            samtools view -@ {config[NUMBER_CPU]} -hf 4 {input[0]} | samtools bam2fq -s {params.out_dir}{params.sample}/singleton_reads.fastq -@ {config[NUMBER_CPU]} - > {params.out_dir}{params.sample}/unaligned_reads.fastq
            cat {params.out_dir}{params.sample}/singleton_reads.fastq >> {params.out_dir}{params.sample}/unaligned_reads.fastq ; gzip {params.out_dir}{params.sample}/unaligned_reads.fastq
        fi
        """

rule identifyviralmaterial:
    input:
        rules.extractnonhumanreads.output.unaligned_fastq,
    output:
        virus_stats = results_dir + "{sample}/{sample}_virus_stats.txt",
        virus_report = reports_dir + "{sample}/{sample}_virus_report.txt"
    conda:
        "envs/alignmentfast.yaml"
    params:
        virus,
        out_dir = results_dir,
        sample = "{sample}"
    log:
        log_dir + "{sample}/identifyvirus.log"
    resources:
        mem_mb = memory
    shell:
        """
        virus="{params[0]"
        if [[ "$virus" == "True" ]]; then
            hisat2 --no-spliced-alignment -p {config[NUMBER_CPU]} -x {config[VIRUS_INDEX]} -U {input[0]} | samtools view -@ {config[NUMBER_CPU]} -hSb - | samtools sort -@ {config[NUMBER_CPU]} -T {params.out_dir}{params.sample}tempvirus.file -o {params.out_dir}{params.sample}/{params.sample}_output_virus.bam -
            samtools index -@ {config[NUMBER_CPU]} {params.out_dir}{params.sample}/{params.sample}_output_virus.bam ; samtools idxstats {params.out_dir}{params.sample}/{params.sample}_output_virus.bam > {output.virus_stats}
            python scripts/non_human_reads_report.py {output.virus_stats} {params.out_dir}{params.sample}/{params.sample}_virus_stats.list {params.out_dir}{params.sample}/{params.sample}_output_virus.bam {params.out_dir}{params.sample}/{params.sample}_virus_coverage_stats.txt {output.virus_report}
        fi
        """

rule identifybacterialmaterial:
    input:
        rules.extractnonhumanreads.output.unaligned_fastq,
    output:
        bacteria_stats = results_dir + "{sample}/{sample}_bacteria_stats.txt",
        bacteria_report = reports_dir + "{sample}/{sample}_bacteria_report.txt"
    conda:
        "envs/alignmentfast.yaml"
    params:
        bacteria,
        out_dir = results_dir,
        sample = "{sample}"
    log:
        log_dir + "{sample}/identifybacteria.log"
    resources:
        mem_mb = memory
    shell:
        """
        bacteria="{params[0]}"
        if [[ "$bacteria" == "True" ]]; then
            hisat2 --no-spliced-alignment -p {config[NUMBER_CPU]} -x {config[BACTERIA_INDEX]} -U {input[0]} | samtools view -@ {config[NUMBER_CPU]} -hSb - | samtools sort -@ {config[NUMBER_CPU]} -T {params.out_dir}{params.sample}tempbacteria.file -o {params.out_dir}{params.sample}/{params.sample}_output_bacteria.bam -
            samtools index -@ {config[NUMBER_CPU]} {params.out_dir}{params.sample}/{params.sample}_output_bacteria.bam ; samtools idxstats {params.out_dir}{params.sample}/{params.sample}_output_bacteria.bam > {output.bacteria_stats}
            python scripts/non_human_reads_report.py {output.bacteria_stats} {params.out_dir}{params.sample}/{params.sample}_bacteria_stats.list {params.out_dir}{params.sample}/{params.sample}_output_bacteria.bam {params.out_dir}{params.sample}/{params.sample}_bacteria_coverage_stats.txt {output.bacteria_report}
        fi
        """

rule identifycustommaterial:
    input:
        rules.extractnonhumanreads.output.unaligned_fastq,
    output:
        microbes_stats = results_dir + "{sample}/{sample}_microbes_stats.txt",
        microbes_report = reports_dir + "{sample}/{sample}_microbes_report.txt"
    conda:
        "envs/alignmentfast.yaml"
    params:
        microbes,
        out_dir = results_dir,
        sample = "{sample}"
    log:
        log_dir + "{sample}/identifycustommicrobes.log"
    resources:
        mem_mb = memory
    shell:
        """
        microbes="{params[0]}"
        if [[ "$microbes" == "True" ]]; then
            hisat2 --no-spliced-alignment -p {config[NUMBER_CPU]} -x {config[CUSTOM_MICROBES_INDEX]} -U {input[0]} | samtools view -@ {config[NUMBER_CPU]} -hSb - | samtools sort -@ {config[NUMBER_CPU]} -T {params.out_dir}{params.sample}tempmicrobes.file -o {params.out_dir}{params.sample}/{params.sample}_output_microbes.bam -
            samtools index -@ {config[NUMBER_CPU]} {params.out_dir}{params.sample}/{params.sample}_output_microbes.bam ; samtools idxstats {params.out_dir}{params.sample}/{params.sample}_output_microbes.bam > {output.microbes_stats}
            python scripts/non_human_reads_report.py {output.micrboes_stats} {params.out_dir}{params.sample}/{params.sample}_microbes_stats.list {params.out_dir}{params.sample}/{params.sample}_output_microbes.bam {params.out_dir}{params.sample}/{params.sample}_microbes_coverage_stats.txt {output.microbes_report}
        fi
        """

rule alignmentreport:
    input:
        bam_file = input_dir + "{sample}.bam" if config["INPUT_FORMAT"] == "bam" else (input_dir + "{sample}.cram" if config["INPUT_FORMAT"] == "cram" else (rules.sam2bam.output.bam_file if config["INPUT_FORMAT"] == "sam" else (rules.alignment.output.bam_file if config["INPUT_FORMAT"] == "fastq" else []))),
    output:
        flagstat = reports_dir + "{sample}/{sample}_alignment_flagstat.txt",
        stats = reports_dir + "{sample}/{sample}_alignment_stats.txt"
    conda:
        "envs/simplereports.yaml"
    params:
        alignmentreport = config["ALIGNMENT_REPORT"]
    log:
        log_dir + "{sample}/alignmentreport.log"
    resources:
        mem_mb = memory
    shell:
        """
        alignment_report="{params.alignmentreport}"
        if [[ "$alignment_report" == "True" ]]; then
            samtools flagstat -@ {config[NUMBER_CPU]} {input[0]} > {output.flagstat}
            samtools stats -@ {config[NUMBER_CPU]} {input[0]} > {output.stats}
        fi
        """

rule sequencingreport:
    input:
        input_dir + "{sample}.1.fq.gz" if format == "fastq" else [],
        input_dir + "{sample}.2.fq.gz" if format == "fastq" and paired == "paired" else []
    output:
        fastqc_report = reports_dir + "{sample}/{sample}_sequencing_report.txt"
    conda:
        "envs/simplereports.yaml"
    params:
        format,
        out_dir = reports_dir,
        sequencingreport = config["SEQUENCING_REPORT"],
    log:
        log_dir + "{sample}sequencingreport.log"
    resources:
        mem_mb = memory
    shell:
        """
        format="{params[0]}"
        sequencing_report="{params.sequencingreport}"
        if [[ "$sequencing_report" == "True" ]] && [[ "$format" == "fastq" ]]; then
            fastqc -o {params.out_dir} -f fastq -t {config[NUMBER_CPU]} {input[0]} {input[1]}
        fi
        """

rule callsreport:
    input:
        variant_file = rules.variantcallingfiltered.output.variant_results_file_filtered if variantcalling == "true" and filter_string == "true" else (rules.variantcalling.output.variant_results_file if variantcalling == "true" and filter_string != "true" else [])
    output:
        vcfstats = reports_dir + "{sample}/{sample}_calls_vcfstats.txt"
    conda:
        "envs/simplereports.yaml"
    params:
        variantcalling,
        out_dir = reports_dir,
        sample = "{sample}",
        callsreport = config["SNP_INDEL_CALLS_REPORT"],
    log:
        log_dir + "{sample}/callsreport.log"
    resources:
        mem_mb = memory
    shell:
        """
        calls_report="{params.callsreport}"
        variantcalling="{params[0]}"
        if [[ "$calls_report" == "True" ]] && and [[ "$variantcalling" == "True" ]]; then
            bcftools stats --threads {config[NUMBER_CPU]} {input.variant_file} > {output.vcfstats}
        fi
        """

rule multireport:
    input:
        reports_dir
    conda:
        "envs/simplereports.yaml"
    output:
        reports_dir + "{sample}_multiqc_report.html"
    params:
        alignmentreport = config["ALIGNMENT_REPORT"],
        sequencingreport = config["SEQUENCING_REPORT"],
        callsreport = config["SNP_INDEL_CALLS_REPORT"]
    log:
        log_dir + "{sample}/multireport.log"
    resources:
        mem_mb = memory
    shell:
        """
        alignment_report="{params.alignmentreport}"
        sequencing_report="{params.sequencingreport}"
        calls_report="{params.callsreport}"
        if [[ "$alignment_report" == "True" || "$sequencing_report" == "True" || "$calls_report" == "True" ]]; then
            multiqc -o {input[0]} {input[0]}
        fi
        """

rule variantresultsreport:
    input:
        rules.variantannotation.output.annotated_variant_results_file,
        path_gene_list,
        rules.variantannotation.output.annotated_variant_results_text if alsgenescanner == "true" else []
    output:
        SNPindel_annotation_report = reports_dir + "{sample}/{sample}_annovar_SNPindel.txt",
        variant_annotation_file_zipped = results_dir + "{sample}/{sample}_SNPindel_annotated.vcf.gz",
        variant_annotation_file_zipped_index = results_dir + "{sample}/{sample}_SNPindel_annotated.vcf.gz.tbi"
    conda:
        "envs/AGS.yaml"
    params:
        annotation,
        variantcalling,
        out_dir = reports_dir,
        results = results_dir,
        sample = "{sample}",
        resultsreport = config["ANNOTATION_RESULTS_REPORT"]
    log:
        log_dir + "{sample}/variantannotationreport.log"
    resources:
        mem_mb = memory
    shell:
        """
        annotation="{params[0]}"
        results_report="{params.resultsreport}"
        variantcalling="{params[1]}"
        if [[ "$variantcalling" == "True" ]] && [[ "$annotation" == "True" ]] && [[ "$results_report" == "True" ]]; then
            gzip -d {input[0]}
            python scripts/simpleannotationreport.py {params.results}{params.sample}/{params.sample}_SNPindel_annotated.vcf {input[1]} {output.SNPindel_annotation_report}
            bgzip {params.results}{params.sample}/{params.sample}_SNPindel_annotated.vcf
            tabix -p vcf {output.variant_annotation_file_zipped}
        fi
        """

rule expansionresultsreport:
    input:
        rules.expansionannotation.output.annotated_expansion_file,
        path_gene_list
    output:
        expansion_annotation_report = reports_dir + "{sample}/{sample}_annovar_expansions.txt",
        expansion_annotation_file_zipped = results_dir + "{sample}/{sample}_expansions_annotated.vcf.gz",
        expansion_annotation_file_zipped_index = results_dir + "{sample}/{sample}_expansions_annotated.vcf.gz.tbi"
    conda:
        "envs/AGS.yaml"
    params:
        annotation,
        expansion,
        out_dir = reports_dir,
        results = results_dir,
        sample = "{sample}",
        resultsreport = config["ANNOTATION_RESULTS_REPORT"]
    log:
        log_dir + "{sample}/expansionannotationreport.log"
    resources:
        mem_mb = memory
    shell:
        """
        annotation="{params[0]}"
        results_report="{params.resultsreport}"
        expansion="{params[1]}"
        if [[ "$expansion" == "True" ]] && [[ "$annotation" == "True" ]] && [[ "$results_report" == "True" ]]; then
            gzip -d {input[0]}
            python scripts/simpleannotationreport.py {params.results}{params.sample}/{params.sample}_expansions_annotated.vcf {input[1]} {output.expansion_annotation_report}
            bgzip {params.results}{params.sample}/{params.sample}_expansions_annotated.vcf
            tabix -p vcf {output.expansion_annotation_file_zipped}
        fi
        """

rule STRresultsreport:
    input:
        rules.STRannotation.output.annotated_EHDN_expansion_file,
        path_gene_list
    output:
        STR_annotation_report = reports_dir + "{sample}/{sample}_annovar_STR.txt",
        STR_annotation_file_zipped = results_dir + "{sample}/{sample}_EHDNexpansions_annotated.vcf.gz",
        STR_annotation_file_zipped_index = results_dir + "{sample}/{sample}_EHDNexpansions_annotated.vcf.gz.tbi"
    conda:
        "envs/AGS.yaml"
    params:
        annotation,
        genotypeSTR,
        out_dir = reports_dir,
        results = results_dir,
        sample = "{sample}",
        resultsreport = config["ANNOTATION_RESULTS_REPORT"],
    log:
        log_dir + "{sample}/STRannotationreport.log"
    resources:
        mem_mb = memory
    shell:
        """
        annotation="{params[0]}"
        results_report="{params.resultsreport}"
        genotype="{params[1]}"
        if [[ "$genotype" == "True" ]] && [[ "$annotation" == "True" ]] && [[ "$results_report" == "True" ]]; then
            gzip -d {input[0]}
            python scripts/simpleannotationreport.py {params.results}{params.sample}/{params.sample}_STR_annotated.vcf {input[1]} {output.STR_annotation_report}
            bgzip {params.results}{params.sample}/{params.sample}_STR_annotated.vcf
            tabix -p vcf {output.STR_annotation_file_zipped}
        fi
        """

rule SVandMEIreport:
    input:
        rules.SVandMEImergingandannotation.output.SV_MEI_annotation_file
    output:
        SV_MEI_report = reports_dir + "{sample}/{sample}_SV_MEI_annotated.html"
    conda:
        "envs/SVannotation.yaml"
    params:
        SV,
        MEI,
        annotation,
        out_dir = reports_dir,
        sample = "{sample}",
        ref_version = annovar_ref_version,
        resultsreport = config["ANNOTATION_RESULTS_REPORT"],
    log:
        log_dir + "{sample}/SVMEIreport.log"
    resources:
        mem_mb = memory
    shell:
        """
        SV="{params[0]}"
        MEI="{params[1]}"
        results_report="{params.resultsreport}"
        annotation="{params[2]}"
        if [[ "$SV" == "True" ]] && [[ "$MEI" == "True" ]] && [[ "$annotation" == "True" ]] && [[ "$results_report" == "True" ]]; then
            cpan YAML::XS
            cpan Sort::Key::Natural
            export ANNOTSV={config[ANNOTSV_DIR]}
            perl {config[KNOTANNOTSV_DIR]}knotAnnotSV.pl --configFile {config[KNOTANNOTSV_DIR]}config_AnnotSV.yaml --annotSVfile {input[0]} --outDir {params.out_dir} --genomeBuild {params.ref_version}
        fi
        """

rule SVreport:
    input:
        rules.SVannotation.output.SV_annotation_file
    output:
        SV_report = reports_dir + "{sample}/{sample}_SV_annotated.html"
    conda:
        "envs/SVannotation.yaml"
    params:
        SV,
        annotation,
        out_dir = reports_dir,
        sample = "{sample}",
        ref_version = annovar_ref_version,
        resultsreport = config["ANNOTATION_RESULTS_REPORT"]
    log:
        log_dir + "{sample}/SVreport.log"
    resources:
        mem_mb = memory
    shell:
        """
        SV="{params[0]}"
        results_report="{params.resultsreport}"
        annotation="{params[1]}"
        if [[ "$SV" == "True" ]] && [[ "$annotation" == "True" ]] && [[ "$results_report" == "True" ]]; then
            cpan YAML::XS
            cpan Sort::Key::Natural
            export ANNOTSV={config[ANNOTSV_DIR]}
            perl {config[KNOTANNOTSV_DIR]}knotAnnotSV.pl --configFile {config[KNOTANNOTSV_DIR]}config_AnnotSV.yaml --annotSVfile {input[0]} --outDir {params.out_dir} --genomeBuild {params.ref_version}
        fi
        """

rule MEIreport:
    input:
        rules.MEIannotation.output.MEI_annotation_file
    output:
        MEI_report = reports_dir + "{sample}/{sample}_MEI_annotated.html"
    conda:
        "envs/SVannotation.yaml"
    params:
        MEI,
        annotation,
        out_dir = reports_dir,
        sample = "{sample}",
        ref_version = annovar_ref_version,
        resultsreport = config["ANNOTATION_RESULTS_REPORT"]
    log:
        log_dir + "{sample}/MEIreport.log"
    resources:
        mem_mb = memory
    shell:
        """
        MEI="{params[0]}"
        results_report="{params.resultsreport}"
        annotation="{params[1]}"
        if [[ "$MEI" == "True" ]] && [[ "$annotation" == "True" ]] && [[ "$results_report" == "True" ]]; then
            cpan YAML::XS
            cpan Sort::Key::Natural
            export ANNOTSV={config[ANNOTSV_DIR]}
            perl {config[KNOTANNOTSV_DIR]}knotAnnotSV.pl --configFile {config[KNOTANNOTSV_DIR]}config_AnnotSV.yaml --annotSVfile {input[0]} --outDir {params.out_dir} --genomeBuild {params.ref_version}
        fi
        """

rule conciseresultsreport:
    input:
        variantresults = rules.variantresultsreport.output.SNPindel_annotation_report,
        expansionresults = rules.expansionresultsreport.output.expansion_annotation_report if expansion == "true" else [],
        STRresults = rules.STRresultsreport.output.STR_annotation_report if (STR and genotypeSTR) == "true" else [],
        SVresults = rules.SVandMEIreport.SV_MEI_report if (SV and MEI) == "true" else (rules.SVreport.output.SV_report if SV == "true" else (rules.MEIreport.output.MEI_report if MEI == "true" else []))
    output:
        concise_report = reports_dir + "{sample}/{sample}_all_variants.txt"
    params:
        annotation,
        expansion,
        variantcalling,
        STR,
        genotypeSTR,
        SV,
        MEI,
        out_dir = reports_dir,
        sample = "{sample}",
        resultsreport = config["ANNOTATION_RESULTS_REPORT"],
        tempSVMEI_variants = reports_dir + "{sample}/{sample}_temp_SVMEI_variants.tsv",
        tempSNPindel_variants = reports_dir + "{sample}/{sample}_temp_SNVindel_variants.tsv",
        tempexpansion_variants = reports_dir + "{sample}/{sample}_temp_expansion_variants.tsv",
        tempSTR_variants = reports_dir + "{sample}/{sample}_temp_STR_variants.tsv"
    log:
        log_dir + "{sample}/conciseresultsreport.log"
    resources:
        mem_mb = memory
    shell:
        """
        results_report="{params.resultsreport}"
        annotation="{params[0]}"
        expansion="{params[1]}"
        variantcalling="{params[2]}"
        STR="{params[3]}"
        genotype="{params[4]}"
        SV="{params[5]}"
        MEI="{params[6]}"
        if [[ "$results_report" == "True" ]] && [[ "$annotation" == "True" ]] && [[ "$variantcalling" == "True" ]] && [[ "$SV" == "True" ]] && [[ "$MEI" == "True" ]]; then
            python scripts/concisereportSNPindelSVandMEI.py {input.SV_results} {params.tempSVMEI_variants} {input.variantresults} {params.tempSNPindel_variants} scripts/all_variants_report_header.txt {output.concise_report}
            if [[ "$expansion" == "True" ]] && [[ "$genotype" = "True" ]]; then
                python scripts/concisereportexpansionandSTR.py {input.expansionresults} {params.tempexpansion_variants} {input.STRresults} {params.tempSTR_variants} {output.concise_report}
                need to make expansion and then add genotype onto that before adding onto snpindel
            elif [[ "$expansion" == "True" ]]; then
                python scripts/concisereportexpansion.py {input.expansionresults} {params.tempexpansion_variants} {output.concise_report}
            elif [[ "$genotype" == "True" ]]; then
                python scripts/concisereportSTR.py {input.STRresults} {params.tempSTR_variants} {output.concise_report}
            fi
        fi
        """

rule AGSscoring:
    input:
        rules.alignment.output.bam_file,
        rules.variantcalling.output.variant_results_file,
        rules.variantannotation.output.annotated_variant_results_text,
        rules.expansion.output.expansion_file,
        rules.SV.output.merged_SV,
        rules.expansionannotation.output.annotated_expansion_file,
        rules.SVannotation.output.SV_annotation_file,
        rules.variantresultsreport.output.SNPindel_annotation_report,
        rules.expansionresultsreport.output.expansion_annotation_report,
        rules.SVreport.output.SV_report
    output:
        alsgenescanner_all = reports_dir + "{sample}/alsgenescanner/{sample}_alsgenescanner_all.txt",
        alsgenescanner_alsod = reports_dir + "{sample}/alsgenescanner/{sample}_alsgenescanner_alsod.txt",
        alsgenescanner_clinvar = reports_dir + "{sample}/alsgenescanner/{sample}_alsgenescanner_clinvar.txt",
        alsgenescanner_manual_review = reports_dir + "{sample}/alsgenescanner/{sample}_alsgenescanner_manual_review.txt",
        alsgenescanner_ranked = reports_dir + "{sample}/alsgenescanner/{sample}_alsgenescanner_all_ranked.txt"
    conda:
        "envs/AGS.yaml"
    params:
        variantcalling,
        alsgenescanner,
        annotation,
        alsod_list = "resources/alsgenescanner/list_genes_alsod.txt",
        clinvar_list = "resources/alsgenescanner/list_genes_clinvar.txt",
        review_list = "resources/alsgenescanner/list_genes_manual_review.txt",
        sample = "{sample}",
    log:
        log_dir + "{sample}/runAGS.log"
    resources:
        mem_mb = memory
    shell:
        """
        variantcalling="{params[0]}"
        alsgenescanner="{params[1]}"
        annotation="{params[2]}"
        if [[ "$variantcalling" == "True" ]] && [[ "$alsgenescanner" == "True" ]] && [[ "$annotation" == "True" ]]; then
            python scripts/alsgenescanner.py {input[0]} {output.alsgenescanner_all}
            python scripts/alsgenescannerreport.py {output.alsgenescanner_all} {output.alsgenescanner_alsod} {params.alsod_list} {output.alsgenescanner_clinvar} {params.clinvar_list} {output.alsgenescanner_manual_review} {params.review_list} {output.alsgenescanner_ranked}
        fi
        """
