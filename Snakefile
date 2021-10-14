configfile: "config.yaml"

import os
import os.path
import sys
import glob
import gzip
import bgzip
#install bedtools

#INPUT OPTIONS
##bam, cram, fastq and vcf options
format = config["INPUT_FORMAT"]
sample_name = config["SAMPLE_NAMES"]
sample = "{sample}"
alsgenescanner = config["ALSGENESCANNER"]
reference = config["REFERENCE_VERSION"]
paired = config["READ_TYPE"]
input_dir = config["INPUT_FILE_DIR"]
vcf = config["ADDITIONAL_VCF_FILE_DIR"]
filter_string = config["FILTER_VARIANTS"]
debug = config["DEBUG"]
exome = config["EXOME"]
BED = config["BED"]
RG = config["ADD_READ_GROUP"]
rm_dup = config["REMOVE_DUPLICATES"]

#MAIN OUTPUT DIRECTORIES
results_dir = config["OUT_DIR"]
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
custom_microbes = config["CUSTOM_MICROBES"]
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
        path_bed = "resources/alsgenescanner/als_gene_scanner_hg19.bed"
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
        path_bed = "resources/alsgenescanner/als_gene_scanner_hg38.bed"
    if reference == "grch38":
        os.system("zcat resources/exome_hg38.bed.gz | sed 's/chr//g' | bgzip -c > resources/exome_grch38.bed.gz")
        os.system("zcat resources/hg38_gene_db.txt.gz | sed 's/chr//g' | bgzip -c > resources/grch38_gene_db.txt.gz")
        os.system("cp resources/hg38_gene_names.txt.gz resources/grch38_gene_names.txt.gz")

variant_results_file = ""

if format == "vcf":
    variant_results_file = input_file

    if not vcf:
        vcf = variant_results_file

#ADAPTING INPUT FILE FORMATS - unsure about how to get this
if format == "fastq":
    input_files = expand(input_dir + "{sample}.1.fq.gz", sample=sample_name) + expand(input_dir + "{sample}.2.fq.gz", sample=sample_name)
    input_file = input_dir + "{sample}.1.fq.gz"
    if paired == "paired":
        input_file2 = input_dir + "{sample}.2.fq.gz"
if format == "sam":
    input_files = expand(input_dir + "{sample}.sam", sample=sample_name)
if format == "bam":
    input_files = expand(input_dir + "{sample}.bam", sample=sample_name)
    input_file = input_dir + "{sample}.bam"
if format = "cram":
    input_files = expand(input_dir + "{sample}.cram", sample=sample_name)
    input_file = input_dir + "{sample}.bam"
if format == "vcf":
    input_files = expand(input_dir + "{sample}.vcf.gz", sample=sample_name)

if alsgenescanner == "true":
    annovar_protocols = "refGene,dbnsfp33a,clinvar_20210501,intervar_20180118"
    annovar_operations = "g,f,f,f"
    path_gene_list = ""
    BED = "true"
    alignment = "true"
    variantcalling = "true"
    annotation = "true"
    expansion = "true"
    SV = "true"

if rm_dup == "true":
    if exome == "true":
        samblaster_cmq = "samblaster --ignoreUnmated |"
    else:
        samblaster_cmq = "samblaster |"
        samblaster_bwa = "samblaster --ignoreUnmated |"

#sorting out the gene list to see if there are any unmatched genes in the reference and to make a custom bed out of that if no bed file is provided
if BED or path_gene_list:
    if len(path_bed) == 0:
        if len(path_gene_list) != 0:
            os.system("zgrep -iwf %s resources/%s_gene_names.txt.gz | awk '{print $2}' > %s/matched_genes.txt" % (path_gene_list, reference, results_dir))
            os.system("zgrep -viwf %s/matched_genes.txt %s > %s/unmatched_genes.txt" % (results_dir, path_gene_list, results_dir))
            os.system("zgrep -iwf %s resources/%s_gene_names.txt.gz | awk '{print $1}' > %s/matched_genes_codes.txt" % (path_gene_list, reference, results_dir))
            os.system("zgrep -wf %s/matched_genes_codes.txt resources/%s_gene_db.txt | awk '{i=1; while (i<= int($8)) {n=split($9,a,/,/);n=split($10,b,/,/); print $2\"\t\"a[i]\"\t\"b[i]; i+=1}}' > %s/custom_tmp.bed" % (results_dir, reference, results_dir))
            os.system("bedtools sort -i %s/custom_tmp.bed > %s/custom_sorted.bed" % (out, out))
            os.system("bedtools merge -i %s/custom_sorted.bed > %s/custom.bed" % (out, out))
            os.system("rm %s/custom_sorted.bed %s/custom_tmp.bed" % (out, out))
            BED = "true"
            path_bed = "%s/custom.bed" % (results_dir)

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

#then do the rule all - need to think about this


#ALIGNMENT
if format == "fastq" and alignment == "true":
    if paired == "paired":
        if variantcalling and not (SV or MEI or STR or genotypeSTR or expansion):
            rule alignment:
                input:
                    fastq1 = input_file,
                    fastq2 = input_file2
                output:
                    bam_file = results_dir + "{sample}/{sample}_sorted.bam",
                    bam_file_index = results_dir + "{sample}/{sample}_sorted.bam.bai"
                conda:
                    "envs/alignmentfast.yaml"
                resources:
                    mem_gb = ?
                log:
                    log_dir + "alignmentpairedSNPindel.log"
                shell:
                    """
                    hisat2 {config[HISAT_CUSTOM_OPTIONS]} --no-spliced-alignment -p {config[NUMBER_CPU]} -x {config[HISAT2_INDEX]} -1 {input.fastq1} -2 {input.fastq2} | {samblaster_cmq} samtools view -@ {config[NUMBER_CPU]} -Sb - | sambamba sort -t {config[NUMBER_CPU]} --tmpdir={tmp_dir} -o {output.bam_file} /dev/stdin
                    samtools index -@ {config[NUMBER_CPU]} {output.bam_file}
                    """

        else:
            if (SV or MEI or STR or genotypeSTR or expansion):
                rule alignment:
                    input:
                        fastq1 = input_file,
                        fastq2 = input_file2
                    output:
                        bam_file = results_dir + "{sample}/{sample}_sorted_merged.bam",
                        bam_file_index = results_dir + "{sample}/{sample}_sorted_merged.bam.bai"
                    conda:
                        "envs/alignmentnormal.yaml"
                    params:
                        hisat2_bam = results_dir + "{sample}/{sample}_hisat2.bam",
                        unaligned_reads = results_dir + "{sample}/{sample}_unaligned_reads.fq",
                        bwa_bam = results_dir + "{sample}/{sample}_bwa.bam",
                        header = results_dir + "{sample}/header.txt",
                        sample = "{sample}",
                        out_dir = results_dir,
                        rg_hisat2 = rg_option_hisat2,
                        rg_bwa = rg_option_bwa
                    resources:
                        mem_gb = ?
                    log:
                        log_dir + "alignmentpaired.log"
                    shell:
                        """
                        hisat2 {config[HISAT_CUSTOM_OPTIONS]} {params.rg_hisat2} --no-softclip --no-spliced-alignment -p {config[NUMBER_CPU]} -x {config[HISAT2_INDEX]} -1 {input.fastq1} -2 {input.fastq2} | {samblaster_cmq} samtools view -@ {config[NUMBER_CPU]} -Sb - | sambamba sort -t {config[NUMBER_CPU]} --tmpdir={tmp_dir} -o {params.hisat2_bam} /dev/stdin
                        samtools index -@ {config[NUMBER_CPU]} {params.hisat2_bam}
                        samtools view -@ {config[NUMBER_CPU]} -bhf 4 {params.hisat2_bam} | samtools bam2fq - > {params.unaligned_reads}
                        bwa mem {config[BWA_CUSTOM_OPTIONS]} {params.rg_bwa} -t {config[NUMBER_CPU]} {config[BWA_INDEX]} {params.unaligned_reads} | {samblaster_bwa} samtools view -@ {config[NUMBER_CPU]} -Sb - | sambamba sort -t {config[NUMBER_CPU]} --tmpdir={tmp_dir} -o {params.bwa_bam} /dev/stdin
                        samtools index -@ {config[NUMBER_CPU]} {params.bwa_bam}
                        samtools view -H {params.hisat2_bam} > {params.header}
                        samtools merge -c -@ {config[NUMBER_CPU]} -f -h {params.header} {output.bam_file} {params.hisat2_bam} {params.bwa_bam}
                        rm {params.unaligned_reads} {params.header} {params.hisat2_bam} {params.bwa_bam} {params.out_dir}{params.sample}/{params.sample}_hisat2.bam.bai {params.out_dir}{params.sample}/{params.sample}_bwa.bam.bai
                        samtools index -@ {config[NUMBER_CPU]} {output.bam_file}
                        """

    else:
        if variantcalling and not (SV or STR or MEI or expansion or genotypeSTR):
            rule alignment:
                input:
                    input_file
                output:
                    bam_file = results_dir + "{sample}/{sample}_sorted.bam",
                    bam_file_index = results_dir + "{sample}/{sample}_sorted.bam.bai"
                conda:
                    "envs/alignmentfast.yaml"
                resources:
                    mem_gb = ?
                log:
                    log_dir + "alignmentsinglefast.log"
                shell:
                    """
                    hisat2 {config[HISAT_CUSTOM_OPTIONS]} --no-spliced-alignment -p {config[NUMBER_CPU]} -x {config[HISAT2_INDEX]} -U {input[0]} | {samblaster_cmq} samtools view -@ {config[NUMBER_CPU]} -Sb - | sambamba sort -t {config[NUMBER_CPU]} --tmpdir={tmp_dir} -o {output.bam_file} /dev/stdin
                    samtools index -@ {config[NUMBER_CPU]} {output.bam_file}
                    """
        else:
            if (SV or MEI or STR or genotypeSTR or expansion):
                rule alignment:
                    input:
                        input_file
                    output:
                        bam_file = results_dir + "{sample}/{sample}_sorted_merged.bam",
                        bam_file_index = results_dir + "{sample}/{sample}_sorted_merged.bam.bai"
                    conda:
                        "envs/alignmentnormal.yaml"
                    params:
                        hisat2_bam = results_dir + "{sample}/{sample}_hisat2.bam",
                        unaligned_reads = results_dir + "{sample}/{sample}_unaligned_reads.fq",
                        bwa_bam = results_dir + "{sample}/{sample}_bwa.bam",
                        header = results_dir + "{sample}/header.txt",
                        sample = "{sample}",
                        out_dir = results_dir,
                        rg_hisat2 = rg_option_hisat2,
                        rg_bwa = rg_option_bwa
                    resources:
                        mem_gb = ?
                    log:
                        log_dir + "alignmentsingle.log"
                    shell:
                        """
                        hisat2 {config[HISAT_CUSTOM_OPTIONS]} {params.rg_hisat2} --no-softclip --no-spliced-alignment -p {config[NUMBER_CPU]} -x {config[HISAT2_INDEX]} -1 {input.fastq1} -2 {input.fastq2} | {samblaster_cmq} samtools view -@ {config[NUMBER_CPU]} -Sb - | sambamba sort -t {config[NUMBER_CPU]} --tmpdir={tmp_dir} -o {params.hisat2_bam} /dev/stdin
                        samtools index -@ {config[NUMBER_CPU]} {params.hisat2_bam}
                        samtools view -@ {config[NUMBER_CPU]} -bhf 4 {params.hisat2_bam} | samtools bam2fq - > {params.unaligned_reads}
                        bwa mem {config[BWA_CUSTOM_OPTIONS]} {params.rg_bwa} -t {config[NUMBER_CPU]} {config[BWA_INDEX]} {params.unaligned_reads} | {samblaster_bwa} samtools view -@ {config[NUMBER_CPU]} -Sb - | sambamba sort -t {config[NUMBER_CPU]} --tmpdir={tmp_dir} -o {params.bwa_bam} /dev/stdin
                        samtools index -@ {config[NUMBER_CPU]} {params.bwa_bam}
                        samtools view -H {params.hisat2_bam} > {params.header}
                        samtools merge -c -@ {config[NUMBER_CPU]} -f -h {params.header} {output.bam_file} {params.hisat2_bam} {params.bwa_bam}
                        rm {params.unaligned_reads} {params.header} {params.hisat2_bam} {params.bwa_bam} {params.out_dir}{params.sample}/{params.sample}_hisat2.bam.bai {params.out_dir}{params.sample}/{params.sample}_bwa.bam.bai
                        samtools index -@ {config[NUMBER_CPU]} {output.bam_file}
                        """

#SAMTOBAM CONVERSION AND DEFINING BAM FILES FOR VARIANT CALLING
if format == "sam":
    rule sam2bam:
        input:
            input_file = input_dir + "{sample}.sam"
        output:
            bam_file = results_dir + "{sample}/{sample}_sorted.bam",
            bam_file_index = results_dir + "{sample}/{sample}_sorted.bam.bai"
        conda:
            "envs/samtools.yaml"
        log:
            log_dir + "samtobam.log"
        shell:
            """
            samtools view -Sb {input[0]} > {output.bam_file}
            samtools index -@ {config[NUMBER_CPU]} {output.bam_file}
            """

if format == "bam" or format == "cram":
    bam_file = input_file

if format == "cram" and SV == "true":
    rule CramToBam:
        input:
            path_reference,
            input_file = input_dir + "{sample}.cram"
        output:
            bam_file = results_dir + "{sample}/{sample}_delly.bam"
            bam_file = results_dir + "{sample}/{sample}_delly.bam.bai"
        conda:
            "envs/cramtobam.yaml"
        log:
            log_dir + "cramtobam.log"
        shell:
            """
            samtools view -b -h -@ {config[NUMBER_CPU]} -T {input[0]} -o {output.delly_bam} {input.input_file}
            samtools index -@ {config[NUMBER_CPU]} {output.delly_bam}
            """

# VARIANT CALLING WITH STRELKA2 (requires paired-end)
if variantcalling == "true":
    if paired == "paired":
        if BED == "true":
            rule variantcalling:
                input:
                    bam_file,
                    path_bed,
                    path_reference
                output:
                    variant_results_file_unfiltered = results_dir + "{sample}/{sample}_sorted.vcf.gz"
                    variant_results_file_unfiltered_index = results_dir + "{sample}/{sample}_sorted.vcf.gz.tbi"
                conda:
                    "envs/variantcalling.yaml"
                params:
                    temp_bed = results_dir + "{sample}/temp.bed.gz",
                    sorted_bed = results_dir + "{sample}/sorted.bed.gz",
                    out_dir = results_dir,
                    sample = "{sample}"
                log:
                    log_dir + "variantcallingwithBED.log"
                shell:
                    """
                    bgzip -c {input[1]} > {params.temp_bed}
                    sortBed -i {params.temp_bed} | bgzip -c > {params.sorted_bed}
                    tabix -p bed {params.sorted_bed}
                    {config[STRELKA_DIR]}bin/configureStrelkaGermlineWorkflow.py --bam {input[0]} --referenceFasta {input[2]} --runDir {params.out_dir}/{params.sample}/strelka --callRegions {params.sorted_bed}
                    {params.out_dir}{params.sample}/strelka/runWorkflow.py -j {config[NUMBER_CPU]} -m local
                    mv {params.out_dir}{params.sample}/strelka/results/variants/genome.S1.vcf.gz {output.variant_results_file_unfiltered}
                    mv {params.out_dir}{params.sample}/strelka/results/variants/genome.S1.vcf.gz.tbi {output.variant_results_file_unfiltered_index}
                    """

        if exome == "true":
            rule variantcalling:
                input:
                    bam_file,
                    path_reference
                output:
                    variant_results_file_unfiltered = results_dir + "{sample}/{sample}_sorted.vcf.gz",
                    variant_results_file_unfiltered_index = results_dir + "{sample}/{sample}_sorted.vcf.gz.tbi"
                conda:
                    "envs/variantcalling.yaml"
                params:
                    out_dir = results_dir,
                    sample = "{sample}"
                log:
                    log_dir + "variantcallingwithexome.log"
                shell:
                    """
                    {config[STRELKA_DIR]}bin/configureStrelkaGermlineWorkflow.py --bam {input[0]} --referenceFasta {input[1]} --runDir {params.out_dir}/{params.sample}/strelka --exome
                    {params.out_dir}{params.sample}/strelka/runWorkflow.py -j {config[NUMBER_CPU]} -m local
                    mv {params.out_dir}{params.sample}/strelka/results/variants/genome.S1.vcf.gz {output.variant_results_file_unfiltered}
                    mv {params.out_dir}{params.sample}/strelka/results/variants/genome.S1.vcf.gz.tbi {output.variant_results_file_unfiltered_index}
                    """
        else:
            rule variantcalling:
                input:
                    bam_file,
                    path_reference
                output:
                    variant_results_file_unfiltered = results_dir + "{sample}/{sample}_sorted.vcf.gz",
                    variant_results_file_unfiltered_index = results_dir + "{sample}/{sample}_sorted.vcf.gz.tbi"
                conda:
                    "envs/variantcalling.yaml"
                params:
                    out_dir = results_dir,
                    sample = "{sample}"
                log:
                    log_dir + "variantcalling.log"
                shell:
                    """
                    {config[STRELKA_DIR]}bin/configureStrelkaGermlineWorkflow.py --bam {input[0]} --referenceFasta {input[1]} --runDir {params.out_dir}/{params.sample}/strelka
                    {params.out_dir}{params.sample}/strelka/runWorkflow.py -j {config[NUMBER_CPU]} -m local
                    mv {params.out_dir}{params.sample}/strelka/results/variants/genome.S1.vcf.gz {output.variant_results_file_unfiltered}
                    mv {params.out_dir}{params.sample}/strelka/results/variants/genome.S1.vcf.gz.tbi {output.variant_results_file_unfiltered_index}
                    rm -r {params.out_dir}{params.sample}/strelka
                    """

#need to find a way to annotate the unfiltered variant file

if filter_string == "true":
    rule variantfilter:
        input:
            rules.variantcalling.output.variant_results_file_unfiltered,
        output:
            variant_results_file_filtered = results_dir + "{sample}/{sample}_sorted_filtered.vcf.gz",
            variant_results_file_filtered_index = results_dir + "{sample}/{sample}_sorted_filtered.vcf.gz.tbi"
        conda:
            "envs/variantfilter.yaml"
        log:
            log_dir + "filtervariants.log"
        shell:
            """
            bcftools filter -i '{config[VARIANT_FILTER_STRING]}' {input[0]} | bgzip -c > {output.variant_results_file_filtered} ; tabix -fp vcf {output.variant_results_file_filtered}
            """

    if annotation == "true":
        rule variantannotation:
            input:
                rules.variantfilter.output.variant_results_file_filtered
            output:
                annotated_variant_results_file = results_dir + "{sample}/{sample}_SNPindel_annotated.vcf.gz",
                annotated_variant_results_file_index = results_dir + "{sample}/{sample}_SNPindel_annotated.vcf.gz.tbi",
                annotated_variant_results_text = results_dir + "{sample}/{sample}_SNPindel_annotated.txt"
            conda:
                "envs/annotation.yaml"
            params:
                out_dir = results_dir,
                sample = "{sample}",
                ref_version = annovar_ref_version,
                operations = annovar_operations,
                protocols = annovar_protocols
            log:
                log_dir + "variantannotation.log"
            shell:
                """
                perl {config[ANNOVAR_DIR]}table_annovar.pl --thread {config[NUMBER_CPU]} --vcfinput {input[0]} {config[ANNOVAR_DB]} -buildver {params.ref_version} -remove -protocol {params.protocols} -operation {params.operations} -nastring . --outfile {params.out_dir}{params.sample}/{sample}_annovar_SNPindel.vcf
                mv {params.out_dir}{params.sample}/{sample}_annovar_SNPindel.vcf.{params.ref_version}_multianno.vcf {params.out_dir}{params.sample}/{params.sample}_SNPindel_annotated.vcf
                mv {params.out_dir}{params.sample}/annovar_SNPindel.vcf.{params.ref_version}_multianno.txt {output.annotated_variant_results_text}
                bgzip -f {params.out_dir}{params.sample}/{params.sample}_SNPindel_annotated.vcf ; tabix -fp vcf {params.out_dir}{params.sample}/{params.sample}_SNPindel_annotated.vcf.gz
                """
else:
    if annotation == "true":
        rule variantannotation:
            input:
                rules.variantfilter.output.variant_results_file_unfiltered
            output:
                annotated_variant_results_file = results_dir + "{sample}/{sample}_SNPindel_annotated.vcf.gz",
                annotated_variant_results_file = results_dir + "{sample}/{sample}_SNPindel_annotated.vcf.gz.tbi",
                annotated_variant_results_text = results_dir + "{sample}/{sample}_SNPindel_annotated.txt"
            conda:
                "envs/annotation.yaml"
            params:
                out_dir = results_dir,
                sample = "{sample}",
                ref_version = annovar_ref_version,
                operations = annovar_operations,
                protocols = annovar_protocols
            log:
                log_dir + "variantannotation.log"
            shell:
                """
                perl {config[ANNOVAR_DIR]}table_annovar.pl --thread {config[NUMBER_CPU]} --vcfinput {input[0]} {config[ANNOVAR_DB]} -buildver {params.ref_version} -remove -protocol {params.protocols} -operation {params.operations} -nastring . --outfile {params.out_dir}{params.sample}/annovar_SNPindel.vcf
                mv {params.out_dir}{params.sample}/annovar_SNPindel.vcf.{params.ref_version}_multianno.vcf {params.out_dir}{params.sample}/{params.sample}_SNPindel_annotated.vcf
                mv {params.out_dir}{params.sample}/annovar_SNPindel.vcf.{params.ref_version}_multianno.txt {output.annotated_variant_results_text}
                bgzip -f {params.out_dir}{params.sample}/{params.sample}_SNPindel_annotated.vcf ; tabix -fp vcf {params.out_dir}{params.sample}/{params.sample}_SNPindel_annotated.vcf.gz
                """

if expansion == "true":
    rule expansion:
        input:
            bam_file,
            path_reference
        output:
            expansion_file = results_dir + "{sample}/{sample}_expansions.vcf.gz"
            expansion_file_index = results_dir + "{sample}/{sample}_expansions.vcf.gz.tbi"
        conda:
            "envs/expansion.yaml"
        params:
            out_dir = results_dir,
            sample = "{sample}",
            variant_catalog = path_expansionHunter_catalog
        log:
            log_dir + "expansion.log"
        shell:
            """
            ExpansionHunter --reads {input[0]} --reference {input[1]} --variant-catalog {params.variant_catalog} --output-prefix {params.out_dir}{params.sample}_expansions
            bgzip {params.out_dir}{params.sample}_expansions.vcf
            tabix -p vcf {output.expansion_file}
            """

    if annotation == "true":
        rule expansionannotation:
            input:
                rules.expansion.output.expansion_file
            output:
                annotated_expansion_file = results_dir + "{sample}/{sample}_expansions_annotated.vcf.gz"
                annotated_expansion_file_index = results_dir + "{sample}/{sample}_expansions_annotated.vcf.gz.tbi"
            conda:
                "envs/annotation.yaml"
            params:
                out_dir = results_dir,
                sample = "{sample}",
                ref_version = annovar_ref_version,
                operations = annovar_operations,
                protocols = annovar_protocols
            log:
                log_dir + "expansionannotation.log"
            shell:
                """
                perl {config[ANNOVAR_DIR]}table_annovar.pl --thread {config[NUMBER_CPU]} --vcfinput {input[0]} {config[ANNOVAR_DB]} -buildver {params.ref_version} -remove -protocol {params.protocols} -operation {params.operations} -nastring . --outfile {params.out_dir}{params.sample}/annovar_expansions.vcf
                mv {params.out_dir}{params.sample}/annovar_expansions.vcf.{params.ref_version}_multianno.vcf {params.out_dir}{params.sample}/{params.sample}_expansions_annotated.vcf
                bgzip -f {params.out_dir}{params.sample}/{params.sample}_expansions_annotated.vcf ; tabix -fp vcf {params.out_dir}{params.sample}/{params.sample}_expansions_annotated.vcf.gz
                """

if STR == "true":
    rule STRprofile:
        input:
            bam_file,
            path_reference
        output:
            STR_profile = results_dir + "{sample}/{sample}_expansiondenovo.str_profile.json"
            genotypeSTR_input = results_dir + "{sample}/{sample}_genotypeSTRinput.txt"
        conda:
            "envs/expansion.yaml"
        params:
            out_dir = results_dir,
            sample = "{sample}"
        log:
            log_dir + "STRprofile.log"
        shell:
            """
            {config[EXPANSIONHUNTERDENOVO_DIR]}bin/ExpansionHunterDenovo profile --reads {input[0]} --reference {input[1]} --output-prefix {params.out_dir}{params.sample}/{params.sample}_expansiondenovo --min-anchor-mapq 50 --max-irr-mapq 40 --log-reads
            cat {params.out_dir}{params.sample}/{params.sample}_expansiondenovo.locus.tsv | sed 's/contig/chr/g' | cut -f1-4 > {output.genotypeSTR_input}
            """

if genotypeSTR == "true":
    rule genotypeSTR:
        input:
            rules.STRprofile.output.genotypeSTR_input,
            path_reference,
            bam_file
        output:
            EHDN_variant_catalog = results_dir + "{sample}/{sample}_EHDN_variant_catalog.json",
            EHDN_expansion_file = results_dir + "{sample}/{sample}_EHDNexpansions.vcf.gz",
            EHDN_expansion_file_index = results_dir + "{sample}/{sample}_EHDNexpansions.vcf.gz"
        conda:
            "envs/expansion.yaml"
        params:
            out_dir = results_dir,
            sample = "{sample}",
            EHDN_excluded = results_dir + "{sample}/{sample}_EHDN_excluded.csv",
            EHDN_unmatched = results_dir + "{sample}/{sample}_EHDN_unmatched.csv"
        log:
            log_dir + "genotypeSTR.log"
        shell:
            """
            python scripts/conversion_EHDN_catalog.py {input[0]} {input[1]} {output.EHDN_variant_catalog} {params.EHDN_unmatched} {params.EHDN_excluded}
            ExpansionHunter --reads {input[2]} --reference {input[1]} --variant-catalog {output.EHDN_variant_catalog} --output-prefix {params.out_dir}{params.sample}_EHDNexpansions
            bgzip {params.out_dir}{params.sample}_EHDNexpansions.vcf
            tabix -p vcf {output.EHDN_expansion_file}
            """

    if annotation == "true":
        rule STRannotation:
            input:
                rules.genotypeSTR.output.EHDN_expansion_file
            output:
                annotated_EHDN_expansion_file = results_dir + "{sample}/{sample}_EHDNexpansions_annotated.vcf.gz",
                annotated_EHDN_expansion_file_index = results_dir + "{sample}/{sample}_EHDNexpansions_annotated.vcf.gz.tbi"
            conda:
                "envs/annotation.yaml"
            params:
                out_dir = results_dir,
                sample = "{sample}",
                ref_version = annovar_ref_version,
                protocols = annovar_protocols,
                operations = annovar_operations
            log:
                log_dir + "STRannotation.log"
            shell:
                """
                perl {config[ANNOVAR_DIR]}table_annovar.pl --thread {config[NUMBER_CPU]} --vcfinput {input[0]} {config[ANNOVAR_DB]} -buildver {params.ref_version} -remove -protocol {params.protocols} -operation {params.operations} -nastring . --outfile {params.out_dir}{params.sample}/annovar_EHDNexpansions.vcf
                mv {params.out_dir}{params.sample}/annovar_EHDNexpansions.vcf.{params.ref_version}_multianno.vcf {params.out_dir}{params.sample}/{params.sample}_EHDNexpansions_annotated.vcf
                bgzip -f {params.out_dir}{params.sample}/{params.sample}_EHDNexpansions_annotated.vcf ; tabix -fp vcf {params.out_dir}{params.sample}/{params.sample}_EHDNexpansions_annotated.vcf.gz
                """

if SV == "true":
    if paired == "paired":
        if BED == "true":
            rule SV:
                input:
                    bam_file,
                    path_reference,
                    path_bed
                output:
                    manta_SV = results_dir + "{sample}/{sample}_manta_SV.vcf",
                    delly_SV = results_dir + "{sample}/{sample}_delly_SV.vcf",
                    merged_SV = results_dir + "{sample}/{sample}_merged_SV.vcf.gz",
                    merged_SV_index = results_dir + "{sample}/{sample}_merged_SV.vcf.gz.tbi"
                conda:
                    "envs/variantcallingwithBED.yaml"
                params:
                    temp_bed = results_dir + "{sample}/temp.bed.gz",
                    sorted_bed = results_dir + "{sample}/sorted.bed.gz",
                    out_dir = results_dir,
                    sample = "{sample}",
                    delly_exclude_regions = path_delly_exclude_regions
                log:
                    log_dir + "SVwithBED.log"
                shell:
                    """
                    bgzip -c {input[2]} > {params.temp_bed}
                    sortBed -i {params.temp_bed} | bgzip -c > {params.sorted_bed}
                    tabix -p bed {params.sorted_bed}
                    {config[MANTA_DIR]}bin/configManta.py --bam {input[0]} --referenceFasta {input[1]} --runDir {params.out_dir}{params.sample}/SV_manta --callRegions {params.sorted_bed}
                    {params.out_dir}{params.sample}SV_manta/runWorkflow.py -j {config[NUMBER_CPU]} -m local
                    mv {params.out_dir}{params.sample}/SV_manta/results/variants/diploidSV.vcf.gz {params.out_dir}{params.sample}/{params.sample}_SV_manta.vcf.gz
                    gzip -d {params.out_dir}{params.sample}/{params.sample}_SV_manta.vcf.gz
                    rm -r {params.out_dir}{params.sample}/SV_manta

                    delly call -g {input[1]} -o {params.out_dir}{params.sample}/{params.sample}_delly_SV.bcf -x {params.delly_exclude_regions} {input[0]}
                    bcftools view {params.out_dir}{params.sample}/{params.sample}_delly_SV.bcf > {output.delly_SV}

                    ls {params.out_dir}{params.sample}/*SV.vcf > {params.out_dir}{params.sample}/survivor_sample_files
                    SURVIVOR merge {params.out_dir}{params.sample}/survivor_sample_files 1000 1 1 1 0 30 {params.out_dir}{params.sample}/{params.sample}_merged_SV.vcf
                    perl scripts/vcf-sort.pl {params.out_dir}{params.sample}/{params.sample}_merged_SV.vcf | bgzip -c > {output.merged_SV}
                    tabix -p vcf {output.merged_SV}

                    rm {params.out_dir}{params.sample}/survivor_sample_files {params.out_dir}{params.sample}/{params.sample}_delly_SV.bcf {params.out_dir}{params.sample}/{params.sample}_delly_SV.bcf.csi {params.out_dir}{params.sample}/{params.sample}_merged_SV.vcf
                    """

            if alsgenescanner == "true":
                if annotation == "true":
                    rule SVannotation:
                        input:
                            rules.SV.output.merged_SV,
                        output:
                            SV_annotation_file = results_dir + "{sample}/{sample}_SV_annotated.tsv"
                        conda:
                            "envs/SVannotation.yaml"
                        params:
                            out_dir = results_dir,
                            sample = "{sample}",
                            genomebuild = annotsv_ref_version
                        log:
                            log_dir + "SVannotation.log"
                        shell:
                            """
                            cpan YAML::XS
                            cpan Sort::Key::Natural
                            export ANNOTSV={config[ANNOTSV_DIR]}
                            {config[ANNOTSV_DIR]}bin/AnnotSV -annotationsDir {config[ANNOTSV_DIR]}/share/AnnotSV -SVinputFile {input[0]} -genomeBuild {params.genomebuild} -candidateGenesFiltering yes -candidateGenesFile resources/alsgenescanner/list_genes_all.txt outputDir {params.out_dir}{params.sample} -outputFile {params.sample}_SV_annotated -SVminSize 30 {config[ANNOTSV_CUSTOM_OPTIONS]}
                            """
            else:
                if annotation == "true":
                     rule SVannotation:
                        input:
                            rules.SV.output.merged_SV
                        output:
                            SV_annotation_file = results_dir + "{sample}/{sample}_SV_annotated.tsv"
                        conda:
                            "envs/SVannotation.yaml"
                        params:
                            out_dir = results_dir,
                            sample = "{sample}",
                            genomebuild = annotsv_ref_version
                        log:
                            log_dir + "SVannotation.log"
                        shell:
                            """
                            cpan YAML::XS
                            cpan Sort::Key::Natural
                            export ANNOTSV={config[ANNOTSV_DIR]}
                            {config[ANNOTSV_DIR]}bin/AnnotSV -annotationsDir {config[ANNOTSV_DIR]}/share/AnnotSV -SVinputFile {input[0]} -genomeBuild {params.genomebuild} outputDir {params.out_dir}{params.sample} -outputFile {params.sample}_SV_annotated -SVminSize 30 {config[ANNOTSV_CUSTOM_OPTIONS]}
                            """
        else:
            rule SV:
                input:
                    bam_file,
                    path_reference,
                output:
                    manta_SV = results_dir + "{sample}/{sample}_manta_SV.vcf",
                    delly_SV = results_dir + "{sample}/{sample}_delly_SV.vcf",
                    merged_SV = results_dir + "{sample}/{sample}_merged_SV.vcf.gz",
                    merged_SV_index = results_dir + "{sample}/{sample}_merged_SV.vcf.gz.tbi"
                conda:
                    "envs/variantcalling.yaml"
                params:
                    out_dir = results_dir,
                    sample = "{sample}",
                    delly_exclude_regions = path_delly_exclude_regions
                log:
                    log_dir + "SV.log"
                shell:
                    """
                    {config[MANTA_DIR]}bin/configManta.py --bam {input[0]} --referenceFasta {input[1]} --runDir {params.out_dir}{params.sample}/SV_manta
                    {params.out_dir}{params.sample}SV_manta/runWorkflow.py -j {config[NUMBER_CPU]} -m local
                    mv {params.out_dir}{params.sample}/SV_manta/results/variants/diploidSV.vcf.gz {params.out_dir}{params.sample}/{params.sample}_SV_manta.vcf.gz
                    gzip -d {params.out_dir}{params.sample}/{params.sample}_SV_manta.vcf.gz
                    rm -r {params.out_dir}{params.sample}/SV_manta

                    delly call -g {input[1]} -o {params.out_dir}{params.sample}/{params.sample}_delly_SV.bcf -x {params.delly_exclude_regions} {input[0]}
                    bcftools view {params.out_dir}{params.sample}/{params.sample}_delly_SV.bcf > {output.delly_SV}

                    ls {params.out_dir}{params.sample}/*SV.vcf > {params.out_dir}{params.sample}/survivor_sample_files
                    {config[SURVIVOR_DIR]}SURVIVOR merge {params.out_dir}{params.sample}/survivor_sample_files 1000 1 1 1 0 30 {params.out_dir}{params.sample}/{params.sample}_merged_SV.vcf
                    perl scripts/vcf-sort.pl {params.out_dir}{params.sample}/{params.sample}_merged_SV.vcf | bgzip -c > {output.merged_SV}
                    tabix -p vcf {output.merged_SV}

                    rm {params.out_dir}{params.sample}/survivor_sample_files {params.out_dir}{params.sample}/{params.sample}_delly_SV.bcf {params.out_dir}{params.sample}/{params.sample}_delly_SV.bcf.csi {params.out_dir}{params.sample}/{params.sample}_merged_SV.vcf
                    """

            if alsgenescanner == "true":
                if annotation == "true":
                    rule SVannotation:
                        input:
                            rules.SV.output.merged_SV,
                        output:
                            SV_annotation_file = results_dir + "{sample}/{sample}_SV_annotated.tsv"
                        conda:
                            "envs/SVannotation.yaml"
                        params:
                            out_dir = results_dir,
                            sample = "{sample}",
                            genomebuild = annotsv_ref_version
                        log:
                            log_dir + "SVannotation.log"
                        shell:
                            """
                            cpan YAML::XS
                            cpan Sort::Key::Natural
                            export ANNOTSV={config[ANNOTSV_DIR]}
                            {config[ANNOTSV_DIR]}bin/AnnotSV -annotationsDir {config[ANNOTSV_DIR]}/share/AnnotSV -SVinputFile {input[0]} -genomeBuild {params.genomebuild} -candidateGenesFiltering yes -candidateGenesFile resources/alsgenescanner/list_genes_all.txt outputDir {params.out_dir}{params.sample} -outputFile {params.sample}_SV_annotated -SVminSize 30 {config[ANNOTSV_CUSTOM_OPTIONS]}
                            """
            else:
                if annotation == "true":
                    rule SVannotation:
                        input:
                            rules.SV.output.merged_SV
                        output:
                            SV_annotation_file = results_dir + "{sample}/{sample}_SV_annotated.tsv"
                        conda:
                            "envs/SVannotation.yaml"
                        params:
                            out_dir = results_dir,
                            sample = "{sample}",
                            genomebuild = annotsv_ref_version
                        log:
                            log_dir + "SVannotation.log"
                        shell:
                            """
                            cpan YAML::XS
                            cpan Sort::Key::Natural
                            export ANNOTSV={config[ANNOTSV_DIR]}
                            {config[ANNOTSV_DIR]}bin/AnnotSV -annotationsDir {config[ANNOTSV_DIR]}/share/AnnotSV -SVinputFile {input[0]} -genomeBuild {params.genomebuild} outputDir {params.out_dir}{params.sample} -outputFile {params.sample}_SV_annotated -SVminSize 30 {config[ANNOTSV_CUSTOM_OPTIONS]}
                            """

if MEI == "true":
    if exome == "true":
        rule MEI:
            input:
                bam_file,
                path_reference,
                melt_bed
            output:
                MEI_file = results_dir + "{sample}/{sample}_MEI.vcf.gz",
                MEI_file_index = results_dir + "{sample}/{sample}_MEI.vcf.gz.tbi"
            conda:
                "envs/MEI.yaml"
            params:
                out_dir = results_dir,
                sample="{sample}",
                melt_zipped_files,
                transposon_list = path_melt + "transposon_reference.list",
                removal_dir=bam_file
            log:
                log_dir + "MEIexome.log"
            shell:
                """
                mkdir {params.out_dir}/{params.sample}/melt
                ls {params[2]} | sed 's/\*//g' > {params.transposon_list}
                java -Xmx{config[MEM_GB]}G -jar {config[MELT_DIR]}MELT.jar Single -bamfile {input[0]} -h {input[1]} -t {params.transposon_list} -n {input[2]} -w {params.out_dir}/{params.sample}/melt -exome {config[MELT_CUSTOM_OPTIONS]}
                cat {params.out_dir}/{params.sample}/melt/SVA.final_comp.vcf | grep '^#' > {params.out_dir}/{params.sample}/melt/{params.sample}.header.txt
                cat {params.out_dir}/{params.sample}/melt/SVA.final_comp.vcf | grep -v '^#' > {params.out_dir}/{params.sample}/melt/{params.sample}.sva.vcf
                cat {params.out_dir}/{params.sample}/melt/LINE1.final_comp.vcf | grep -v '^#' > {params.out_dir}/{params.sample}/melt/{params.sample}.line1.vcf
                cat {params.out_dir}/{params.sample}/melt/ALU.final_comp.vcf | grep -v '^#' > {params.out_dir}/{params.sample}/melt/{params.sample}.alu.vcf
                cat {params.out_dir}/{params.sample}/melt/{params.sample}.header.txt {params.out_dir}/{params.sample}/melt/{params.sample}.sva.vcf {params.out_dir}/{params.sample}/melt/{params.sample}.line1.vcf {params.out_dir}/{params.sample}/melt/{params.sample}.alu.vcf | perl scripts/vcf-sort.pl -c | bgzip -c > {output.MEI_file}
                tabix -p vcf {output.MEI_file}
                rm -r {params.out_dir}/{params.sample}/melt
                rm {params.removal_dir}.disc {params.removal_dir}.disc.bai {params.removal_dir}.fq
                """

    if exome == "false":
        rule MEI:
            input:
                bam_file,
                path_reference,
                melt_bed
            output:
                MEI_file = results_dir + "{sample}/{sample}_MEI.vcf.gz",
                MEI_file_index = results_dir + "{sample}/{sample}_MEI.vcf.gz.tbi"
            conda:
                "envs/MEI.yaml"
            params:
                out_dir = results_dir,
                sample="{sample}",
                melt_zipped_files,
                transposon_list = path_melt + "transposon_reference.list",
                removal_dir=bam_file
            log:
                log_dir + "MEI.log"
            shell:
                """
                mkdir {params.out_dir}/{params.sample}/melt
                ls {params[2]} | sed 's/\*//g' > {params.transposon_list}
                java -Xmx{config[MEM_GB]}G -jar {config[MELT_DIR]}MELT.jar Single -bamfile {input[0]} -h {input[1]} -t {params.transposon_list} -n {input[2]} -w {params.out_dir}/{params.sample}/melt {config[MELT_CUSTOM_OPTIONS]}
                cat {params.out_dir}/{params.sample}/melt/SVA.final_comp.vcf | grep '^#' > {params.out_dir}/{params.sample}/melt/{params.sample}.header.txt
                cat {params.out_dir}/{params.sample}/melt/SVA.final_comp.vcf | grep -v '^#' > {params.out_dir}/{params.sample}/melt/{params.sample}.sva.vcf
                cat {params.out_dir}/{params.sample}/melt/LINE1.final_comp.vcf | grep -v '^#' > {params.out_dir}/{params.sample}/melt/{params.sample}.line1.vcf
                cat {params.out_dir}/{params.sample}/melt/ALU.final_comp.vcf | grep -v '^#' > {params.out_dir}/{params.sample}/melt/{params.sample}.alu.vcf                    cat {params.out_dir}/{params.sample}/melt/{params.sample}.header.txt {params.out_dir}/{params.sample}/melt/{params.sample}.sva.vcf {params.out_dir}/{params.sample}/melt/{params.sample}.line1.vcf {params.out_dir}/{params.sample}/melt/{params.sample}.alu.vcf | perl scripts/vcf-sort.pl -c | bgzip -c > {output.MEI_file}
                tabix -p vcf {output.MEI_file}
                rm -r {params.out_dir}/{params.sample}/melt
                rm {params.removal_dir}.disc {params.removal_dir}.disc.bai {params.removal_dir}.fq
                """

    if alsgenescanner == "true":
        if annotation == "true":
            rule MEIannotation:
                input:
                    rules.MEI.output.MEI_file,
                output:
                    MEI_annotation_file = results_dir + "{sample}/{sample}_MEI_annotated.tsv"
                conda:
                    "envs/SVannotation.yaml"
                params:
                    out_dir = results_dir,
                    sample = "{sample}",
                    genomebuild = annotsv_ref_version
                log:
                    log_dir + "MEIannotation.log"
                shell:
                    """
                    cpan YAML::XS
                    cpan Sort::Key::Natural
                    export ANNOTSV={config[ANNOTSV_DIR]}
                    {config[ANNOTSV_DIR]}bin/AnnotSV -annotationsDir {config[ANNOTSV_DIR]}/share/AnnotSV -SVinputFile {input[0]} -genomeBuild {params.genomebuild} -candidateGenesFiltering yes -candidateGenesFile resources/alsgenescanner/list_genes_all.txt outputDir {params.out_dir}{params.sample} -outputFile {params.sample}_MEI_annotated -SVminSize 30 {config[ANNOTSV_CUSTOM_OPTIONS]}
                    """
    else:
        if annotation == "true":
            rule MEIannotation:
                input:
                    rules.MEI.output.MEI_file
                output:
                    MEI_annotation_file = results_dir + "{sample}/{sample}_MEI_annotated.tsv"
                conda:
                    "envs/SVannotation.yaml"
                params:
                    out_dir = results_dir,
                    sample = "{sample}",
                    genomebuild = annotsv_ref_version
                log:
                    log_dir + "MEIannotation.log"
                shell:
                    """
                    cpan YAML::XS
                    cpan Sort::Key::Natural
                    export ANNOTSV={config[ANNOTSV_DIR]}
                    {config[ANNOTSV_DIR]}bin/AnnotSV -annotationsDir {config[ANNOTSV_DIR]}/share/AnnotSV -SVinputFile {input[0]} -genomeBuild {params.genomebuild} outputDir {params.out_dir}{params.sample} -outputFile {params.sample}_MEI_annotated -SVminSize 30 {config[ANNOTSV_CUSTOM_OPTIONS]}
                    """

if SV == "true" and MEI == "true":
    rule SVandMEImerging:
        input:
            rules.SV.output.merged_SV,
            rules.MEI.output.MEI_file
        output:
            merged_SV_MEI = results_dir + "{sample}/{sample}_SV_MEI.merged.vcf.gz",
            merged_SV_MEI_index = results_dir + "{sample}/{sample}_SV_MEI.merged.vcf.gz.tbi"
        conda:
            "envs/variantcalling.yaml"
        params:
            merged_dir = results_dir + "{sample}/merging"
            out_dir = results_dir,
            sample = "{sample}",
            SV_vcf = results_dir + "{sample}/merging/{sample}_merged_SV.vcf",
            MEI_vcf = results_dir + "{sample}/merging/{sample}_MEI.vcf"
        log:
            log_dir + "SVandMEImerging.log"
        shell:
            """
            mkdir {params.merged_dir}
            bgzip -d {input[0]} > {params.SV_vcf}
            bgzip -d {input[1]} > {params.MEI_vcf}
            ls {params.merged_dir}*.vcf > {params.out_dir}{params.sample}survivor_sample_files
            {config[SURVIVOR_DIR]}SURVIVOR merge {params.out_dir}{params.sample}survivor_sample_files 1000 1 1 1 0 30 {params.out_dir}{params.sample}/{params.sample}_SV_MEI.merged.vcf
            perl scripts/vcf-sort.pl {params.out_dir}{params.sample}/{params.sample}_SV_MEI.merged.vcf | bgzip -c > {output.merged_SV_MEI}
            tabix -p vcf {output.merged_SV_MEI}
            rm {params.out_dir}{params.sample}/survivor_sample_files
            """

    if alsgenescanner == "true":
        if annotation == "true":
            rule SVMEIannotation:
                input:
                    rules.SVandMEImerging.output.merged_SV_MEI,
                output:
                    SV_MEI_annotation_file = results_dir + "{sample}/{sample}_SV_MEI_annotated.tsv"
                conda:
                    "envs/SVannotation.yaml"
                params:
                    out_dir = results_dir,
                    sample = "{sample}",
                    genomebuild = annotsv_ref_version
                log:
                    log_dir + "SVandMEIannotation.log"
                shell:
                    """
                    cpan YAML::XS
                    cpan Sort::Key::Natural
                    export ANNOTSV={config[ANNOTSV_DIR]}
                    {config[ANNOTSV_DIR]}bin/AnnotSV -annotationsDir {config[ANNOTSV_DIR]}/share/AnnotSV -SVinputFile {input[0]} -genomeBuild {params.genomebuild} -candidateGenesFiltering yes -candidateGenesFile resources/alsgenescanner/list_genes_all.txt outputDir {params.out_dir}{params.sample} -outputFile {params.sample}_SV_MEI_annotated -SVminSize 30 {custom[ANNOTSV_CUSTOM_OPTIONS]}
                    """
    else:
        if annotation == "true":
            rule SVMEIannotation:
                input:
                    rules.SVandMEImerging.output.merged_SV_MEI
                output:
                    SV_MEI_annotation_file = results_dir + "{sample}/{sample}_SV_MEI_annotated.tsv"
                conda:
                    "envs/SVannotation.yaml"
                params:
                    out_dir = results_dir,
                    sample = "{sample}",
                    genomebuild = annotsv_ref_version
                log:
                    log_dir + "SVandMEIannotation.log"
                shell:
                    """
                    cpan YAML::XS
                    cpan Sort::Key::Natural
                    export ANNOTSV={config[ANNOTSV_DIR]}
                    {config[ANNOTSV_DIR]}bin/AnnotSV -annotationsDir {config[ANNOTSV_DIR]}/share/AnnotSV -SVinputFile {input[0]} -genomeBuild {params.genomebuild} outputDir {params.out_dir}{params.sample} -outputFile {params.sample}_SV_MEI_annotated -SVminSize 30 {config[ANNOTSV_CUSTOM_OPTIONS]}
                    """

if (virus or bacteria or custom_microbes) == "true":
    rule extractnonhumanreads:
        input:
            bam_file
        output:
            unaligned_fastq = results_dir + "{sample}/unaligned_reads.fastq.gz"
        conda:
            "envs/samtools.yaml"
        params:
            out_dir = results_dir,
            sample = "{sample}"
        log:
            log_dir + "extractnonhumanreads.log"
        shell:
            """
            samtools view -@ {config[NUMBER_CPU]} -hf 4 {input[0]} | samtools bam2fq -s {params.out_dir}{params.sample}/singleton_reads.fastq -@ {config[NUMBER_CPU]} - > {params.out_dir}{params.sample}/unaligned_reads.fastq
            cat {params.out_dir}{params.sample}/singleton_reads.fastq >> {params.out_dir}{params.sample}/unaligned_reads.fastq ; gzip {params.out_dir}{params.sample}/unaligned_reads.fastq
            """
#need to do the separate reports in python scripts
    if virus == "true":
        rule identifyvirus:
            input:
                rules.extractnonhumanreads.output.unaligned_fastq,
            output:
                virus_stats = results_dir + "{sample}/{sample}_virus_stats.txt",
                virus_report = reports_dir + "{sample}/{sample}_virus_report.txt"
            conda:
                "envs/alignmentfast.yaml"
            params:
                out_dir = results_dir,
                sample = "{sample}"
            log:
                log_dir + "identifyvirus.log"
            shell:
                """
                hisat2 --no-spliced-alignment -p {config[NUMBER_CPU]} -x {config[VIRUS_INDEX]} -U {input[0]} | samtools view -@ {config[NUMBER_CPU]} -hSb - | samtools sort -@ {config[NUMBER_CPU]} -T {params.out_dir}{params.sample}tempvirus.file -o {params.out_dir}{params.sample}/{params.sample}_output_virus.bam -
                samtools index -@ {config[NUMBER_CPU]} {params.out_dir}{params.sample}/{params.sample}_output_virus.bam ; samtools idxstats {params.out_dir}{params.sample}/{params.sample}_output_virus.bam > {output.virus_stats}
                python scripts/non_human_reads_report.py {output.virus_stats} {params.out_dir}{params.sample}/{params.sample}_virus_stats.list {params.out_dir}{params.sample}/{params.sample}_output_virus.bam {params.out_dir}{params.sample}/{params.sample}_virus_coverage_stats.txt {output.virus_report}
                """

    if bacteria == "true":
        rule identifybacteria:
            input:
                rules.extractnonhumanreads.output.unaligned_fastq,
            output:
                bacteria_stats = results_dir + "{sample}/{sample}_bacteria_stats.txt"
                bacteria_report = reports_dir + "{sample}/{sample}_bacteria_report.txt"
            conda:
                "envs/alignmentfast.yaml"
            params:
                out_dir = results_dir,
                sample = "{sample}"
            log:
                log_dir + "identifybacteria.log"
            shell:
                """
                hisat2 --no-spliced-alignment -p {config[NUMBER_CPU]} -x {config[BACTERIA_INDEX]} -U {input[0]} | samtools view -@ {config[NUMBER_CPU]} -hSb - | samtools sort -@ {config[NUMBER_CPU]} -T {params.out_dir}{params.sample}tempbacteria.file -o {params.out_dir}{params.sample}/{params.sample}_output_bacteria.bam -
                samtools index -@ {config[NUMBER_CPU]} {params.out_dir}{params.sample}/{params.sample}_output_bacteria.bam ; samtools idxstats {params.out_dir}{params.sample}/{params.sample}_output_bacteria.bam > {output.bacteria_stats}
                python scripts/non_human_reads_report.py {output.bacteria_stats} {params.out_dir}{params.sample}/{params.sample}_bacteria_stats.list {params.out_dir}{params.sample}/{params.sample}_output_bacteria.bam {params.out_dir}{params.sample}/{params.sample}_bacteria_coverage_stats.txt {output.bacteria_report}
                """

    if custom_microbes == "true":
        rule identifycustommicrobes:
            input:
                rules.extractnonhumanreads.output.unaligned_fastq,
            output:
                microbes_stats = results_dir + "{sample}/{sample}_microbes_stats.txt"
                microbes_report = reports_dir + "{sample}/{sample}_microbes_report.txt"
            conda:
                "envs/alignmentfast.yaml"
            params:
                out_dir = results_dir,
                sample = "{sample}"
            log:
                log_dir + "identifycustommicrobes.log"
            shell:
                """
                hisat2 --no-spliced-alignment -p {config[NUMBER_CPU]} -x {config[CUSTOM_MICROBES_INDEX]} -U {input[0]} | samtools view -@ {config[NUMBER_CPU]} -hSb - | samtools sort -@ {config[NUMBER_CPU]} -T {params.out_dir}{params.sample}tempmicrobes.file -o {params.out_dir}{params.sample}/{params.sample}_output_microbes.bam -
                samtools index -@ {config[NUMBER_CPU]} {params.out_dir}{params.sample}/{params.sample}_output_microbes.bam ; samtools idxstats {params.out_dir}{params.sample}/{params.sample}_output_microbes.bam > {output.microbes_stats}
                python scripts/non_human_reads_report.py {output.microbes_stats} {params.out_dir}{params.sample}/{params.sample}_microbes_stats.list {params.out_dir}{params.sample}/{params.sample}_output_microbes.bam {params.out_dir}{params.sample}/{params.sample}_microbes_coverage_stats.txt {output.microbes_report}
                """

if alignment_report == "true":
    rule alignmentreport:
        input:
            bam_file
        output:
            flagstat = reports_dir + "{sample}/{sample}_alignment_flagstat.txt",
            stats = reports_dir + "{sample}/{sample}_alignment_stats.txt"
        conda:
            "envs/simplereports.yaml"
        log:
            log_dir + "alignmentreport.log"
        shell:
            """
            samtools flagstat -@ {config[NUMBER_CPU]} {input[0]} > {output.flagstat}
            samtools stats -@ {config[NUMBER_CPU]} {input[0]} > {output.stats}
            """

if sequencing_report == "true" and format == "fastq":
    rule sequencingreport:
        input:
            input_file,
            input_file2
        output:
            fastqc_report = reports_dir + "{sample}/{sample}_sequencing_report.txt"
        conda:
            "envs/simplereports.yaml"
        params:
            out_dir = reports_dir
        log:
            log_dir + "sequencingreport.log"
        shell:
            """
            fastqc -o {params.out_dir} -f fastq -t {config[NUMBER_CPU]} {input[0]} {input[1]}
            """

if calls_report == "true" and variantcalling == "true":
    rule callsreport:
        input:
            variant_results_file = rules.variantfilter.output.variant_results_file_filtered if filter_string == "true" else rules.variantcalling.output.variant_results_file_filtered
        output:
            vcfstats = reports_dir + "{sample}/{sample}_calls_vcfstats.txt"
        conda:
            "envs/simplereports.yaml"
        params:
            out_dir = reports_dir,
            sample = "{sample}"
        log:
            log_dir + "callsreport.log"
        shell:
            """
            bcftools stats --threads {config[NUMBER_CPU]} {input.variant_results_file} > {output.vcfstats}
            """

if (alignment_report or calls_report or sequencing_report) == "true":
    rule multireport:
        input:
            reports_dir
        conda:
            "envs/simplereports.yaml"
        log:
            log_dir + "multireport.log"
        shell:
            """
            multiqc -o {input[0]} {input[0]}
            """

if results_report == "true":
    if annotation == "true":
        if variantcalling == "true":
            rule variantannotationreport:
                input:
                    rules.variantannotation.output.annotated_variant_results_file,
                    path_gene_list
                output:
                    SNPindel_annotation_report = reports_dir + "{sample}/{sample}_annovar_SNPindel.txt"
                    variant_annotation_file_zipped = results_dir + "{sample}/{sample}_SNPindel_annotated.vcf.gz",
                    variant_annotation_file_zipped_index = results_dir + "{sample}/{sample}_SNPindel_annotated.vcf.gz.tbi"
                conda:
                    "envs/AGS.yaml"
                params:
                    out_dir = reports_dir,
                    results = results_dir,
                    sample = "{sample}"
                log:
                    log_dir + "variantannotationreport.log"
                shell:
                    """
                    gzip -d {input[0]}
                    python scripts/simpleannotationreport.py {params.results}{params.sample}/{params.sample}_SNPindel_annotated.vcf {input[1]} {output.SNPindel_annotation_report}
                    bgzip {params.results}{params.sample}/{params.sample}_SNPindel_annotated.vcf
                    tabix -p vcf {output.variant_annotation_file_zipped}
                    """

        if expansion == "true":
            rule expansionannotationreport:
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
                    out_dir = reports_dir,
                    results = results_dir,
                    sample = "{sample}"
                log:
                    log_dir + "expansionannotationreport.log"
                shell:
                    """
                    gzip -d {input[0]}
                    python scripts/simpleannotationreport.py {params.results}{params.sample}/{params.sample}_expansions_annotated.vcf {input[1]} {output.expansion_annotation_report}
                    bgzip {params.results}{params.sample}/{params.sample}_expansions_annotated.vcf
                    tabix -p vcf {output.expansion_annotation_file_zipped}
                    """

        if genotypeSTR == "true":
            rule STRannotationreport:
                input:
                    rules.STRannotation.output.annotated_EHDN_expansion_file,
                    path_gene_list
                output:
                    STR_annotation_report = reports_dir + "{sample}/{sample}_annovar_STR.txt",
                    STR_annotation_file_zipped = results_dir + "{sample}/{sample}_EHDNexpansions_annotated.vcf.gz",
                    STR_annotation_file_zipped_index = results_dir + "{sample}/{sample}_EHDNexpansions_annotated.vcf.gz.tbi"
                conda:
                    "envs/runAGS.yaml"
                params:
                    out_dir = reports_dir,
                    results = results_dir,
                    sample = "{sample}"
                log:
                    log_dir + "STRannotationreport.log"
                shell:
                    """
                    gzip -d {input[0]}
                    python scripts/simpleannotationreport.py {params.results}{params.sample}/{params.sample}_STR_annotated.vcf {input[1]} {output.STR_annotation_report}
                    bgzip {params.results}{params.sample}/{params.sample}_STR_annotated.vcf
                    tabix -p vcf {output.STR_annotation_file_zipped}
                    """

        if (SV or MEI) == "true":
            if os.path.isfile("%s%s/%s_SV_MEI_merged.vcf.gz" % (results_dir, sample, sample)) == True:
                rule SVMEIreport:
                    input:
                        rules.SVMEIannotation.output.SV_MEI_annotation_file
                    output:
                        SV_MEI_report = reports_dir + "{sample}/{sample}_SV_MEI_annotated.html"
                    conda:
                        "envs/SVannotation.yaml"
                    params:
                        out_dir = reports_dir,
                        sample = "{sample}",
                        ref_version = annovar_ref_version
                    log:
                        log_dir + "SVMEIreport.log"
                    shell:
                        """
                        cpan YAML::XS
                        cpan Sort::Key::Natural
                        export ANNOTSV={config[ANNOTSV_DIR]}
                        perl {config[KNOTANNOTSV_DIR]}knotAnnotSV.pl --configFile {config[KNOTANNOTSV_DIR]}config_AnnotSV.yaml --annotSVfile {input[0]} --outDir {params.out_dir} --genomeBuild {params.ref_version}
                        """

            else:
                if SV:
                    rule SVreport:
                        input:
                            rules.SVannotation.output.SV_annotation_file
                        output:
                            SV_report = reports_dir + "{sample}/{sample}_SV_annotated.html"
                        conda:
                            "envs/SVannotation.yaml"
                        params:
                            out_dir = reports_dir,
                            sample = "{sample}",
                            ref_version = annovar_ref_version
                        log:
                            log_dir + "SVreport.log"
                        shell:
                            """
                            cpan YAML::XS
                            cpan Sort::Key::Natural
                            export ANNOTSV={config[ANNOTSV_DIR]}
                            perl {config[KNOTANNOTSV_DIR]}knotAnnotSV.pl --configFile {config[KNOTANNOTSV_DIR]}config_AnnotSV.yaml --annotSVfile {input[0]} --outDir {params.out_dir} --genomeBuild {params.ref_version}
                            """

                if MEI:
                    rule MEIreport:
                        input:
                            rules.MEIannotation.output.MEI_annotation_file
                        output:
                            MEI_report = reports_dir + "{sample}/{sample}_MEI_annotated.html"
                        conda:
                            "envs/MEIannotation.yaml"
                        params:
                            out_dir = reports_dir,
                            sample = "{sample}",
                            ref_version = annovar_ref_version
                        log:
                            log_dir + "MEIreport.log"
                        shell:
                            """
                            cpan YAML::XS
                            cpan Sort::Key::Natural
                            export ANNOTSV={config[ANNOTSV_DIR]}
                            perl {config[KNOTANNOTSV_DIR]}knotAnnotSV.pl --configFile {config[KNOTANNOTSV_DIR]}config_AnnotSV.yaml --annotSVfile {input[0]} --outDir {params.out_dir} --genomeBuild {params.ref_version}
                            """

#slap in concise report and write a script
if (results_report and annotation and variantcalling and SV and MEI) == "true":
    rule conciseresultsreport:
        input:
            rules.variantannotationreport.output.SNPindel_annotation_report,
            rules.SVMEIannotation.output.SV_MEI_annotation_file
        output:
            temp_SV_MEI_variants = reports_dir + "{sample}/{sample}_temp_SV_MEI_variants.tsv",
            temp_SNVindel_variants = reports_dir + "{sample}/{sample}_temp_SNVindel_variants.tsv",
            concise_report = reports_dir + "{sample}/{sample}_all_variants.tsv"
        params:
            report_header = "scripts/all_variants_report_header.txt"
        log:
            log_dir + "conciseresultsreport.log"
        shell:
            """
            python scripts/concisereportSVandMEI.py {input[1]} {output.temp_SV_MEI_variants} {input[0]} {output.temp_SNVindel_variants} {params.report_header} {output.concise_report}
            """

    if expansion or genotypeSTR == "true":
        if expansion == "true":
            rule addexpansiontoreport:
                input:
                    rules.conciseresultsreport.output.concise_report,
                    rules.expansionannotationreport.output.expansion_annotation_report
                output:
                    temp_expansion_variants = reports_dir + "{sample}/{sample}_temp_expansion_variants.tsv"
                log:
                    log_dir + "conciseresultsreportexpansion.log"
                shell:
                    """
                    python concisereportexpansion.py {input[1]} {output.temp_expansion_variants} {input[0]}
                    """
        if genotypeSTR == "true":
            rule addSTRtoreport:
                input:
                    rules.conciseresultsreport.output.concise_report,
                    rules.STRannotationreport.output.STR_annotation_report
                output:
                    temp_STR_variants = reports_dir + "{sample}/{sample}_temp_STR_variants.tsv"
                log:
                    log_dir + "conciseresultsreportSTR.log"
                shell:
                    """
                    python concisereportSTR.py {input[1]} {output.temp_STR_variants} {input[0]}
                    """

if (results_report and annotation and variantcalling and SV) == "true":
    rule conciseresultsreport:
        input:
            rules.variantannotationreport.output.SNPindel_annotation_report,
            rules.SVannotation.output.SV_annotation_file
        output:
            temp_SV_variants = reports_dir + "{sample}/{sample}_temp_SV_variants.tsv",
            temp_SNVindel_variants = reports_dir + "{sample}/{sample}_temp_SNVindel_variants.tsv",
            concise_report = reports_dir + "{sample}/{sample}_all_variants.tsv"
        params:
            report_header = "scripts/all_variants_report_header.txt"
        log:
            log_dir + "conciseresultsreport.log"
        shell:
            """
            python scripts/concisereportSVorMEI.py {input[1]} {output.temp_SV_variants} {input[0]} {output.temp_SNVindel_variants} {params.report_header} {output.concise_report}
            """

    if expansion or genotypeSTR == "true":
        if expansion == "true":
            rule addexpansiontoreport:
                input:
                    rules.conciseresultsreport.output.concise_report,
                    rules.expansionannotationreport.output.expansion_annotation_report
                output:
                    temp_expansion_variants = reports_dir + "{sample}/{sample}_temp_expansion_variants.tsv"
                log:
                    log_dir + "conciseresultsreportexpansion.log"
                shell:
                    """
                    python concisereportexpansion.py {input[1]} {output.temp_expansion_variants} {input[0]}
                    """
        if genotypeSTR == "true":
            rule addSTRtoreport:
                input:
                    rules.conciseresultsreport.output.concise_report,
                    rules.STRannotationreport.output.STR_annotation_report
                output:
                    temp_STR_variants = reports_dir + "{sample}/{sample}_temp_STR_variants.tsv"
                log:
                    log_dir + "conciseresultsreportSTR.log"
                shell:
                    """
                    python concisereportSTR.py {input[1]} {output.temp_STR_variants} {input[0]}
                    """

if (results_report and annotation and variantcalling and MEI) == "true":
    rule conciseresultsreport:
        input:
            rules.variantannotationreport.output.SNPindel_annotation_report,
            rules.MEIannotation.output.MEI_annotation_file
        output:
            temp_MEI_variants = reports_dir + "{sample}/{sample}_temp_MEI_variants.tsv",
            temp_SNVindel_variants = reports_dir + "{sample}/{sample}_temp_SNVindel_variants.tsv",
            concise_report = reports_dir + "{sample}/{sample}_all_variants.tsv"
        params:
            report_header = "scripts/all_variants_report_header.txt"
        log:
            log_dir + "conciseresultsreport.log"
        shell:
            """
            python scripts/concisereportSVorMEI.py {input[1]} {output.temp_MEI_variants} {input[0]} {output.temp_SNVindel_variants} {params.report_header} {output.concise_report}
            """

    if expansion or genotypeSTR == "true":
        if expansion == "true":
            rule addexpansiontoreport:
                input:
                    rules.conciseresultsreport.output.concise_report,
                    rules.expansionannotationreport.output.expansion_annotation_report
                output:
                    temp_expansion_variants = reports_dir + "{sample}/{sample}_temp_expansion_variants.tsv"
                log:
                    log_dir + "conciseresultsreportexpansion.log"
                shell:
                    """
                    python concisereportexpansion.py {input[1]} {output.temp_expansion_variants} {input[0]}
                    """
        if genotypeSTR == "true":
            rule addSTRtoreport:
                input:
                    rules.conciseresultsreport.output.concise_report,
                    rules.STRannotationreport.output.STR_annotation_report
                output:
                    temp_STR_variants = reports_dir + "{sample}/{sample}_temp_STR_variants.tsv"
                log:
                    log_dir + "conciseresultsreportSTR.log"
                shell:
                    """
                    python concisereportSTR.py {input[1]} {output.temp_STR_variants} {input[0]}
                    """

if alsgenescanner == "true":
    rule runAGS:
        input:
            rules.variantannotation.output.annotated_variant_results_text
        output:
            alsgenescanner_all = reports_dir + "{sample}/alsgenescanner/{sample}_alsgenescanner_all.txt",
            alsgenescanner_alsod = reports_dir + "{sample}/alsgenescanner/{sample}_alsgenescanner_alsod.txt",
            alsgenescanner_clinvar = reports_dir + "{sample}/alsgenescanner/{sample}_alsgenescanner_clinvar.txt",
            alsgenescanner_manual_review = reports_dir + "{sample}/alsgenescanner/{sample}_alsgenescanner_manual_review.txt",
            alsgenescanner_ranked = reports_dir + "{sample}/alsgenescanner/{sample}_alsgenescanner_all_ranked.txt"
        conda:
            "envs/AGS.yaml"
        params:
            out_dir = reports_dir,
            alsod_list = "resources/alsgenescanner/list_genes_alsod.txt",
            clinvar_list = "resources/alsgenescanner/list_genes_clinvar.txt",
            review_list = "resources/alsgenescanner/list_genes_manual_review.txt",
            sample = "{sample}"
        shell:
            """
            python scripts/alsgenescanner.py {input[0]} {output.alsgenescanner_all}
            python scripts/alsgenescannerreport.py {output.alsgenescanner_all} {output.alsgenescanner_alsod} {params.alsod_list} {output.alsgenescanner_clinvar} {params.clinvar_list} {output.alsgenescanner_manual_review} {params.review_list} {output.alsgenescanner_ranked}
            """
