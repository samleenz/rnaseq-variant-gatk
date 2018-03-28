# 2018-03-26
# Sam Lee
# Snakefile for rnaseq-variant-gatk pipeline

from os.path import join
import re

# Import config file for conda options
configfile: "config.yaml"
## extract resource variables
star_ref    = config["star_ref"]
ref         = config["ref"]
KGsnps      = config["KGsnps"]
millsIndels = config["millsIndels"]
dbSNP       = config["dbSNP"]
hcArgs      = config["hcArgs"]
## Path to GATK executable
gatk        = config["gatkPath"]

# List of "{sample}.g.vcf.gz"
# used for rule "combineGVCFs"
gvcfLst = expand(
    os.path.join(
      "out",
      "haploCaller",
      "{sample}",
      "{sample}.g.vcf.gz" 
      ),
    sample=config["samples"]
    )

####################
# Rule definitions
####################

rule all:
    input:
      os.path.join(
        "out",
        "haplocaller",
        "all_samples.genotyped.filtered.vcf.gz"
        )   


rule filterVCF:
    input:
      os.path.join(
          "out",
          "haplocaller",
          "all_samples.genotyped.vcf.gz"
      )
    output:
      os.path.join(
          "out",
          "haplocaller",
          "all_samples.genotyped.filtered.vcf.gz"
      )
    params:
      ref = ref
    log:
      "logs/haploCaller/filterVariants.log"
    shell:
      """
      {gatk} -T VariantFiltration -R {params.ref} -V {input} \
        -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" \
        -filterName QD -filter "QD < 2.0" -o {output} 2> {log}
      """


rule genotypeGVCFs:
    input:
      os.path.join(
          "out",
          "haplocaller",
          "all_samples.g.vcf.gz"
      )
    output:
      os.path.join(
          "out",
          "haplocaller",
          "all_samples.genotyped.vcf.gz"
      )
    params:
      ref = ref
    log:
      "logs/genotypeGVCF/{sample}.log"
    shell:
      """
      {gatk} -T GenotypeGVCFs -R {params.ref} \
        --variant {input} -o {output} 2> {log}
      """


rule combineGVCFs:
    input:
      gvcfLst
    output:
      temp(os.path.join(
          "out",
          "haplocaller",
          "all_samples.g.vcf.gz"
      ))
    params:
      lst = lambda wildcards : " --variant ".join(gvcfLst),
      ref = ref
    log:
      "logs/haploCall/combineGVCF.log"
    shell:
      """
      {gatk} -T CombineGVCFs -R {params.ref} \
        --variant {params.lst} -o {output} 2> {log}
      """


rule haplotypeCaller:
    input:
      bam = os.path.join(
          "out",
          "bam",
          "{sample}",
          "{sample}_Aligned.sortedByCoord.dupMarked.split.bsqr.out.bam"
          )
    output:
      temp(os.path.join(
          "out",
          "haploCaller",
          "{sample}",
          "{sample}.g.vcf.gz"
          ))
    params:
      ref = ref
      hcArgs = hcArgs
    log:
      "logs/haploCall/{sample}.log"
    shell:
      """
      {gatk} -T HaplotypeCaller -R {params.ref} -I {input.bam} \
      -dontUseSoftClippedBases -stand_call_conf 20.0 \
      -stand_emit_conf 20.0 -ERC GVCF \
      {params.hcArgs} -o {output} 2> {log}
      """


rule printBsqr:
    input:
      bam = os.path.join(
          "out", 
          "bam", 
          "{sample}", 
          "{sample}_Aligned.sortedByCoord.dupMarked.split.out.bam"
          ),
      table = "out/bsqr/{sample}_recal.table"
    output:
    temp(os.path.join(
        "out",
        "bam",
        "{sample}",
        "{sample}_Aligned.sortedByCoord.dupMarked.split.bsqr.out.bam"
        ))
    params:
      ref = ref
    log:
      "logs/printBsqr/{sample}.log"
    shell:
      """
      {gatk} -T PrintReads -R {params.ref} -I {input.bam} \
       -nct 50 -BQSR {input.table} -o {output} 2> {log}
      """


rule bsqr:
    input:
      os.path.join(
          "out", 
          "bam", 
          "{sample}", 
          "{sample}_Aligned.sortedByCoord.dupMarked.split.out.bam"
          )
    output:
      os.path.join(
        "out",
        "bsqr",
        "{sample}_recal.table"
        )
    params:
      ref = ref,
      KGsnps = KGsnps,
      millsIndels = millsIndels,
      dbSNP = dbSNP
    log:
      "logs/bsqr/{sample}.log"
    shell:
      """
      {gatk} -T BaseRecalibrator -I {input} -R {param.ref} \
      -knownSites {params.KGsnps} -knownSites {params.millsIndels} \
      -knownSites {params.dbSNP} -o {output} 2> {log}
      """



rule splitNcigar:
    input:
      bam = os.path.join(
          "out", 
          "bam", 
          "{sample}", 
          "{sample}_Aligned.sortedByCoord.dupMarked.out.bam"
          )
    output:
      temp(os.path.join(
          "out", 
          "bam", 
          "{sample}", 
          "{sample}_Aligned.sortedByCoord.dupMarked.split.out.bam"
          ))
    params:
      ref = ref
    log:
      "logs/splitNcigar/{sample}.log"
    shell:
      """
      {gatk} -T SplitNCigarReads -R {params.ref} -I {input.bam} \
      -rf ReassignOneMappingQuality -RMQF 255 \
      -RMQT 60 -U ALLOW_N_CIGAR_READS 2> {log}
      """



rule markDuplicates:
    input:
      bam = os.path.join(
          "out", 
          "bam", 
          "{sample}", 
          "{sample}_Aligned.sortedByCoord.out.bam"
          ),
      bai = os.path.join(
          "out", 
          "bam", 
          "{sample}", 
          "{sample}_Aligned.sortedByCoord.out.bam.bai"
          )
    output:
      temp(os.path.join(
          "out", 
          "bam", 
          "{sample}", 
          "{sample}_Aligned.sortedByCoord.dupMarked.out.bam"
          ))
    params:
      metrics = os.path.join(
          "out",
          "picard-tools-marked-dup-metrics.txt"
          )
    conda:
      "env/picard.yaml"
    log:
      "logs/markDuplicates/{sample}.log"
    shell:
      """
      picard MarkDuplicates I={input.bam} O={output} \
      M={params.metrics} CREATE_INDEX=true \
      VALIDATION_STRINGENCY=SILENT 2> {log}
      """


rule indexReads:
    input:
      os.path.join(
          "out", 
          "bam", 
          "{sample}", 
          "{sample}_Aligned.sortedByCoord.out.bam"
          )
    output:
      os.path.join(
          "out", 
          "bam", 
          "{sample}", 
          "{sample}_Aligned.sortedByCoord.out.bam.bai"
          )
    conda:
      "env/samtools.yaml"
    log:
      "logs/indexReads/{sample}.log"
    shell:
      """
      samtools index -b {input} 2> {log}
      """


rule alignReads:
    input:
      fq1 = "{sample}_1_trimmed.fastq.gz,
      fq2 = "{sample}_2_trimmed.fastq.gz
    output:
      os.path.join(
          "out", 
          "bam", 
          "{sample}", 
          "{sample}_Aligned.sortedByCoord.out.bam"
          )
    params:
      prefix  = "star/{sample}/{sample}_",
      starRef = star_ref,
      rg = "ID:{sample} SM:{sample}"
    conda:
      "env/star.yaml"
    log:
      "logs/alignReads/{sample}.log"
    threads: 8 # this is arbitrary...
    shell:
      """
      STAR --twopassMode Basic --genomeDir {params.starRef} \
      --readFilesIn {input.fq1} {input.fq2} --readFilesCommand zcat \
      --outFileNamePrefix {params.prefix} --outSAMtype BAM SortedByCoordinate \
      SortedByCoordinate {params.rg} --runThreadN {threads} 2> {log}
      """
      