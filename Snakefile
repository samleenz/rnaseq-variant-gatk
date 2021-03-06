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
## output directory name
outDir      = config["outDir"]
## fastq file directory 
fastqDir    = config["fastqDir"]

# List of "{sample}.g.vcf.gz"
# used for rule "combineGVCFs"
gvcfLst = expand(
    os.path.join(
      outDir,
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
        outDir,
        "haplocaller",
        "all_samples.genotyped.filtered.vcf.gz"
        )   


rule filterVCF:
    input:
      os.path.join(
          outDir,
          "haplocaller",
          "all_samples.genotyped.vcf.gz"
      )
    output:
      os.path.join(
          outDir,
          "haplocaller",
          "all_samples.genotyped.filtered.vcf.gz"
      )
    params:
      ref = ref
    log:
      "logs/haploCaller/filterVariants.log"
    shell:
      """
      {gatk} VariantFiltration -R {params.ref} -V {input} \
        -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" \
        -filterName QD -filter "QD < 2.0" -o {output} 2> {log}
      """


rule genotypeGVCFs:
    input:
      os.path.join(
          outDir,
          "haplocaller",
          "all_samples.g.vcf.gz"
      )
    output:
      os.path.join(
          outDir,
          "haplocaller",
          "all_samples.genotyped.vcf.gz"
      )
    params:
      ref = ref
    log:
      "logs/haploCall/genotypeGVCFs.log"
    shell:
      """
      {gatk} GenotypeGVCFs -R {params.ref} \
        --variant {input} -o {output} 2> {log}
      """


rule combineGVCFs:
    input:
      gvcfLst
    output:
      temp(os.path.join(
          outDir,
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
      {gatk} CombineGVCFs -R {params.ref} \
        --variant {params.lst} -o {output} 2> {log}
      """


rule haplotypeCaller:
    input:
      bam = os.path.join(
          outDir,
          "bam",
          "{sample}",
          "{sample}_Aligned.sortedByCoord.dupMarked.split.bsqr.out.bam"
          )
    output:
      temp(os.path.join(
          outDir,
          "haploCaller",
          "{sample}",
          "{sample}.g.vcf.gz"
          ))
    params:
      ref = ref,
      hcArgs = hcArgs
    log:
      "logs/haploCall/{sample}.log"
    shell:
      """
      {gatk} HaplotypeCaller -R {params.ref} -I {input.bam} \
      -dontUseSoftClippedBases -stand_call_conf 20.0 \
      --output_mode EMIT_ALL_CONFIDENT_SITES -stand_emit_conf 20.0 \
      -ERC GVCF {params.hcArgs} -o {output} 2> {log}
      """


rule printBsqr:
    input:
      bam = os.path.join(
          outDir, 
          "bam", 
          "{sample}", 
          "{sample}_Aligned.sortedByCoord.dupMarked.split.out.bam"
          ),
      table = os.path.join(
          outDir,
          "bsqr/{sample}_recal.table"
          )
    output:
      temp(os.path.join(
        outDir,
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
      {gatk} PrintReads -R {params.ref} -I {input.bam} \
       -nct 50 -BQSR {input.table} -o {output} 2> {log}
      """


rule bsqr:
    input:
      os.path.join(
          outDir, 
          "bam", 
          "{sample}", 
          "{sample}_Aligned.sortedByCoord.dupMarked.split.out.bam"
          )
    output:
      os.path.join(
        outDir,
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
      {gatk} BaseRecalibrator -I {input} -R {params.ref} \
      -knownSites {params.KGsnps} -knownSites {params.millsIndels} \
      -knownSites {params.dbSNP} -o {output} 2> {log}
      """



rule splitNcigar:
    input:
      bam = os.path.join(
          outDir, 
          "bam", 
          "{sample}", 
          "{sample}_Aligned.sortedByCoord.dupMarked.out.bam"
          )
    output:
      temp(os.path.join(
          outDir, 
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
      {gatk} SplitNCigarReads -R {params.ref} -I {input.bam} \
      -O {output} -rf ReassignOneMappingQuality -RMQF 255 \
      -RMQT 60 -U ALLOW_N_CIGAR_READS 2> {log}
      """



rule markDuplicates:
    input:
      bam = os.path.join(
          outDir, 
          "bam", 
          "{sample}", 
          "{sample}_Aligned.sortedByCoord.out.bam"
          ),
      bai = os.path.join(
          outDir, 
          "bam", 
          "{sample}", 
          "{sample}_Aligned.sortedByCoord.out.bam.bai"
          )
    output:
      temp(os.path.join(
          outDir, 
          "bam", 
          "{sample}", 
          "{sample}_Aligned.sortedByCoord.dupMarked.out.bam"
          ))
    params:
      metrics = os.path.join(
          outDir,
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
          outDir, 
          "bam", 
          "{sample}", 
          "{sample}_Aligned.sortedByCoord.out.bam"
          )
    output:
      os.path.join(
          outDir, 
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
      fq1 = os.path.join(
          fastqDir,
          "{sample}_1_trimmed.fastq.gz",
          ),
      fq2 = os.path.join(
          fastqDir,
          "{sample}_2_trimmed.fastq.gz",
          ),
    output:
      os.path.join(
          outDir, 
          "bam", 
          "{sample}", 
          "{sample}_Aligned.sortedByCoord.out.bam"
          )
    params:
      prefix  = os.path.join(
          outDir,
          "star/{sample}/{sample}_"
          ),
      starRef = star_ref,
      rg = "ID:{sample} SM:{sample} PL:Illumina LB:PairedEnd"
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
      --outSAMattrRGline {params.rg} --runThreadN {threads} 2> {log}
      """
      
