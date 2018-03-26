# 2018-03-26
# Sam Lee

configfile: "config.yaml"
from os.path import join
import re
#


rule all:
    input:
      ...
    output:
      ...
    conda:
      "env/AAA.yaml"
    log:
      "logs/AAA/{sample}.log"
    shell:
      """

      """

rule haplotypeCaller:
    input:
      ...
    output:
      ...
    conda:
      "env/AAA.yaml"
    log:
      "logs/AAA/{sample}.log"
    shell:
      """

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
    os.path.join(
        "out",
        "bam",
        "{sample}",
        "{sample}_Aligned.sortedByCoord.dupMarked.split.bsqr.out.bam"
    )
    conda:
      "env/gatk.yaml"
    log:
      "logs/printBsqr/{sample}.log"
    shell:
      """
      java -d64 -jar $gatk -T PrintReads -R $ref -I $opdir/$bn"_processed.bam" -nct 50 -BQSR $opdir/$bn"_recal.table" -o {output}
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
      "out/bsqr/{sample}_recal.table"
    params:
      ref = "",
      KGIndels = "",
      millsIndels = "",
      dbSNP = ""
    conda:
      "env/gatk.yaml"
    log:
      "logs/bsqr/{sample}.log"
    shell:
      """
      gatk -T BaseRecalibrator -I {input} -R {param.ref} \
      -knownSites {params.KGIndels} -knownSites {params.millsIndels} \
      -knownSites {params.dbSNP} -o {output}
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
      os.path.join(
          "out", 
          "bam", 
          "{sample}", 
          "{sample}_Aligned.sortedByCoord.dupMarked.split.out.bam"
          )
    conda:
      "env/gatk.yaml"
    log:
      "logs/splitNcigar/{sample}.log"
    shell:
      """
      gatk -T SplitNCigarReads -R {input.ref} -I {input.bam} \
      -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
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
      os.path.join(
          "out", 
          "bam", 
          "{sample}", 
          "{sample}_Aligned.sortedByCoord.dupMarked.out.bam"
          )
    conda:
      "env/picard.yaml"
    log:
      "logs/markDuplicates/{sample}.log"
    shell:
      """
      MarkDuplicates I={input.bam} O={output} M=$opdir/$bn"_dup.metrics" CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT 
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
      samtools index -b {input}
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
      starRef = "{config.star_ref}" ##### CHECK THIS
    conda:
      "env/star.yaml"
    log:
      "logs/alignReads/{sample}.log"
    threads: 99
    shell:
      """
      STAR --twopassMode Basic --genomeDir {params.starRef} \
      --readFilesIn {input.fq1} {input.fq2} \
      --readFilesCommand zcat --outFileNamePrefix {params.prefix}  \
      --outSAMtype BAM SortedByCoordinate --runThreadN {threads}
      """
      