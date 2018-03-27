# 2018-03-26
# Sam Lee

configfile: "config.yaml"
from os.path import join
import re
#

# Path to GATK executable
gatk = "/home/samlee/downloads/gatk-4.0.2.1/gatk"

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
      ref = "ref.ref"
    conda:
      # "env/gatk.yaml"
    log:
      "logs/haploCaller/filterVariants.log"
    shell:
      """
      gatk -T VariantFiltration -R {params.ref} -V {input} \
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
      ref = "ref.ref"
    conda:
    #   "env/gatk.yaml"
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
      os.path.join(
          "out",
          "haplocaller",
          "all_samples.g.vcf.gz"
      )
    params:
      lst = lambda wildcards : " --variant ".join(gvcfLst),
      ref = "ref.ref",
    # conda:
      # "env/gatk.yaml"
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
      os.path.join(
          "out",
          "haploCaller",
          "{sample}",
          "{sample}.g.vcf.gz"
          )
    params:
      ref = "ref.ref",
      reg = "XYZ.vcf",
      extraArgs = "config.xxxtra" ### obv fix this
    # conda:
      # "env/gatk.yaml"
    log:
      "logs/haploCall/{sample}.log"
    shell:
      """
      {gatk} -T HaplotypeCaller -R {params.ref} -I {input.bam} \
      -dontUseSoftClippedBases -stand_call_conf 20.0 \
      -stand_emit_conf 20.0 -ERC GVCF \
      --arguments_file {params.extraArgs} -o {output} 2> {log}
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
    # conda:
      # "env/gatk.yaml"
    log:
      "logs/printBsqr/{sample}.log"
    shell:
      """
      java -d64 -jar ${gatk} -T PrintReads -R $ref -I $opdir/$bn"_processed.bam" -nct 50 -BQSR $opdir/$bn"_recal.table" -o {output} 2> {log}
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
      ref = "",
      KGIndels = "",
      millsIndels = "",
      dbSNP = ""
    # conda:
      # "env/gatk.yaml"
    log:
      "logs/bsqr/{sample}.log"
    shell:
      """
      {gatk} -T BaseRecalibrator -I {input} -R {param.ref} \
      -knownSites {params.KGIndels} -knownSites {params.millsIndels} \
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
      os.path.join(
          "out", 
          "bam", 
          "{sample}", 
          "{sample}_Aligned.sortedByCoord.dupMarked.split.out.bam"
          )
    # conda:
      # "env/gatk.yaml"
    log:
      "logs/splitNcigar/{sample}.log"
    shell:
      """
      {gatk} -T SplitNCigarReads -R {input.ref} -I {input.bam} \
      -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS 2> {log}
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
      picard MarkDuplicates I={input.bam} O={output} M=$opdir/$bn"_dup.metrics" CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT 2> {log}
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
      --outSAMtype BAM SortedByCoordinate --runThreadN {threads} 2> {log}
      """
      