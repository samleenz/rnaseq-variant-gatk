# rnaseq-variant-gatk
Variant calling + processing pipline using GATK based in Snakemake

## Usage

```
source activate gatk # activate conda enviornment
```
See [here](#anchors-in-markdown) to create the conda environment

```
snakemake --dag --configfile "XXX-config.yaml" | dot -Tpng > "XXX-worflow.png"
snakemake --configfile "XXX-config.yaml" --cores "N"
```

Additional parameters can be passed to the Haplotype caller step by way of the `hcArgs` variable in the config file. For example to restrict genotyping to specific positions on the genome we would add the following line to the config file. (ip provides padding arond each site, recommended by GATK site.)

```
hcArgs: "-L my_variants.vcf -ip 100"
```

### To do

* Test that calls work with GATK 4


### Build process for GATK 4 conda environment

This current process is very much not ideal, but the conda package for gatk4 is not standalone so is unavoidable. 

1. Clone the gatk 4 repo

Check if this version is that latest.

```
wget https://github.com/broadinstitute/gatk/releases/download/4.0.2.1/gatk-4.0.2.1.zip
unzip gatk-4.0.2.1.zip
cd gatk-4.0.2.1.zip
```

2. Make the conda environment

```
conda env create -n gatk4 -f gatkcondaenv.yml
```

3. Add Snakemake and java 1.8 to the conda environment (graphviz is for making snakemake DAGs)

```
source ~/miniconda/bin/activate gatk4
conda install -c bioconda java-jdk snakemake graphviz
```

Then update the `gatkPath` variable in `config.yaml` with the path to the `gatk` executable from the gatk downlaod.