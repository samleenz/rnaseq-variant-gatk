# rnaseq-variant-gatk
Variant calling + processing pipline using GATK based in Snakemake

### To do

* Test that calls work with GATK 4


### Build process for GATK 4

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

3. Add Snakemake and java 1.8 to the conda environment

```
source ~/miniconda/bin/activate gatk4
conda install -c bioconda java-jdk snakemake
```

Then update the `gatk` variable in `Snakemake` with the path to the `gatk` executable from the gatk downlaod.