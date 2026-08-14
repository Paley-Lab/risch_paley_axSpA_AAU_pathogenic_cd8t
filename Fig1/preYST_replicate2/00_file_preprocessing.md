
## Trim and remove low-quality sequences in the fastq files with fastp
```bash
cd /storage1/fs1/paleym/Active/SeqData/MP163/20260706_preYST_Library_Rep3/
mkdir trimmed_fastq

bsub -Is -M 64GB -R 'rusage[mem=64GB]' -G compute-mgriffit -q general-interactive -a 'docker(biocontainers/fastp:v0.20.1_cv1)' /bin/bash

for Sample in *_L00*_R1_001.fastq.gz; do
  S1="${Sample%_L00*_R1_001.fastq.gz}";
  S2="${Sample%_R1_001.fastq.gz}";
  echo $S1
  fastp -i ${S2}_R1_001.fastq.gz \
    -I ${S2}_R2_001.fastq.gz \
    -o trimmed_fastq/${S2}_R1_001.fastq.gz \
    -O trimmed_fastq/${S2}_R2_001.fastq.gz \
    -l 25 \
    --trim_front1 15 --trim_front2 15 \
    --html trimmed_fastq/${S1}.fastp.html \
    2>trimmed_fastq/${S1}.fastp.log
done

```

## generate a starting csv file displaying the matched paired R1/R2 fastq file locations
```bash
cd /storage1/fs1/paleym/Active/SeqData/MP163/20260706_preYST_Library_Rep3/trimmed_fastq
for Sample in *_L00*_R1_001.fastq.gz; do
  S1="${Sample%_L00*_R1_001.fastq.gz}";
  S2="${Sample%_R1_001.fastq.gz}";
  echo "${S1},/storage1/fs1/paleym/Active/SeqData/MP163/20260706_preYST_Library_Rep3/trimmed_fastq/${S2}_R1_001.fastq.gz,/storage1/fs1/paleym/Active/SeqData/MP163/20260706_preYST_Library_Rep3/trimmed_fastq/${S2}_R2_001.fastq.gz" >> samples_trimmed.csv
done
```

## now edit THAT csv file in R, adding the actual sample names based on the GTAC metadata, and output a tsv (what the next step expects). You might need to double-check this tsv file and edit it manually so the sample names don't contain spaces or other weird symbols
```r
meta2 <- read.csv("../../SeqData/MP163/20260706_preYST_Library_Rep3/Samplemap2.csv")
tmp <- read.csv("../../SeqData/MP163/20260706_preYST_Library_Rep3/trimmed_fastq/samples_trimmed.csv", header = F)
meta <- meta2[,c(1,11)]
meta[[1]] <- sapply(strsplit(meta[[1]], "_L0"), "[", 1)
meta[[2]] <- stringr::str_replace_all(meta[[2]], " ", "_")
meta <- merge(meta, tmp, by=1)
meta <- unique(meta)
# output a file whose sample names will need to be edited manually 
write.table(meta[,c(2:4)], file = "inputs/samples_trimmed.tsv",col.names = F, row.names = F, sep="\t", quote=F)
```

## now run the 01_amplicon_count_matrix_IR_v2.py script using the following instructions for RIS1 setup (note that this only works interactively for some reason):
```bash
## make sure to copy the variants.csv file into your "inputs" folder manually 

# get the env set up
cd /storage1/fs1/paleym/Active/Isabel_Risch/MP035/
export CONDA_ENVS_DIRS="/storage1/fs1/mgriffit/Active/griffithlab/adhoc/irisch/conda/envs/"
export CONDA_PKGS_DIRS="/storage1/fs1/mgriffit/Active/griffithlab/adhoc/irisch/conda/pkgs/"
bsub -Is -M 64GB -R 'rusage[mem=64GB]' -n 1 -q general-interactive -G compute-paleym -a 'docker(continuumio/anaconda3:2024.10-1)' /bin/bash
conda activate eml-cluster-beta
# run the script
python3 01_amplicon_count_matrix_IR_v2.py
```

## now run the 02a_plot_sample_qc.py script using the following instructions for RIS1 setup (note that this only works interactively for some reason):
```bash
# get the env set up
export CONDA_ENVS_DIRS="/storage1/fs1/mgriffit/Active/griffithlab/adhoc/irisch/conda/envs/"
export CONDA_PKGS_DIRS="/storage1/fs1/mgriffit/Active/griffithlab/adhoc/irisch/conda/pkgs/"
bsub -Is -M 64GB -R 'rusage[mem=64GB]' -n 1 -q general-interactive -G compute-mgriffit -a 'docker(continuumio/anaconda3:2024.10-1)' /bin/bash
conda activate eml-cluster-beta
# run the script
python3 03a_plot_sample_qc.py
```


