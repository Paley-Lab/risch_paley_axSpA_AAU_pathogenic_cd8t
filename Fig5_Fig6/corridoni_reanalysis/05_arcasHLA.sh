bsub -Is -M 16GB -R 'rusage[mem=16GB]' -n 8 -G compute-mgriffit -q general-interactive -a 'docker(continuumio/anaconda3:2024.10-1)' /bin/bash

export CONDA_ENVS_DIRS="/storage1/fs1/mgriffit/Active/griffithlab/adhoc/irisch/conda/envs/"
export CONDA_PKGS_DIRS="/storage1/fs1/mgriffit/Active/griffithlab/adhoc/irisch/conda/pkgs/"

conda activate arcasHLA


for sample in S21 S22 S23 S24 S33 S34; do
	outputdirectory=/storage1/fs1/paleym/Active/Isabel_Risch/MP004/arcasHLA/${sample}/
	mkdir $outputdirectory
	arcasHLA extract --single --unmapped --o $outputdirectory -t 8 -v /storage1/fs1/paleym/Active/SeqData/public/GSE148837_corridoni_2020/cellranger_v8.0.1/experiment1/${sample}/outs/per_sample_outs/${sample}/count/sample_alignments.bam
	arcasHLA genotype --single -g A,B,C,DPB1,DQB1,DRB1 -o $outputdirectory -t 8 -v $outputdirectory/sample_alignments.extracted.fq.gz
done
