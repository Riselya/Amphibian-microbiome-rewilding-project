############################## QIIME2 ANALYSIS 

cd UOW_frogs

conda activate qiime2-2020.2
conda activate qiime2-amplicon-2024.5

#gzip *.fastq

##import sequence data

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path paired_end \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-paired-end.qza

  ##how many sequences were obtained per sample?


qiime demux summarize \
--i-data demux-paired-end.qza \
--o-visualization demux-paired-end.qzv

##view wtih qiime2

#he nuclear ribosomal internal transcribed spacer (ITS) fungal gene was amplified using the primers 
#ITS1F (CTTGGTCATTTAGAGGAAGTAA) and ITS2 (GCTGCGTTCTTCATCGATGC).  

## dada2 quality control, joins sequences, filters out chimeras

qiime dada2 denoise-paired \
--i-demultiplexed-seqs demux-paired-end.qza \
--p-trim-left-f 22 \
--p-trim-left-r 20 \
--p-trunc-len-f 230 \
--p-trunc-len-r 220 \
--p-n-threads 8 \
--o-representative-sequences rep-seqs-dada2.qza \
--o-table table-dada2.qza \
--o-denoising-stats stats-dada2.qza


##visualize denoising results per sample

qiime metadata tabulate \
--m-input-file stats-dada2.qza \
--o-visualization stats-dada2.qzv 

#########  summarize count table 

qiime feature-table summarize \
--i-table table-dada2.qza \
--o-visualization table-dada2.qzv \
--m-sample-metadata-file frog_metadata_2022.txt


##############################################################


qiime feature-table tabulate-seqs \
--i-data rep-seqs-dada2.qza \
--o-visualization rep-seqs-dada2.qzv


######train silva classifier for ITS1-ITS2 region
######train silva classifier for ITS1-ITS2 region
######train silva classifier for ITS1-ITS2 region

# import UNITE database for classification

qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path sh_refs_qiime_ver10_99_04.04.2024.fasta \
--output-path unite_ref_seqs.qza


qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path sh_taxonomy_qiime_ver10_99_04.04.2024.txt \
--output-path unite_taxonomy.qza

## Derelicate the database


qiime rescript dereplicate \
  --i-sequences unite_ref_seqs.qza \
  --i-taxa unite_taxonomy.qza \
  --p-mode uniq \
  --o-dereplicated-sequences unite_derep_seqs.qza \
  --o-dereplicated-taxa unite_derep_taxonomy.qza


# Create classifier

  #ITS1F - ITS2R
# ITS1 CTTGGTCATTTAGAGGAAGTAA
# ITS2 GCTGCGTTCTTCATCGATGC

# trim 

qiime feature-classifier extract-reads \
--i-sequences unite_derep_seqs.qza \
--p-f-primer CTTGGTCATTTAGAGGAAGTAA \
--p-r-primer GCTGCGTTCTTCATCGATGC \
--o-reads unite-trimmed-ITS1-ITS2.qza

# create primer-specific classifier

qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads unite-trimmed-ITS1-ITS2.qza \
--i-reference-taxonomy unite_derep_taxonomy.qza \
--o-classifier unite.IT1.ITS2.classifier.qza



############### assign taxonomy to the rerepsentative sequences


qiime feature-classifier classify-sklearn \
--i-classifier unite.IT1.ITS2.classifier.qza \
--i-reads rep-seqs-dada2.qza \
--p-n-jobs 8 \
--o-classification taxonomyUnite.qza


#visualize 

qiime metadata tabulate \
--m-input-file taxonomyUnite.qza \
--o-visualization taxonomyUnite.qzv


#### create barplot


qiime taxa barplot \
--i-table table-dada2.qza \
--i-taxonomy taxonomyUnite.qza \
--m-metadata-file frog_metadata_2022.txt \
--o-visualization taxa-bar-plots.qzv

# export ASVs

qiime tools export \
--input-path rep-seqs-dada2.qza \
--output-path exported-representative-sequences

#lets align our fasta file using mafft

mafft dna-sequences.fasta > SequencesAligned.fasta


#We next build a tree using the software FastTree:

FastTree < SequencesAligned.fasta > corroboree.ITS.tree 


#################### export files

qiime tools export \
--input-path table-dada2.qza \
--output-path exported-count-table


qiime tools export \
--input-path taxonomyUnite.qza \
--output-path exported-taxonomy





