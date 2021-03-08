kraken2-build --download-library archaea --db krakendb
kraken2-build --download-library bacteria --db krakendb
kraken2-build --download-library plasmid --db krakendb
kraken2-build --download-library viral --db krakendb
kraken2-build --download-library fungi --db krakendb

ls *R1* | cut -d . -f 1 |cut -d _ -f 1| sort | uniq \
    | while read sample; do
python /home/said/anaconda3/envs/qiime2-2020.2/bin/metaphlan2.py \
 qc/"$sample"_1_val_1.fq.gz,qc/"$sample"_2_val_2.fq.gz \
 --bowtie2out bowtie/"$sample".bowtie2.bz2 \
 --nproc 32 \
 --input_type fastq  > metaphlan2/"$sample".profiled_metagenome.txt;
done;


ls *R1* | cut -d . -f 1 |cut -d _ -f 1| sort | uniq \
    | while read sample; do
spades.py --meta \
 -1 qc/"$sample"_1_val_1.fq.gz  -2 qc/"$sample"_2_val_2.fq.gz \
 -t 32 \
 -o spades/"$sample" ;
done;

qiime tools import \
  --type FeatureData[Taxonomy] \
  --input-path taxonomy.tsv \
  --output-path taxonomy.qza
  
biom convert --to-hdf5 --table-type="OTU table" -i otu.tsv  -o otu.biom

qiime tools import \
  --input-path otu.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path otu.qza

qiime taxa barplot --i-table otu.qza \
  --i-taxonomy taxonomy.qza \
  --o-visualization otu.qzv \
  --m-metadata-file metadata.tsv 

qiime micom db \
  --m-meta-file embl.tsv \
  --p-folder embl_models/ \
  --p-rank species \
  --p-threads 32 \
  --o-metabolic-models embl_models.qza \
  --verbose

qiime micom build --i-abundance otu.qza \
	--i-taxonomy taxonomy.qza \
	--i-models embl_models.qza \
	--p-cutoff 0.0001 \
	--p-threads 32 \
	--o-community-models models.qza \
	--verbose

qiime micom minimal-medium --i-models models.qza \
    --p-min-growth 0.01 \
    --p-threads 32 \
    --o-medium minimal_medium.qza \
    --verbose

qiime micom grow --i-models models.qza \
    --i-medium minimal_medium.qza \
    --p-tradeoff 0.3 \
    --p-threads 32 \
    --o-results growth.qza \
    --verbose

qiime micom tradeoff --i-models models.qza \
	--i-medium minimal_medium.qza \
	--p-threads 32 \
	--o-results tradeoff.qza \
	--verbose

qiime micom plot-tradeoff --i-results tradeoff.qza \
	--o-visualization tradeoff.qzv

qiime micom plot-growth --i-results growth.qza \
	--o-visualization growth_rates.qzv

qiime micom exchanges-per-sample --i-results growth.qza \
	--o-visualization exchanges.qzv

qiime micom exchanges-per-taxon --i-results growth.qza \
	--o-visualization niche.qzv


qiime micom fit-phenotype --i-results growth.qza \
	--m-metadata-file metadata.tsv \
	--m-metadata-column status \
	--o-visualization fit.qzv
	
# Check contamination in samples	
seqtk sample data/SRR1950722_1.fastq.gz 20 