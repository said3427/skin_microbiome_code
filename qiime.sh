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