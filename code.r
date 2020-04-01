library(tidyr)
library(taRifx)
library(taxonomizr)
library(data.table)


data<-c()
for(i in dir(path="slimm/")){
dataSlimm<-data.table::fread(paste("./slimm/",i,sep=""))

dataSlimm$sample<-unlist(strsplit(i,"_profile.tsv"))[1]
data<-rbind(data,dataSlimm)
}

taxonomy <- data[,c("taxa_id","linage")]
taxonomy<- taxonomy[!duplicated(taxonomy),]
taxonomy$linage <- gsub("[|]",";",taxonomy$linage)
colnames(taxonomy)<-c("Feature ID","Taxon")

taxonomy$Taxon<-gsub("k__","D_0__",taxonomy$Taxon)
taxonomy$Taxon<-gsub("p__","D_1__",taxonomy$Taxon)
taxonomy$Taxon<-gsub("c__","D_2__",taxonomy$Taxon)
taxonomy$Taxon<-gsub("o__","D_3__",taxonomy$Taxon)
taxonomy$Taxon<-gsub("f__","D_4__",taxonomy$Taxon)
taxonomy$Taxon<-gsub("g__","D_5__",taxonomy$Taxon)
taxonomy$Taxon<-gsub("s__","D_6__",taxonomy$Taxon)
#taxonomy$Taxon<-sub("D_6__.+? ","D_6__",taxonomy$Taxon)


taxonomy$`Feature ID`<-gsub('\\*',"a",taxonomy$`Feature ID`)

metadata <- data.table::fread("metadata_PRJNA277905.txt")

dataMetagenomics <- subset(
  data,sample%in%subset(metadata,library_selection=="RANDOM")$run_accession)

metadata <- metadata[,c(5,1:4,6:19)]
colnames(metadata)[1]<-"sample-id"
write.table(metadata,file="metadata.tsv",sep="\t",quote=F,row.names=F)

dataLong<-dataMetagenomics[,c("taxa_id","read_count","sample")]
dataWide<-spread(dataLong,value="read_count",key = "sample",fill=0)
dataWide<-data.frame(dataWide)
dataWide[dataWide=="9218868437227407266"]<-0
dataWide<-japply( dataWide, which(sapply(dataWide, class)!="character"), as.numeric )
colnames(dataWide)[1]<-"OTU ID"

dataFrac<-(t(t(dataWide[,-1])/colSums(dataWide[,-1],na.rm=T)))
dataFrac<-cbind(dataWide[,1],dataFrac)
dataModel<-(dataFrac[rowSums(dataFrac[,-1]>0.00)>=5,])

taxonomy<-taxonomy[taxonomy$`Feature ID` %in% dataModel$`OTU ID`,]

write.table(taxonomy,file="taxonomy.tsv",sep="\t",quote=F,row.names=F)
write.table(file="otu.tsv",dataWide,row.names=F,quote=F,sep="\t")

# biom convert --to-hdf5 --table-type="OTU table" -i otu.tsv  -o otu.biom
# qiime tools import --input-path otu.biom --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path otu.qza
# qiime tools import --type FeatureData[Taxonomy] --input-path taxonomy.tsv --output-path taxonomy.qza

modelsMetadata<-data.table::fread("embl_gems/model_list.tsv")
modelsMetadata$file_path<-sapply(strsplit(modelsMetadata$file_path,"/"),`[`,4)
modelsMetadata$file_path<-gsub(".xml.gz",".xml",modelsMetadata$file_path)
modelsMetadata$id<-gsub(".xml","",modelsMetadata$file_path)

prepareDatabase('accessionTaxa.sql')
taxaId<-modelsMetadata$taxid
modelsTaxonomy<-getTaxonomy(taxaId,'accessionTaxa.sql')
colnames(modelsTaxonomy)[1]<-"kingdom"


models<-cbind(modelsMetadata,modelsTaxonomy)

models<-cbind(modelsMetadata,modelsTaxonomy)
colnames(models)[5]<-"file"
colnames(models)[2]<-"taxa_id"
colnames(models)[3]<-"organism"

models<-models[,c("id","organism","kingdom","phylum","class","order","family","genus","species","file","taxa_id","assembly_accession","infraspecific_name")]
models<-models[-2005,]
write.table(models,file="embl.tsv",sep="\t",quote=F,row.names=F)