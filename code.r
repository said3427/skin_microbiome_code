library(tidyr)
library(taRifx)
library(taxonomizr)

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

write.table(taxonomy,file="taxonomy.tsv",sep="\t",quote=F,row.names=F)

metadata <- data.table::fread("metadata_PRJNA277905.txt")
dataMetagenomics <- subset(
  data,sample%in%subset(metadata,library_selection=="RANDOM")$run_accession)

dataLong<-dataMetagenomics[,c("taxa_id","read_count","sample")]
dataWide<-spread(dataLong,value="read_count",key = "sample",fill=0)
dataWide<-data.frame(dataWide)
dataWide[dataWide=="9218868437227407266"]<-0
dataWide<-japply( dataWide, which(sapply(dataWide, class)!="character"), as.numeric )
colnames(dataWide)[1]<-"OTU ID"
write.table(file="otu.tsv",dataWide,row.names=F,quote=F,sep="\t")



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