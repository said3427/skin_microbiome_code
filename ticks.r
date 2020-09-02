library(qiime2R)
library(phyloseq)
library(ggplot2)
library(DESeq2)

ps<-qza_to_phyloseq(features = "table_ticks.qza",tree = "tree_ticks.qza",metadata = "manifest_ticks.tsv",taxonomy = "tax.qza")
ps.rarefied = rarefy_even_depth(ps, rngseed=1, sample.size=0.9*min(sample_sums(ps)), replace=F)

plot_bar(ps.rarefied, fill="Phylum") + facet_wrap(~Host, scales="free_x", nrow=1)
plot_bar(ps.rarefied, fill="Phylum") + facet_wrap(Sex~Host, scales="free_x", nrow=2)

ps.phylum<-tax_glom(ps.rarefied,taxrank="Phylum",NArm=F)

plot_bar(ps.phylum,fill="Phylum")
plot_bar(ps.phylum,fill="Phylum",x="Host")
plot_richness(ps.rarefied,x="State",color="Host",measures = c("Shannon"))
plot_richness(ps.rarefied,x="Host",color="Sex",measures = c("Shannon"))
plot_richness(ps.rarefied,x="Host",color="Sex",measures = c("Shannon"))+geom_boxplot()
rich = estimate_richness(ps.rarefied)
pairwise.wilcox.test(rich$Observed, sample_data(ps.rarefied)$Sex)
pairwise.wilcox.test(rich$Shannon, sample_data(ps.rarefied)$Host)

wunifrac_dist = phyloseq::distance(ps.rarefied, method="unifrac", weighted=F)
ordination = ordinate(ps.rarefied, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps.rarefied, ordination, color="Sex") + theme(aspect.ratio=1)
plot_ordination(ps.rarefied, ordination, color="Host") + theme(aspect.ratio=1)
ds = phyloseq_to_deseq2(ps, ~ Season)
ds = phyloseq_to_deseq2(ps, ~ Host)
ds = DESeq(ds)
ds = DESeq(ds)
alpha = 0.01
res = results(ds, contrast=c("Host","Bos taurus","Equus caballus"),alpha = alpha)
res = res[order(res$padj, na.last=NA), ]
res_sig = res[(res$padj < alpha), ]
res_sig = cbind(as(res_sig, "data.frame"), as(tax_table(ps)[rownames(res_sig), ], "matrix"))

ggplot(res_sig, aes(x=Genus, y=log2FoldChange, color=Phylum)) +
geom_jitter(size=3, width = 0.2) +
theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

ggplot(res_sig, aes(x=Species, y=log2FoldChange, color=Genus)) +
geom_jitter(size=3, width = 0.2) +
theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

ggplot(res_sig, aes(x=Genus, y=log2FoldChange, color=Phylum)) +
geom_jitter(size=3, width = 0.2) +
theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

write.table(res_sig,file="Ticks_stats.txt",sep=",",row.names = F,quote=F)

alpha=0.05
res = results(ds, contrast=c("Host","Bos taurus","Equus caballus"),alpha = alpha)
res = res[order(res$padj, na.last=NA), ]
res_sig = cbind(as(res, "data.frame"), as(tax_table(ps)[rownames(res), ], "matrix"))

write.table(res_sig,file="Ticks_stats.txt",sep=",",row.names = F,quote=F)

ggplot(res_sig, aes(x=Genus, y=log2FoldChange, color=Phylum)) +
geom_jitter(size=3, width = 0.2) +
theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

res_sig = res_sig[(res_sig$padj < alpha), ]

ggplot(res_sig, aes(x=Genus, y=log2FoldChange, color=Phylum)) +
geom_jitter(size=3, width = 0.2) +
theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
