# immunpeptidome prediction files selection of HLA allele for drawing histogram of IC50
## binding groups: 0-50 nM: strong binder, 50-500 nM: regular binder, 500-5000 nM: weak binder

library(ggpubr)
selection <- seq_complete[which(seq_complete$HLA=='HLA-B*52:01'),]
dim(selection)
gghistogram(selection, x = "ic50",add = "mean", rug = TRUE,color = "binding_group", fill = "binding_group",palette = c("#00AFBB", "#E7B800", "#FC4E07"),bins = 50)

# heatmap on proportion of peptide groups by HLA alleles

library(magrittr)
library(dplyr)
library(purrr)

## clean imported immunopeptidome data
uniq<-unique(sequence)
seq_complete <- sequnence[complete.cases(uniq), ]
dim(sequence)
dim(seq_complete)

## cross table and Chisq test
table <- xtabs(~HLA+binding_group, data=seq_complete)
print(table)
chisq.test(table, correct=FALSE)

## prepare frequence tables by HLA alleles
mat<-table(seq_complete$conca,seq_complete$binding_group)
mat<-data.frame(mat)
library(tidyr)
wide <- pivot_wider(mat, names_from = Var2, values_from = Freq) 
wide%>%rowwise(.)%>%mutate(total_line=sum(weak,regular,strong))->result
result%>%mutate(weak_pct=weak/total_line*100, regular_pct=regular/total_line*100,strong_pct=strong/total_line*100)->result_pct
result_pct %>% arrange(desc(strong_pct))
results<-as.data.frame(result_pct)
row.names(results)<-results$Var1
small_pct<-results[,6:8]

## heatmap on percents by row
pheatmap(small_pct, color = colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize = 10,cutree_rows=5,
cutree_col=3,,clustering_method = "ward.D2",clustering_distance_cols = "euclidean",clustering_distance_rows = "euclidean", show_colnames=T,
annotation_names_col=F)

# cross barplot between two categorial variables
library(ggplot2)
df <- data.frame(x = set$HLAABC, z = set$strong_binder)
df <- as.data.frame(with(df, prop.table(table(x, z), margin = NULL)))

plot <- ggplot(data = df, aes(x = x, y = Freq, fill = z)) + 
   geom_bar(width = 0.9, position = "fill", stat = "identity") + 
   scale_fill_brewer(palette = "Set1") + 
   scale_y_continuous(expand = c(0.01, 0), labels = 
   scales::percent_format()) + 
   xlab("HLAABC") + 
   ylab("Percent") + 
   labs(fill = "strong_binder") + 
   theme_classic(base_size = 14, base_family = "sans") + 
   theme(legend.position = "right")
