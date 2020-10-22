# expression boxplot stratified by groups for two markers
library(ggpubr)
ggboxplot(set, x = "SARSCOV2",
          y = c("CALR","CANX"),
          combine = TRUE,
          color = "SARSCOV2", palette = c("black","red","green"),
          ylab = "Expression", 
          add = "jitter",                              
          add.params = list(size = 2, jitter = 0.35)  
          )


# plot of correlation analysis of markers in transcriptome assay
library(corrgram)

corrgram(data, order=TRUE, lower.panel=panel.shade,
  upper.panel=panel.pie,   #diag.panel=panel.minmax,
  main="COVID19 lung immunoproteasome") 

