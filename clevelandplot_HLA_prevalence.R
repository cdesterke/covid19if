# create cleveland plot on HLA prevalence

library(ggpubr)


ggdotchart(set, x = "name", y = "B.52.01", color = "group",                               
           palette = c("darkorchid3","darkorange3","aquamarine4"),  						
           sorting = "descending", rotate = TRUE, dot.size = 4,                                 
           y.text.col = TRUE,ggtheme = theme_pubr())+  theme_cleveland()

## set : data input
### column x=name, names of the world regions of ethnies
### column y="B.52.01, prevalences by region for this HLA allele
### group, group of world region predicted by machine learning

