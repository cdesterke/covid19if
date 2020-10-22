library(randomForest)

# randomForest 4.6-14 on prevalence data with group column comprising predict region group


complete_data<-input_data[complete.cases(input_data),]

set<-complete_data[,2:ncol(complete_data)]

rf <- randomForest(group ~ ., data=set, importance=TRUE, proximity=TRUE)

MDSplot(rf, set$group)

varImpPlot(rf)
