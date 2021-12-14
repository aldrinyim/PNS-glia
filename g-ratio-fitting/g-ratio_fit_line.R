
library(ggplot2)
library(reshape2)

data <- read.csv("g-ratio-analysis/Sciatic-g-ratio_data.txt",header=T,sep="\t")
data <- read.csv("g-ratio-analysis/Tibial-g-ratio_data.txt",header=T,sep="\t")
data <- read.csv("g-ratio-analysis/human_data.txt",header=T,sep="\t")
#m.data <- melt(data)
#m.data

head(data)
others <- subset(data, Type == "Others")
pmp2 <- subset(data, Type == "Pmp2+")


ggplot(data,aes(x=AxonDia,y=G.ratio,shape=Type,color=Type)) + geom_point() + geom_smooth(method=lm)


# Creo i bin
others$Groups <- cut(x=others$AxonDia, breaks=seq(from=0, to=ceiling(max(data$AxonDia)), by = 1))
others$Groups
# Raggruppo i dati nei bin
Bygroup = tapply(others$AxonDia, others$Groups, length)
Bygroup
# Faccio il plot dei bin
barplot(height = Bygroup, xlab = "axon.size", ylab = "occurence")


# Creo i bin
pmp2$Groups <- cut(x=pmp2$AxonDia, breaks=seq(from=0, to=ceiling(max(data$AxonDia)), by = 1))
#pmp2$Groups <- cut(x=pmp2$AxonDia, breaks=seq(from=0, to=13, by = 1))
pmp2$Groups
# Raggruppo i dati nei bin
#Bygroup = tapply(pmp2$AxonDia, pmp2$Groups, sum)
Bygroup = tapply(pmp2$AxonDia, pmp2$Groups, length)
Bygroup
# Faccio il plot dei bin
barplot(height = Bygroup, xlab = "axon.size", ylab = "occurence")


ggplot(data, aes(x=Groups, y=Bygroup, group=Type, fill=factor(Type))) + geom_bar(stat="identity", width=.5, position = "dodge")

