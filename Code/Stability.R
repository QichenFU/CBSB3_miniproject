library(readxl)
library(ggplot2)
library(Seurat)
library(sctransform)
library(tidyverse)

# double mutant

data_double <- read_excel("../Raw_wt_TetR_3fk6_Repair_1_double.xlsx",sheet = "Ordered_by_mutation")

data_double$Site <- sub('^[^0-9]*(\\d+).*', '\\1', data_double$mutation)
data_double$From <- substr(data_double$mutation,1,1)
data_double$To <- sub('^.*[0-9]*([A-Z]),.*', '\\1', data_double$mutation)
#data_double$Site <- as.factor(as.numeric(data_double$Site))
data_double$Site <- as.factor(paste0(data_double$From,data_double$Site))
data_double$To <- factor(data_double$To,levels = c("D","E","H","K","R","Y","N","Q","S","T","A","F","I","L","M","V","G","P","W","C"))
labels = c("D","E","H","K","R","Y","N","Q","S","T","A","F","I","L","M","V","G","P","W","C")

g <- ggplot(data_double,aes(x=Site,y=delta_delta_G))
g <- g + geom_boxplot()
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
g <- g + scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0)
g <- g + labs(y = "FoldX ∆∆G (kcal/mol)")
g

g <- ggplot(data_double,aes(x=delta_delta_G))
g <- g + geom_density()
g <- g + labs(x = "Total energy ∆∆G (kcal/mol)")
g

temp <- aggregate(delta_delta_G ~ Site,data=data_double,mean)
top_n(temp,10)
g <- ggplot(temp, aes(x=Site,y=delta_delta_G))
g <- g + geom_point()
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
g <- g + labs(y = "Mean FoldX ∆∆G (kcal/mol)")
g

# To_AA
temp_2 <- aggregate(delta_delta_G ~ To,data=data_double,mean)
g <- ggplot(temp_2, aes(x=To,y=delta_delta_G))
g <- g + geom_point()
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
g <- g + labs(y = "Mean FoldX ∆∆G (kcal/mol)")
g

g <- ggplot(data_double, aes(x=To,y=delta_delta_G))
g <- g + geom_boxplot()
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
g <- g + labs(y = "FoldX ∆∆G (kcal/mol)")
g

domain <- c()
for (i in 1:nrow(data_double)){
  if (as.numeric(data_double[i,"Site"]) <= 50){
    domain <- c(domain,"DNA")
  }
  if (as.numeric(data_double[i,"Site"]) > 50){
    domain <- c(domain,"Effector")
  }
}
data_double$domain <- domain

library(ggpubr)
g <- ggplot(data_double,aes(x=domain,y=delta_delta_G,group=domain))
g <- g + geom_boxplot()
g <- g + geom_signif(comparisons = list(c("DNA","Effector"))
                    ,map_signif_level = F,test = wilcox.test)
g <- g + labs(y="Total energy ∆∆G (kcal/mol)")
g

for (i in 1:nrow(data_double)){
  if (data_double[i,"delta_delta_G"]>=10){
    data_double[i,"delta_delta_G"] <- 10
  }
  if (data_double[i,"delta_delta_G"]<=-10){
    data_double[i,"delta_delta_G"] <- -10
  }
}

# is_outlier(data_double$delta_delta_G,data_single_A$delta_delta_G,th=10)
data_double$Site <- paste0(data_double$From,as.factor(as.numeric(data_double$Site)))
g <- ggplot(data_double, aes(x=Site,y=To,fill=delta_delta_G))
g <- g + geom_raster()
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
g <- g + scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0)
g <- g + labs(y = "Mutant aa",fill="Total energy ∆∆G (kcal/mol)",title="C")
g

# https://www.it1352.com/793543.html

g <- ggplot(data_double, aes(x=Site,y=To,fill=Van_der_Waals_clashes))
g <- g + geom_raster()
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
g

g <- ggplot(data_double,aes(x=delta_delta_G,y=Van_der_Waals_clashes))
g <- g + geom_point()
g <- g + labs(x = "Total energy ∆∆G (kcal/mol)",title="A")
g

# single mutant on chain A

data_single_A <- read_excel("../Raw_wt_TetR_3fk6_Repair_1_single_A.xlsx",sheet = "Ordered_by_mutation")

data_single_A$Site <- sub('^[^0-9]*(\\d+).*', '\\1', data_single_A$mutation)
data_single_A$From <- substr(data_single_A$mutation,1,1)
data_single_A$To <- sub('^.*[0-9]*([A-Z]);', '\\1', data_single_A$mutation)
data_single_A$Site <- as.factor(as.numeric(data_single_A$Site))

for (i in 1:nrow(data_single_A)){
  if (data_single_A[i,"delta_delta_G"]>=10){
    data_single_A[i,"delta_delta_G"] <- 10
  }
  if (data_single_A[i,"delta_delta_G"]<=-10){
    data_single_A[i,"delta_delta_G"] <- -10
  }
}

g <- ggplot(data_single_A, aes(x=Site,y=To,fill=delta_delta_G))
g <- g + geom_raster()
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
g <- g + scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0)
g <- g + labs(y = "Mutant aa")
g



g <- ggplot(data_single_A, aes(x=Site,y=To,fill=Van_der_Waals_clashes))
g <- g + geom_raster()
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
g

g <- ggplot(data_single_A,aes(x=delta_delta_G,y=Van_der_Waals_clashes))
g <- g + geom_point()
g

# single mutant on chain B

data_single_B <- read_excel("../Raw_wt_TetR_3fk6_Repair_1_single_B.xlsx",sheet = "Ordered_by_mutation")

data_single_B$Site <- sub('^[^0-9]*(\\d+).*', '\\1', data_single_B$mutation)
data_single_B$From <- substr(data_single_B$mutation,1,1)
data_single_B$To <- sub('^.*[0-9]*([A-Z]);', '\\1', data_single_B$mutation)
data_single_B$Site <- as.factor(as.numeric(data_single_B$Site))

for (i in 1:nrow(data_single_B)){
  if (data_single_B[i,"delta_delta_G"]>=10){
    data_single_B[i,"delta_delta_G"] <- 10
  }
  if (data_single_B[i,"delta_delta_G"]<=-10){
    data_single_B[i,"delta_delta_G"] <- -10
  }
}

g <- ggplot(data_single_B, aes(x=Site,y=To,fill=delta_delta_G))
g <- g + geom_raster()
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
g <- g + scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0)
g <- g + labs(y = "Mutant aa")
g

g <- ggplot(data_single_B, aes(x=Site,y=To,fill=Van_der_Waals_clashes))
g <- g + geom_raster()
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
g

g <- ggplot(data_single_B,aes(x=delta_delta_G,y=Van_der_Waals_clashes))
g <- g + geom_point()
g

g <- ggplot(epistatic_data,aes(x=Site,y=To,fill=sum_AB_delta_delta_G))
g <- g + geom_raster()
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
g <- g + scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0)
g <- g + labs(y = "Mutant aa",fill="Heterzygous vs Homozygous 
Total energy ∆∆G (kcal/mol)")
g


g <- ggplot(epistatic_data, aes(x=Site,y=To,fill=delta_delta_G))
g <- g + geom_point()
g <- g + geom_label(label=epistatic_data$mutation)
g

g <- ggplot(epistatic_data,aes(x=B_delta_delta_G,y=AB_delta_delta_G))
g <- g + geom_point()
g <- g + geom_label(label=epistatic_data$mutation)
g

g <- ggplot(epistatic_data,aes(x=A_delta_delta_G,y=B_delta_delta_G))
g <- g + geom_point()
g <- g + geom_label(label=epistatic_data$mutation)
g <- g + labs(x="Chain A Total energy ∆∆G (kcal/mol)",y="Chain B Total energy ∆∆G (kcal/mol)")
g

g <- ggplot(epistatic_data,aes(x=sum,y=AB_delta_delta_G))
g <- g + geom_point()
g <- g + labs(x="Chain A Total energy ∆∆G (kcal/mol)",y="Chain B Total energy ∆∆G (kcal/mol)")
g
# TetR-DNA complex

data_DNA <- read_excel("../Raw_1qpi_Repair.xlsx",sheet = "Ordered_by_mutation")

data_DNA$Site <- sub('^[^0-9]*(\\d+).*', '\\1', data_DNA$mutation)
data_DNA$From <- substr(data_DNA$mutation,1,1)
data_DNA$To <- sub('^.*[0-9]*([A-Z]),.*', '\\1', data_DNA$mutation)
data_DNA$Site <- as.factor(as.numeric(data_DNA$Site))

for (i in 1:nrow(data_DNA)){
  if (data_DNA[i,"delta_delta_G"]>=10){
    data_DNA[i,"delta_delta_G"] <- 10
  }
  if (data_DNA[i,"delta_delta_G"]<=-10){
    data_DNA[i,"delta_delta_G"] <- -10
  }
}

g <- ggplot(data_DNA, aes(x=Site,y=To,fill=delta_delta_G))
g <- g + geom_raster()
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
g <- g + scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0)
g <- g + labs(y = "Mutant aa")
g

# correlation between stability and interaction

data_correlation <- data.frame(data_combine[,"mutation"])
data_correlation$stability <- data_double[,"delta_delta_G"]$delta_delta_G
data_correlation$interaction <- data_combine[,"delta_delta_G"]
names(data_correlation) <- c("mutation","stability","interaction")

g <- ggplot(data_correlation,aes(x=stability,y=interaction))
g <- g + geom_point()
g <- g + labs(x="Total energy ∆∆G (kcal/mol)",y="Interaction energy ∆∆G (kcal/mol)",title="")
g

