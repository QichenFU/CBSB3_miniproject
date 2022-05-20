# data cleaning

data <- read.csv("../all_interaction_3fk6_double.csv")

names(data) <- data[8,]

data_combine <- data.frame(data[9,])

for (i in 2:3648){
  data_combine <- rbind(data_combine,data.frame(data[i*10-1,]))
}

data_combine$order <- as.numeric(sub('^.*_(\\d+).pdb', '\\1', data_combine$Pdb))
data_combine <- data_combine[order(data_combine$order),]

write.csv(data_combine,"../all_interaction_cleaned_3fk6_double.csv")

##############
data <- read.csv("../all_interaction_1qpi.csv")

names(data) <- data[8,]

data_combine <- data.frame(data[9,])

for (i in 2:3902){
  data_combine <- rbind(data_combine,data.frame(data[i*10-1,]))
}

data_combine$order <- as.numeric(sub('^.*_(\\d+).pdb', '\\1', data_combine$Pdb))
data_combine <- data_combine[order(data_combine$order),]

write.csv(data_combine,"all_interaction_cleaned_1qpi_DNA.csv")

# 48, 2319 (Chain B)

# plot interaction

data_combine <- read.csv("all_interaction_cleaned_3fk6_B.csv")


data_combine$Site <- sub('^[^0-9]*(\\d+).*', '\\1', data_combine$mutation)
data_combine$From <- substr(data_combine$mutation,1,1)
data_combine$To <- sub('^.*[0-9]*([A-Z]);', '\\1', data_combine$mutation)
data_combine$Site <- as.factor(as.numeric(data_combine$Site))

for (i in 1:nrow(data_combine)){
  if (data_combine[i,"delta_delta_G"]>=10){
    data_combine[i,"delta_delta_G"] <- 10
  }
  if (data_combine[i,"delta_delta_G"]<=-10){
    data_combine[i,"delta_delta_G"] <- -10
  }
}

library(ggplot2)
g <- ggplot(data = data_combine, aes(x=Site,y=To,fill=delta_delta_G))
g <- g + geom_tile()
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
g <- g + scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0)
g

# Chain A&B
# plot interaction

data_combine <- read.csv("../all_interaction_cleaned_3fk6_double.csv")


data_combine$Site <- sub('^[^0-9]*(\\d+).*', '\\1', data_combine$mutation)
data_combine$From <- substr(data_combine$mutation,1,1)
data_combine$To <- sub('^.*[0-9]*([A-Z]);', '\\1', data_combine$mutation)
data_combine$Site <- as.factor(as.numeric(data_combine$Site))
data_combine$To <- factor(data_combine$To,levels = c("D","E","H","K","R","Y","N","Q","S","T","A","F","I","L","M","V","G","P","W","C"))

for (i in 1:nrow(data_combine)){
  if (data_combine[i,"delta_delta_G"]>=10){
    data_combine[i,"delta_delta_G"] <- 10
  }
  if (data_combine[i,"delta_delta_G"]<=-10){
    data_combine[i,"delta_delta_G"] <- -10
  }
}

data_combine$order <- as.numeric(sub('^.*_(\\d+).pdb', '\\1', data_combine$Pdb))
data_combine <- data_combine[order(data_combine$order),]

library(ggplot2)
g <- ggplot(data = data_combine, aes(x=Site,y=To,fill=delta_delta_G))
g <- g + geom_tile()
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
g <- g + labs(y = "Mutant aa",fill="Interaction energy ∆∆G (kcal/mol)")
g <- g + scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0)
g

#########
data_combine <- read.csv("../all_interaction_cleaned_1qpi_DNA.csv")


data_combine$Site <- sub('^[^0-9]*(\\d+).*', '\\1', data_combine$mutation)
data_combine$From <- substr(data_combine$mutation,1,1)
data_combine$To <- sub('^.*[0-9]*([A-Z]);', '\\1', data_combine$mutation)
data_combine$Site <- as.factor(as.numeric(data_combine$Site))

for (i in 1:nrow(data_combine)){
  if (data_combine[i,"delta_delta_G"]>=10){
    data_combine[i,"delta_delta_G"] <- 10
  }
  if (data_combine[i,"delta_delta_G"]<=-10){
    data_combine[i,"delta_delta_G"] <- -10
  }
}

library(ggplot2)
g <- ggplot(data = data_combine, aes(x=Site,y=To,fill=delta_delta_G))
g <- g + geom_tile()
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
g <- g + scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0)
g <- g + labs(y = "Mutant aa",fill="Interaction energy ∆∆G (kcal/mol)")
g

