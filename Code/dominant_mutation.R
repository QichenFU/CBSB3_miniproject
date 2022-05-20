library(readxl)
library(ggplot2)

epistatic_data <- read_excel("../Raw_wt_TetR_3fk6_Repair_1_single_A&B.xlsx",sheet = "Sheet1")

epistatic_data$Site <- sub('^[^0-9]*(\\d+).*', '\\1', epistatic_data$mutation)
epistatic_data$From <- substr(epistatic_data$mutation,1,1)
epistatic_data$To <- sub('^.*[0-9]*([A-Z]),.*', '\\1', epistatic_data$mutation)
epistatic_data$Site <- as.factor(as.numeric(epistatic_data$Site))

hist(epistatic_data$sum_AB_delta_delta_G)

# histogram (S8a)
g <- ggplot(epistatic_data,aes(x=sum_AB_delta_delta_G))
g <- g + geom_histogram(binwidth = 1)
g
# heatmap
g <- ggplot(epistatic_data, aes(x=Site,y=To,fill=sum_AB_delta_delta_G))
g <- g + geom_tile()
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
g <- g + scale_fill_gradientn(colours=c("blue","white","red"))
g
# cumulative distribution
g <- ggplot(epistatic_data,aes(x=sum_AB_delta_delta_G))
g <- g + stat_ecdf()
g

