library(readxl)
library(ggplot2)

#O. Scholz et al.

data_paper <- read_excel("../data_from_paper.xlsx",sheet = "O. Scholz et al.")

delta_delta_G <- c()
DNA_binding <- c()
for (i in 1:nrow(data_paper)){
  for (j in 1:nrow(data_combine)){
    if (grepl(data_paper[i,"mutation"]$mutation,data_combine[j,"mutation"])){
      delta_delta_G <- c(delta_delta_G,data_double[j,"delta_delta_G"]$delta_delta_G)
      DNA_binding <- c(DNA_binding,data_combine[j,"delta_delta_G"])
    }
  }
}
data_paper$delta_delta_G <- delta_delta_G
data_paper$DNA_binding <- DNA_binding

str(data_DNA)
str(data_combine)
str(data_paper)

g <- ggplot(data_paper,aes(x=gal,y=delta_delta_G))
g <- g + geom_point()
g <- g + labs(y="Total Energy ∆∆G (kcal/mol)",x="ß-gal activity (%)")
g

g <- ggplot(data_paper,aes(x=gal,y=DNA_binding))
g <- g + geom_point()
g <- g + labs(y="Interaction ∆∆G (kcal/mol)",x="ß-gal activity (%)")
g

cor(data_paper[,2:4],method = 'spearman')

library(Hmisc)
res_spearman <- rcorr(as.matrix(data_paper[,2:4]),type = "spearman")
res_spearman$P

library(corrplot)
corrplot(res_spearman$P,type="upper",tl.col ="black",tl.srt = 45)

data_paper_2 <- read_excel("../data_from_paper.xlsx",sheet = "G.J. Palm, et al.")

delta_delta_G <- c()
DNA_binding <- c()
for (i in 1:nrow(data_paper_2)){
  for (j in 1:nrow(data_combine)){
    if (grepl(data_paper_2[i,"mutation"]$mutation,data_combine[j,"mutation"])){
      delta_delta_G <- c(delta_delta_G,data_DNA[j,"delta_delta_G"]$delta_delta_G)
      DNA_binding <- c(DNA_binding,data_combine[j,"delta_delta_G"])
    }
  }
}
data_paper_2$delta_delta_G <- delta_delta_G
data_paper_2$DNA_binding <- DNA_binding

g <- ggplot(data_paper_2,aes(x=delta_delta_G,y=`∆∆G`))
g <- g + geom_point()
g <- g + geom_label()
g

g <- ggplot(data_paper_2,aes(x=DNA_binding,y=`∆∆G`))
g <- g + geom_point()
g
