# Relatedness

Calculating relatedness for all individuals. 
(by Vilhelmiina)

```
setwd("C:/Users/Vilhelmiina/Desktop/BO")
library(ggplot2)
BV.names <- read.table(Bombina71.info, header=T)
relatedness2 <- read.table("bombina71.relatedness2", header=T)

#plot histogram of relatedness coefficients
ggplot(data=relatedness2, aes(relatedness2$RELATEDNESS_PHI)) + 
+     geom_histogram(aes(fill=cut(..x.., breaks=c(-1.2, 0.0442, 0.0884, 0.177, 0.354))),binwidth=.01, col="black", alpha = 1, show.legend = FALSE) + 
+     labs(title="Frequency distribution of B. orientalis relatedness2", x="Relatedness2 coefficient", y="Frequency") + xlim(c(-0.5, 0.51)) + ylim(c(0,160)) + scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits=c(-1.2, 0.51))

#filter the dataset to only include coefficients above 0.0442 (threshold for 3rd degree relatedness) and draw another histogram
rel2.filter <- subset(relatedness2, relatedness2$RELATEDNESS_PHI > 0.0442, stringsasFactors=F)
ggplot(data=rel2.filter, aes(rel2.filter$RELATEDNESS_PHI)) + 
+     geom_histogram(aes(fill=cut(..x.., breaks=c(-1.2, 0.0442, 0.0884, 0.177, 0.354))),binwidth=.01, col="black", alpha = 1, show.legend = FALSE) + 
+     labs(title="Frequency distribution of B. orientalis relatedness2", x="Relatedness2 coefficient", y="Frequency") + xlim(c(-0.5, 0.49)) + ylim(c(0,30)) + scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits=c(-1.2, 0.50))

#adjusted axes for the smaller dataset
ggplot(data=rel2.filter, aes(rel2.filter$RELATEDNESS_PHI)) + 
+     geom_histogram(aes(fill=cut(..x.., breaks=c(-1.2, 0.0442, 0.0884, 0.177, 0.354))),binwidth=.01, col="black", alpha = 1, show.legend = FALSE) + 
+     labs(title="Frequency distribution of B. orientalis relatedness2", x="Relatedness2 coefficient", y="Frequency") + xlim(c(-0.7, 0.49)) + ylim(c(0,30)) + scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits=c(-0.6, 0.50))

#order the filtered dataset for editing in excel and write it as a .csv
rel2.order <- rel2.filter[order(rel2.filter$RELATEDNESS_PHI),]
write.csv(rel2.order, "relatedness2_ordered.csv")
