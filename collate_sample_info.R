setwd("~/Dropbox/popprocesses_diversity/ILS_latitude_suboscines/")
getwd()

library(rgeos)

# Get Dstat data
data = read.csv("results/dstat_af_subsample250.csv", stringsAsFactors = F)

# Get sample info from Science SI
sample.info <- read.csv("data/Sample_Info_from_Science.csv", stringsAsFactors = F)
sp.map <- sample.info$AOS...Howard.Moore.Name
names(sp.map) <- sample.info$Tip.Name

sample.info[1,]

data[1,]

out <- cbind(data$sp1, as.character(sp.map[data$sp1]), data$sp2, as.character(sp.map[data$sp2]), data$sp3, as.character(sp.map[data$sp3]), data$out, as.character(sp.map[data$out]), data$prefilter_num_snps, data[,7:9])
colnames(out) <- c("Tip1", "Species 1", "Tip2", "Species 2", "Tip3", "Species 3", "Outgroup Tip", "Outgroup Species", "Number SNPs (Pre-Filter)", "D", "Z", "p")

write.table(out, "./figures/Table_S1.csv", sep=",")
