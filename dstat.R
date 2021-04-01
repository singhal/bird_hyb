library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

calc_dstat <- function(df) {
  d = table(df$sitetype)
  Dstat = (d['ABBA'] - d['BABA']) / (d['ABBA'] + d['BABA'])
  return(Dstat)
}

calc_dstat_af <- function(df) {
  Dstat = (sum(df$abba) - sum(df$baba)) / (sum(df$abba) + sum(df$baba))
  return(Dstat)
}

read_file <- function(x, subsample = FALSE, af = FALSE) {
  x = read.csv(x, stringsAsFactors = F)
  cts = table(x$locus)
  drop = names(cts[cts > 5])
  x1 = x[!x$locus %in% drop, ]
  
  locnum = length(unique(x1$locus))
  snpnum = nrow(x1)
  
  # get species
  sps = strsplit(x$species[1], ":")[[1]]
  
  # calc number of snps
  if (subsample) {
    if (nrow(x1) > subsample) {
      x1 = x1[sample(1:nrow(x1), subsample), ]
    }
  }
  

  # calc bootstrap
  boots = lapply(1:100, function(x) {x1[sample(1:nrow(x1), nrow(x1), replace =  T), ]})
  if (af) {
    dboots = unlist(lapply(boots, calc_dstat_af))
  } else {
    dboots = unlist(lapply(boots, calc_dstat))
  }
  
  if (snpnum > 0) {
    if (af) {
      D = calc_dstat_af(x1)
    } else {
      D = calc_dstat(x1)
    }
    
    # should this be standard error?
    # see both in the literature
    d_sd = sd(dboots)
    d_z = abs(D) / d_sd
    d_p <- 2*pnorm(-abs(d_z))
    vals = c(sps, nrow(x), snpnum, locnum, D, d_z, d_p)
  } else {
    vals = c(sps, nrow(x), snpnum, NA, NA, NA, NA)
  }
  names(vals) = NULL
  return(vals)
}

make_data_frame <- function(res) {
  res2 = as.data.frame(do.call(rbind, res))
  names(res2) = c("sp1", "sp2", "sp3", "out", 
                  "prefilter_num_snps",
                  "num_snps", "num_loci", "D", "Z", "p")
  for (i in 5:10) {
    res2[ , i] =  as.numeric(res2[ , i])  
  }
  return(res2)
}

d = list.files("~/Dropbox/popprocesses_diversity/ILS_latitude_suboscines/data/dstat_shared_PRG/", full.names = T)
nrows = unlist(lapply(d, function(x) { nrow(read.csv(x)) }))

# some files are empty
# everything worked fine
# just no SNPs retained after filter ...
d[ which(nrows == 0) ] 
d1 = d[ which(nrows != 0) ] 

res = lapply(d1, read_file)
res2 = make_data_frame(res)

# these are the triads 
# where many SNPs were filtered out
res2 %>% filter(num_snps / prefilter_num_snps < 0.5)

# so need to subsample down the number of SNPs
# to compare apples to apples
# subsample down to 100 SNPs only
res100 = lapply(d1, read_file, subsample = 100)
res100_2 = make_data_frame(res100)
res100_2 %>% filter(num_snps == 100, p < 0.05)

# subsample down to 200 SNPs only
res200 = lapply(d1, read_file, subsample = 200)
res200_2 = make_data_frame(res200)
res200_2 %>% filter(num_snps == 200, p < 0.05)

# how correlated are these subsamples
res3 = inner_join(res100_2, res200_2, by = c("sp1" = "sp1", 
                                             "sp2" = "sp2",
                                             "sp3" = "sp3"))
# very correlated
res4 = res3 %>% filter(num_snps.y == 200)
cor.test(res4$D.x, res4$D.y)
cor.test(res4$Z.x, res4$Z.y)
cor.test(res4$p.x, res4$p.y)

##############################
# allele frequency approach
##############################
d = list.files("~/Dropbox/popprocesses_diversity/ILS_latitude_suboscines/data/dstat_shared_PRG_af", full.names = T)
# nrows = unlist(lapply(d, function(x) { nrow(read.csv(x)) }))

res = lapply(d, read_file, af = TRUE)
res2 = make_data_frame(res)

res3 = res2 %>% filter(complete.cases(res2$Z), is.finite(res2$Z))
# there is a strong correlation of Z score & pval
# with number of SNPs included
cor.test(res3$num_snps, res3$Z, method = "spearman")
res3$sig = ifelse(res3$p < 0.05, TRUE, FALSE)
aa = aov(res3$sig ~ res3$num_snps)

a = ggplot(res3, aes(num_snps, Z)) + 
  geom_point(col = "gray60", alpha = 0.5) + 
  xlab("num. of sites") +
  geom_smooth(fill = NA, col = "navyblue")
b = ggplot(res3, aes(sig, num_snps)) + 
  geom_boxplot() +
  xlab("significant D-stat") +
  ylab("num. of sites")
ab = plot_grid(a, b, nrow = 1, rel_widths = c(0.6, 0.4),
               labels = c("A", "B"))
save_plot("figures/effect_of_SNP_num.pdf", ab, 
          base_width = 7, base_height = 3)

a = ggplot(res3) + 
  geom_histogram(aes(num_snps)) +
  geom_vline(xintercept = 250, col = "red")
b = ggplot(res3) + 
  geom_histogram(aes(num_loci))
ab = plot_grid(a, b, nrow = 1, labels = c("A", "B"))
save_plot("~/Dropbox/popprocesses_diversity/ILS_latitude_suboscines/figures/snp_loci.png",
          ab, base_width = 7, base_height = 3)

num_snps = c(100, 250, 500)
for (i in num_snps) {
  res = lapply(d, read_file, af = TRUE, subsample = i)
  res3 = make_data_frame(res)
  res3 = res3 %>% filter(num_snps == i)
  outf = paste0("~/Dropbox/popprocesses_diversity/ILS_latitude_suboscines/results/", "dstat_af_subsample", i, ".csv")
  write.csv(res3, outf, row.names = F, quote = F)
}


res23 = inner_join(res2, res3, by = c("sp1" = "sp1",
                                      "sp2" = "sp2",
                                      "sp3" = "sp3"))
cor.test(res23$D.x, res23$D.y)
cor.test(res23$Z.x, res23$Z.y)
res23 %>% filter(p.x > 0.005, p.y < 0.005)

res3$sig = FALSE
res3[which(res3$p < 0.0005), "sig"] = TRUE
dval = ggplot(res3, aes(D)) + geom_histogram(aes(fill = sig))
save_plot("~/Dropbox/popprocesses_diversity/ILS_latitude_suboscines/plots/Dstat.pdf", 
          dval)

# compare to ind PRG prev method
# d1 = list.files("~/Dropbox/popprocesses_diversity/ILS_latitude_suboscines/dstat_individual_PRG/", full.names = T)
# indprg1 = lapply(d1, read_file)
# indprg2 = make_data_frame(indprg1)
# xx = inner_join(res2, indprg2, by = c("sp1" = "sp1", "sp2" = "sp2", "sp3" = "sp3"))
# methods return different results
# looks like shared PRG has marginally more power
# but perhaps more importantly shared PRG 
# does not seem to have as many technical isseus