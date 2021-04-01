library(phylopath)

setwd("~/Dropbox/popprocesses_diversity/ILS_latitude_suboscines/")
t = read.tree("data/T400F_AOS_HowardMoore.tre")

n = 250
d = read.csv(paste0("results/dstat_af_subsample", n, ".expanded.csv"), stringsAsFactors = F)
dim(d)
vars = c("Z", "D", "p",
         "divtime2", "divtime3", "internode",
         "length_to_out", "ils2", "ils3",
         "ils4", "rangesize", "overlap",
         "distance", "latitude",
         "all_stability", "hyb_stability")

# check that species are unique
sps = table(d$sp1)
sps[ sps > 1]

d$abslatitude = abs(d$latitude)
d$abs_hyb_stability = abs(d$hyb_stability)
d$absD = abs(d$D)

d3 = d
logvals = c("rangesize", "internode",
        "divtime2", "Z", "divtime3",
        "abslatitude", "abs_hyb_stability",
        "distance")
for (i in 1:length(logvals)) {
  if (logvals[i] == "distance") {
    d3[, logvals[i]] =  log(d3[, logvals[i]] + 1)
  } else {
    d3[, logvals[i]] =  log(d3[, logvals[i]] )
  }
  
}

m <- define_model_set(
  null = c(),
  ils_only = c(ils2 ~ rangesize + internode),
  hyb_only = c(ils2 ~ Z,
               Z ~ divtime3 + abslatitude + abs_hyb_stability + distance),
  all = c(ils2 ~ rangesize + internode,
          ils2 ~ Z,
          Z ~ divtime3 + abslatitude + abs_hyb_stability  + distance)
)

pdf("figures/model_set.pdf", height = 4, width = 10)
plot_model_set(m, algorithm = "kk")
dev.off()

d4 = d3[complete.cases(d3[, vars]), ]
d4 <- d4[!is.infinite(rowSums(d4[, vars])),]
# make downsampled tree
t1 = keep.tip(t, d4$sp1)

d5 = d4[match(t1$tip.label, d4$sp1), ]
rownames(d5) = d5$sp1

p <- phylo_path(m, d5, t1)
p
s <- summary(p)
s
plot(s)
(best_model <- best(p))


# plot model choice
s1 = as.data.frame(s)
models = c("all", "ils_only", "hyb_only", "null")
modelnames = c("ILS & hyb.", "ILS only", "hyb. only", "null")
names(modelnames) = models
s1$modelnames = modelnames[ s1$model ]

s1$weight = s1$w
s1[s1$weight < 0.001, "weight"] = 0.001
wplt = ggplot(s1, aes(x = reorder(modelnames, w), y = weight)) + 
  geom_bar(stat = "identity") +
  coord_flip() + xlab("") +
  ylab("model weight")
save_plot(paste0("figures/model_weight", n, ".pdf"), wplt, 
          base_height = 3, base_width = 4)

# plot best model
pdf(paste0("figures/best_model", n, ".pdf"), height = 4, width = 8)
plot(best_model, text_size=3.5)
dev.off()

# plot correlations
dv = c("Z", "ils2") 
coef = best_model$coef[, dv]
se = best_model$se[ , dv]
sig = abs(coef) > se
sig2 = sig[apply(sig, 1, function(x) {sum(x)}) > 0, ]

a = ggplot(d4, aes(distance, Z)) +
    geom_point(alpha = 0.4, color = "gray40") +
    geom_smooth(method = "lm", fill = NA, col = "navyblue") + 
    ylab("Log Z-score") +
    xlab("Log geographic distance")
b = ggplot(d4, aes(abs_hyb_stability, Z)) +
  geom_point(alpha = 0.4, color = "gray40") +
  geom_smooth(method = "lm", fill = NA, col = "navyblue") + 
  ylab("Log Z-score") +
  xlab("Log range instability")
c = ggplot(d4, aes(abslatitude, Z)) +
  geom_point(alpha = 0.4, color = "gray40") +
  geom_smooth(method = "lm", fill = NA, col = "navyblue") + 
  ylab("Log Z-score") +
  xlab("Log latitude")
d = ggplot(d4, aes(internode, ils2)) +
  geom_point(alpha = 0.4, color = "gray40") +
  geom_smooth(method = "lm", fill = NA, col = "navyblue") + 
  ylab("gene tree discordance") +
  xlab("Log internode distance")
abc = plot_grid(a, b, c, d, ncol = 2, 
                labels = c("B", "C", "D", "E"))
save_plot(paste0("figures/correlations", n, ".pdf"), abc,
          base_height = 5, base_width = 6.5)


# coefficients, if desired
a = coef_plot(best_model, error_bar = "se", 
          order_by="strength", 
          to="Z")
b = coef_plot(best_model, error_bar = "se", 
          order_by="strength", 
          to="ils2")
ab = plot_grid(a, b, ncol = 2)
save_plot(paste0("figures/coef", n, ".pdf"), ab,
          base_height = 4, base_width = 14)

all = p$model_set[["all"]]
dots <- phylopath:::combine_dots(p$dots)
do.call(phylopath:::est_DAG, c(list(all, p$data, 
                        p$tree, 
                        p$model, p$method),
                   dots))

