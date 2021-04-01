library(ape)

# https://gist.github.com/darmitage/2484454
collapse_node <- function(tree, cutoff) {
  tree$node.label = as.numeric(tree$node.label)
  h <- which(tree$node.label < cutoff)
  nodes <- matrix(h + Ntip(tree))
  nodes <- nodes[-1]
  k <- matrix(tree$edge[,2])
  
  j <- matrix(nrow = length(i), ncol = 1)
  
  for (i in 1:length(nodes)){
    j[i] = which(k == nodes[i]) 
  }
  tree$edge.length[j] <- 0
  
  newtree <- di2multi(tree)
}

t = read.tree("~/Dropbox/popprocesses_diversity/ILS_latitude_suboscines/data/phyparts/iqtree_files/iqtree_genetrees.tre")
sh = c(20, 50, 80)
res = vector("list", length(sh))
for (i in 1:length(sh)) {
  res[[i]] = lapply(t, collapse_node, sh[i])
}


for (i in 1:length(sh)) {
 tt = res[[i]]
 class(tt) = "multiPhylo"
 write.tree(tt, paste0("~/Dropbox/popprocesses_diversity/ILS_latitude_suboscines/data/phyparts/iqtree_files/iqtree_genetrees.SH", sh[i], ".trees"))
}
