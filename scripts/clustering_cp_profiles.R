# Load the controls (3 replicates)
ctrl1 = read.csv('Documents/zarkolab/publications/bmc_spiroindane_pyrrolidines/Ctrl1.csv')
ctrl1$Name <- mapvalues(ctrl1$Molecule, from=unique(ctrl1$Molecule), names)
rownames(ctrl1) = paste(ctrl1$Name, round(ctrl1$Concentration,4))
nums <- unlist(lapply(ctrl1, is.numeric))  
ctrl1.num = ctrl1[,nums]
ctrl1.num.scaled <- scale(ctrl1.num)
ctrl1.cor = cor(t(ctrl1.num.scaled), use = 'pairwise.complete.obs')

ctrl2 = read.csv('Documents/zarkolab/publications/bmc_spiroindane_pyrrolidines/Ctrl2.csv')
rownames(ctrl2) = paste(ctrl2$Molecule, round(ctrl2$Concentration, 5))
nums <- unlist(lapply(ctrl2, is.numeric))  
ctrl2.num = ctrl2[,nums]
ctrl2.num.scaled <- scale(ctrl2.num)
ctrl2.cor = cor(t(ctrl2.num.scaled), use = 'pairwise.complete.obs')

ctrl3 = read.csv('Documents/zarkolab/publications/bmc_spiroindane_pyrrolidines/Ctrl3.csv')
rownames(ctrl3) = paste(ctrl3$Molecule, round(ctrl3$Concentration, 5))
nums <- unlist(lapply(ctrl3, is.numeric))  
ctrl3.num = ctrl3[,nums]
ctrl3.num.scaled <- scale(ctrl3.num)
ctrl3.cor = cor(t(ctrl3.num.scaled), use = 'pairwise.complete.obs')

# Load the test compounds (3 replicates)
rep1 = read.csv('Documents/zarkolab/publications/bmc_spiroindane_pyrrolidines/Rep1.csv')
rownames(rep1) = paste(rep1$Molecule, round(rep1$Concentration, 5))
nums <- unlist(lapply(rep1, is.numeric))  
rep1.num = rep1[,nums]
rep1.num.scaled <- scale(rep1.num)
rep1.cor = cor(t(rep1.num), use = 'pairwise.complete.obs')

rep2 = read.csv('Documents/zarkolab/publications/bmc_spiroindane_pyrrolidines/Rep2.csv')
rownames(rep2) = paste(rep2$Molecule, round(rep2$Concentration, 5))
nums <- unlist(lapply(rep2, is.numeric))  
rep2.num = rep2[,nums]
rep2.num.scaled <- scale(rep2.num)
rep2.cor = cor(t(rep2.num.scaled), use = 'pairwise.complete.obs')

rep3 = read.csv('Documents/zarkolab/publications/bmc_spiroindane_pyrrolidines/Rep3.csv')
rownames(rep3) = paste(rep3$Molecule, round(rep3$Concentration, 5))
nums <- unlist(lapply(rep3, is.numeric))  
rep3.num = rep3[,nums]
rep3.num.scaled <- scale(rep3.num)
rep3.cor = cor(t(rep3.num.scaled), use = 'pairwise.complete.obs')

# Aggregate all concentrations of individual compounds
rownames(rep3.agg) = rep3.agg$Group.1
nums <- unlist(lapply(rep3.agg, is.numeric))
rep3.agg.num = rep3.agg[,nums]
rep3.agg.num.scaled <- scale(rep3.agg.num)
rep3.agg.cor = cor(t(rep3.agg.num.scaled), use ='pairwise.complete.obs')

rep1.agg = aggregate(rep1, list(rep1$Molecule),mean)
rownames(rep1.agg) = rep1.agg$Group.1
nums <- unlist(lapply(rep1.agg, is.numeric))
rep1.agg.num = rep1.agg[,nums]
rep1.agg.num.scaled <- scale(rep1.agg.num)
rep1.agg.cor = cor(t(rep1.agg.num.scaled), use ='pairwise.complete.obs')

rep2.agg = aggregate(rep2, list(rep2$Molecule),mean)
rownames(rep2.agg) = rep2.agg$Group.1
nums <- unlist(lapply(rep2.agg, is.numeric))
rep2.agg.num = rep2.agg[,nums]
rep2.agg.num.scaled <- scale(rep2.agg.num)
rep2.agg.cor = cor(t(rep2.agg.num.scaled), use ='pairwise.complete.obs')

ctrl1.agg = aggregate(ctrl1, list(ctrl1$Name),mean)
rownames(ctrl1.agg) = ctrl1.agg$Group.1
nums <- unlist(lapply(ctrl1.agg, is.numeric))
ctrl1.agg.num = ctrl1.agg[,nums]
ctrl1.agg.num.scaled <- scale(ctrl1.agg.num)
ctrl1.agg.cor = cor(t(ctrl1.agg.num.scaled), use ='pairwise.complete.obs')

ctrl2.agg = aggregate(ctrl2, list(ctrl2$Molecule),mean)
rownames(ctrl2.agg) = ctrl2.agg$Group.1
nums <- unlist(lapply(ctrl2.agg, is.numeric))
ctrl2.agg.num = ctrl2.agg[,nums]
ctrl2.agg.num.scaled <- scale(ctrl2.agg.num)
ctrl2.agg.cor = cor(t(ctrl2.agg.num.scaled), use ='pairwise.complete.obs')

ctrl3.agg = aggregate(ctrl3, list(ctrl3$Molecule),mean)
rownames(ctrl3.agg) = ctrl3.agg$Group.1
nums <- unlist(lapply(ctrl3.agg, is.numeric))
ctrl3.agg.num = ctrl1.agg[,nums]
ctrl3.agg.num.scaled <- scale(ctrl3.agg.num)
ctrl3.agg.cor = cor(t(ctrl3.agg.num.scaled), use ='pairwise.complete.obs')

# Extract fingerprints of active compounds from replicate 1 of test compounds
actives = rep1[rep1$Molecule == 'XK00001384-001' & rep1$Concentration == 50, ]
rbind(actives, rep1[rep1$Molecule == 'XK00001401-001' & rep1$Concentration == 50, ])
rbind(actives, rep1[rep1$Molecule == 'XK00001401-001' & rep1$Concentration == 50, ])
actives = rbind(actives, rep1[rep1$Molecule == 'XK00001401-001' & rep1$Concentration == 50, ])
actives = rbind(actives, rep1[rep1$Molecule == 'XK00001401-001' & rep1$Concentration == 16.67, ])
actives = rbind(actives, rep1[rep1$Molecule == 'XK00001401-001' & rep1$Concentration == 5.560, ])
actives = rbind(actives, rep1[rep1$Molecule == 'XK00001387-001' & rep1$Concentration == 50, ])
actives = rbind(actives, rep1[rep1$Molecule == 'XK00001386-001' & rep1$Concentration == 50, ])

# Extract fingerprints of active compounds from replicate 1 of test compounds
actives = rep2[rep2$Molecule == 'XK00001384-001' & rep2$Concentration == 50, ]
rbind(actives, rep2[rep2$Molecule == 'XK00001401-001' & rep2$Concentration == 50, ])
rbind(actives, rep2[rep2$Molecule == 'XK00001401-001' & rep2$Concentration == 50, ])
actives = rbind(actives, rep2[rep2$Molecule == 'XK00001401-001' & rep2$Concentration == 50, ])
actives = rbind(actives, rep2[rep2$Molecule == 'XK00001401-001' & rep2$Concentration == 16.67, ])
actives = rbind(actives, rep2[rep2$Molecule == 'XK00001401-001' & rep2$Concentration == 5.560, ])
actives = rbind(actives, rep2[rep2$Molecule == 'XK00001387-001' & rep2$Concentration == 50, ])
actives = rbind(actives, rep2[rep2$Molecule == 'XK00001386-001' & rep2$Concentration == 50, ])


#append actives to controls and cluster
names_with_actives = c( "etoposide", "berberine HCl", "CA074 Me ester", "SB203580", "rapamycin", 
                        "taxol", "latrunculine B", "cytochalasine B", "rotenone", "NPPB", 
                        "tetrandrine","vincristine", "fenbendazole", 'XK00001384', 'XK00001401', 
                        'XK00001387', 'XK00001386')

map = setNames(unique(ctrl1$Molecule), names)



ctrl1_plus_actives = rbind (ctrl3, actives)
ctrl1_plus_actives$Name = mapvalues(ctrl1_plus_actives$Molecule, from=unique(ctrl1_plus_actives$Molecule), to=names_with_actives)
ctrl1_plus_actives.agg = aggregate(ctrl1_plus_actives, list(ctrl1_plus_actives$Name),mean)
rownames(ctrl1_plus_actives.agg) = ctrl1_plus_actives.agg$Group.1
nums <- unlist(lapply(ctrl1_plus_actives.agg, is.numeric))  
ctrl1_plus_actives.agg.num = ctrl1_plus_actives.agg[,nums]
ctrl1_plus_actives.agg.num.scaled <- scale(ctrl1_plus_actives.agg.num)
ctrl1_plus_actives.agg.cor = cor(t(ctrl1_plus_actives.agg.num.scaled), use = 'pairwise.complete.obs')
plot(hclust(dist(ctrl1_plus_actives.agg.num)), main = 'Actives from rep2')


