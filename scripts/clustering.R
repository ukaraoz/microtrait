library(dplyr)
library(stringr)
library(tidyr)
library(vegan)
library(pheatmap)
library(viridis)
library(foreach)

base = "~/Work/microtrait/code/microtrait-scratch/1"
base.out = "/Users/ukaraoz/Work/microtrait/code/AGU-microtrait-vignettes"

####### read metadata
load("~/Work/microtrait/code/microtrait/data/all_public_wmetadata.rda")
all_public_wmetadata = as_tibble(all_public_wmetadata)
cols.select = c("taxon_oid", "domain", "phylum", "ir_class","ir_order","family","genus", "species", "Strain", 
                "Taxon DOI", "Ecosystem", "Ecosystem Category","Ecosystem Type", "Ecosystem Subtype", "Specific Ecosystem",
                "Gene Count", "Estimated Size", "LATITUDE", "LONGITUDE", "GEOGRAPHIC_LOCATION", 
                "Cell Diameter",  "Cell Shape", "Color","Gram Stain", "Motility", "Oxygen Requirement", 
                "Ph","Salinity","Pressure","Sporulation","Carbon Source", "Symbiont Name","Symbiont Taxon Id",
                "Symbiotic Physical Interaction ", "Symbiotic Relationship","Temperature Optimum","Temperature Range")
#empty_counts = all_public_wmetadata %>%
#  summarise_all(funs(sum(is.na(.)|.=="")))/nrow(all_public_wmetadata)
all_public_wmetadata = all_public_wmetadata %>% select(cols.select)
#######
####### read growthpred
growthpred = readRDS("/Users/ukaraoz/Work/microtrait/growthpred/all_public.tbl.taxonoid.bacarc.growthpred.rds")
colnames(growthpred)[1] = "taxon_oid"

# non metabolic
traitmatrix_nonmetabolic = readRDS(file.path(base, "traitmatrix.rds")) %>% as.data.frame
taxon_oid = rownames(traitmatrix_nonmetabolic)
traitmatrix_nonmetabolic = traitmatrix_nonmetabolic %>%
  mutate_all(funs(case_when(. > 0 ~ 1, TRUE ~ 0))) #%>%
  #summarise_all(funs(sum(.)/1594)) %>% as.data.frame
traitmatrix_nonmetabolic = data.frame(taxon_oid, traitmatrix_nonmetabolic, check.names = F, stringsAsFactors = F)

# metabolic
load(file.path(base, "traitmatrix_binarytraits.RData")) # 1594 x 90
traitmatrix_binarytraits = traitmatrix_binarytraits[, c(1, 9:46, 48:90)]
traitmatrix_metabolic = traitmatrix_binarytraits  %>% tbl_df

alltraits = traitmatrix_nonmetabolic %>%
  inner_join(traitmatrix_metabolic, by = c("taxon_oid" = "taxon_oid")) %>%
  mutate_at(vars(-taxon_oid), funs(as.integer))

#rownames = alltraits[, "taxon_oid"]
#alltraits = data.frame(alltraits[, 2:ncol(alltraits)], check.names = F)
#rownames(alltraits) = rownames

#saveRDS(alltraits, file = "/Users/ukaraoz/Work/microtrait/code/AGU-microtrait-vignettes/alltraits.rds")
alltraits = readRDS("/Users/ukaraoz/Work/microtrait/code/AGU-microtrait-vignettes/alltraits.rds")
#braydist = as.dist(vegdist(traitmatrixnormperc, method = "bray"))
#alltraits.dist = as.dist(dist(alltraits,method = "binary"))
#alltraits = alltraits[which(apply(alltraits, 1, sum) != 0),]
alltraits.dist = vegdist(alltraits,method = "jaccard")

all = alltraits %>% 
  #select(-domain, -phylum, -ir_class, -ir_order, -family, -genus, -species, -sdgentime, -nHEG, -nNonHEG) %>%
  inner_join(all_public_wmetadata) %>%
  inner_join(growthpred, by = c("taxon_oid" = "taxon_oid")) %>%
  mutate(`Ecosystem Subtype` = replace(`Ecosystem Subtype`, `Ecosystem Subtype`== "", "Unclassified")) %>%
  mutate(`Specific Ecosystem` = replace(`Specific Ecosystem`, `Specific Ecosystem`== "", "Unclassified")) %>%
  mutate(`Specific Ecosystem` = replace(`Specific Ecosystem`, `Specific Ecosystem`== "Forest Soil", "Forest soil")) %>%
  mutate(Ecosystem = paste0(`Ecosystem Subtype`, ": ", `Specific Ecosystem`)) %>%
  select("taxon_oid", cols.select, everything())

allsampled = all %>%
  group_by(species) %>% 
  top_n(n=1,wt=taxon_oid) %>% as.data.frame

# Generate annotations for rows and columns
library(pals)
#colors = glasbey(n=32)
colors = polychrome(n=36)
allphyla = all[,"phylum"]%>%unique
colorsphyla = colors[1:length(allphyla)]
names(colorsphyla) = allphyla

allclass = all[,"ir_class"]%>%unique
colorsclass = colors[1:length(allclass)]
names(colorsclass) = allclass

allecosystems = all[,"Ecosystem"]%>%unique
colorsecosystems = colors[1:length(allecosystems)]
names(colorsecosystems) = allecosystems

annotation_row = all%>%select(domain, phylum, ir_class, Ecosystem) %>% as.data.frame
rownames(annotation_row) = all %>% select("taxon_oid") %>% pull %>% as.character
annotation_colors = list(
  domain = c(Bacteria = "red", Archaea = "blue"),
  phylum = colorsphyla,
  class = colorsclass,
  ecosystem = colorsecosystems
)

colnames2remove = c("transport:osmolyte", "Nitrogen: ammonia assimilation to glutamate-GS_GOGAT", "Nutrient Utilization: fructose utilization", "Nitrogen: ammonia assimilation to glutamine-GS",
"Methanotrophy: PHB cycle", "Nitrogen: ammonia assimilation to glutamate-GDH", "transport:nitrogen compound", "uptake:organophosphorus compound",
"transport:siderophore", "PRPP biosynthesis", "drug/toxin resistance",  "transport:ion",
"uptake:oligosaccharide", "transport:lipid (cell membrane)", "uptake:monosaccharide", "uptake:vitamin",
"transport:carboxylic acid", "transport:metal", "Nitrogen: glutamine to aspartate", "Methanotrophy: serine pathway")
rownames = allsampled[, "taxon_oid"]
colnames = setdiff(colnames(allsampled)[38:174], colnames2remove)
toplot = allsampled %>% select(colnames)
rownames(toplot) = rownames
pheatmap(toplot, 
         color = c("white", "red"), show_rownames = F, show_colnames = T,
         #cellwidth = 1, cellheight = 1, labels_row = NA, border_color = NA, treeheight_row = 50,
         width = 20*(16/9), height = 20, fontsize_row = 6, fontsize_col = 20, angle_col = 90, border = "black",
         cluster_cols = F, clustering_distance_rows = "binary",
         #annotation_row = annotation_row, annotation_colors = annotation_colors,
         silent = T, filename = "/Users/ukaraoz/Work/microtrait/code/AGU-microtrait-vignettes/alltraits.heatmap.pdf")

results.binarydist.adonis = readRDS("/Users/ukaraoz/Work/microtrait/code/AGU-microtrait-vignettes/results.binarydist.adonis.rds")
results.binarydist.adonis = data.frame(`number of clusters` = results.binarydist.adonis[,1],
                                       `percent variance` = round(results.binarydist.adonis[,2]/results.binarydist.adonis[,3]*100,2),
                                       check.names = F)
library(ggplot2)
ggplot(results.binarydist.adonis,aes(x=`number of clusters`, y=`percent variance`)) +
  geom_line() + geom_point(size=0.8) + 
  #xlim(1, 500) + 
  geom_hline(yintercept = 80, color = "red", size = 0.2) + 
  geom_vline(xintercept = 276, color = "red", size = 0.2) +
  geom_hline(yintercept = 70, color = "red", size = 0.2) + 
  geom_vline(xintercept = 153, color = "red", size = 0.2) +
  geom_hline(yintercept = 60, color = "red", size = 0.2) + 
  geom_vline(xintercept = 84, color = "red", size = 0.2) +
  scale_y_continuous(breaks = seq(0, 100, by = 10))
  

# define dendrogram object to play with:
hc <- hclust(alltraits.dist, "average")
#dend <- as.dendrogram(hc)
#results = matrix(nrow = 0, ncol =3)
cl <- parallel::makeForkCluster(30)
doParallel::registerDoParallel(cl)
results = foreach(i = 2:attr(alltraits.dist, "Size"), .combine = 'rbind') %dopar% {
  cat(i, "\n")
  o.relgr = cutree(hc, k=i)
  o.adonis = adonis(alltraits.dist ~ as.factor(o.relgr))
  c(i, o.adonis$aov.tab$SumsOfSqs[c(1,3)])
}
saveRDS(results, file = "results.jaccarddist.adonis.rds")
parallel::stopCluster(cl)

for(i in 2:attr(alltraits.dist, "Size")){
  cat(k, "\n")
  o.relgr = cutree(hc, k=i)
  o.adonis = adonis(alltraits.dist ~ as.factor(o.relgr))
  #o.adonis$aov.tab$SumsOfSqs[c(1,3)]
  results = rbind(results,
                  c(k, o.adonis$aov.tab$SumsOfSqs[c(1,3)]))
}

########
#### tradeoffs between degradation enzymes and transporters
########
all = readRDS("/Users/ukaraoz/Work/microtrait/code/microtrait-scratch/1/alltraitswmetadata.terrestrial.rds")
traitlist1 = grep("degradation:|Nutrient Utilization:", colnames(all), value = T, perl = T)
traitlist2 = grep("transport:|uptake:", colnames(all), value = T, perl = T)

all = all %>%
  #select(c(traitlist1, traitlist2)) %>%
  mutate_at(vars(`degradation:cellulose`:`Chlorate: chlorate reduction`), 
            funs(case_when(.>=1 ~ 1, TRUE ~ 0))) %>%
  group_by(species) %>% 
  top_n(n=1,wt=taxon_oid) %>% as.data.frame
temp = all %>% 
  select(c(traitlist1, traitlist2)) %>% 
  #select(-adhesion) %>%                # remove unnecessary columns
  cor(method = "pearson") %>%                      # get all correlations (even ones you don't care about)
  data.frame(check.names = F) %>%               # save result as a dataframe
  mutate(v1 = row.names(.)) %>%  # add row names as a column
  gather(v2,cor, -v1) %>%        # reshape data
  filter(v1 != v2) %>%
  arrange(desc(-cor))
#write.table(temp, file = "~/Work/microtrait/code/AGU-microtrait-vignettes/corr.binary.xls",
#            row.names = F, col.names = T, sep = "\t", quote = F)

# Generate annotations for rows and columns
colors = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00","#cab2d6",
           "#fbb4ae", "#b3cde3", "#ccebc5", "#decbe4", "#fed9a6", "#ffffcc", "#e5d8bd", "#fddaec", "#f2f2f2", 
           "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999", 
           "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9",
           "#b2182b", "#d6604d", "#f4a582")
allphyla = all[,"phylum"]%>%unique
colorsphyla = colors[1:length(allphyla)]
names(colorsphyla) = allphyla

allclass = all[,"ir_class"]%>%unique
colorsclass = colors[1:length(allclass)]
names(colorsclass) = allclass

annotation_row = all%>%select(domain, phylum, ir_class) %>% as.data.frame
rownames(annotation_row) = all %>% select("taxon_oid") %>% pull %>% as.character
annotation_colors = list(
  domain = c(Bacteria = "red", Archaea = "blue"),
  phylum = colorsphyla,
  class = colorsclass
)

rownames = all[, "taxon_oid"]
toplot = all%>%select(c(traitlist1, traitlist2))
rownames(toplot) = rownames
pheatmap(toplot, 
         color = c("white", "red"), show_rownames = F,
         #cellwidth = 1, cellheight = 1, labels_row = NA, border_color = NA, treeheight_row = 50,
         width = 18, height = 12, fontsize_row = 6, fontsize_col = 6, angle_col = 45,
         cluster_cols = T, clustering_distance_rows = "binary",
         annotation_row = annotation_row, annotation_colors = annotation_colors,
         silent = T, filename = "/Users/ukaraoz/Work/microtrait/code/AGU-microtrait-vignettes/binary.nonmetabolic.heatmap.pdf")
######################
# Hi Ulas (cc Gianna and Jinyun) - some quick questions and I'll follow up with you later.
# 
# In genomes, what covaries with predicted Umax?
# 
# I'm interested in:
# Occurrence of EMP glycolysis versus ED glycolysis - we might expect that higher Umax have higher frequency of ED (power versus yield)
# This paper - https://www.pnas.org/content/110/24/10039 - attributes the cost-benefit to differences in 
# protein investment which is only part of the picture as differences in heat loss should be much higher 
# per unit glucose for ED versus EMP.
# 
# Genome size - simply put does a larger genome (or parts associated with energy transduction) 
# have more steps in thermodynamic dissipation than smaller genomes and therefore more yield optimized due to lower heat loss.
all = readRDS("/Users/ukaraoz/Work/microtrait/code/microtrait-scratch/1/alltraitswmetadata.terrestrial.rds")
traitlist = colnames(all)[38:174]
results = matrix(nrow = 0, ncol = 5)
colnames(results) = c("trait", "pvalue", "logodds", "conf_2.5%", "conf_97.5%")
for(i in 1:length(traitlist)){
  cat(traitlist[i], "\n")
  tempbin = all[, traitlist[i]]
  tempbin[tempbin>0] = 1
  if(sum(tempbin) == 0){next}
  mydata = data.frame(trait = tempbin, 
                    mingentime = all[, "mingentime"])
  mylogit <- glm(trait ~ mingentime, data = mydata, family = "binomial")
  pvalue = summary(mylogit)$coefficients[2,"Pr(>|z|)"]
  logodds = summary(mylogit)$coefficients[2,"Estimate"]
  conf_2.5 = confint(mylogit)[2,"2.5 %"]
  conf_97.5 = confint(mylogit)[2,"97.5 %"]
  results = rbind(results, c(traitlist[i], pvalue, logodds, conf_2.5, conf_97.5))
}
results.df = data.frame(trait = results[,"trait"],
                        pvalue = format(as.numeric(results[, "p-value"]), scientific=TRUE),
                        logodds = format(as.numeric(results[, "logodds"]), scientific=TRUE),
                        conf_2.5 = format(as.numeric(results[, "conf_2.5%"]), scientific=TRUE),
                        conf_97.5 = format(as.numeric(results[, "conf_97.5%"]), scientific=TRUE),
                        stringsAsFactors = F)
write.table(results.df, file = "/Users/ukaraoz/Work/microtrait/code/microtrait-scratch/1/mgt-traits-logodds.xls",
            row.names = F, col.names = T, sep = "\t", quote = F)

grep("degradation:|Nutrient Utilization:", colnames(all), value = T, perl = T)

mydata <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")
mydata$rank <- factor(mydata$rank)
mylogit <- glm(admit ~ gre, data = mydata, family = "binomial")


######################

#  select(-taxon_oid) %>%
#  dist(method = "binary") %>% as.dist
pheatmap(alltraits, scale = "none",clustering_distance_rows = alltraits.dist,
         color =  magma(100, direction = -1), cluster_cols=T)
fviz_nbclust(alltraits, FUN = hcut, diss = alltraits.dist, method = "wss", k.max = 1594)

# Dissimilarity matrix


# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "complete" )

library(cati)
library(mice)
library(hypervolume)
library(e1071)
data(finch.ind)

growthpred = readRDS("/Users/ukaraoz/Work/microtrait/growthpred/all_public.tbl.taxonoid.bacarc.growthpred.rds")
