library(ggplot2)
library(phyloseq)
library(extrafont)
library(phyloseq)
library(vegan)
#library(RColorBrewer)
#Load phyloseq to amvis2 converter
devtools::source_gist("8d0ca4206a66be7ff6d76fc4ab8e66c6")
#Set working directory to script directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#load format
loadfonts()
mytheme<- theme(plot.title = element_text(hjust=0.5, family = "Arial", size=12),
                legend.position ="right",
                legend.text = element_text(family = "Arial", size = 10),
                legend.background = element_blank(),
                strip.background = element_blank(),
                strip.placement = "outside",
                strip.text = element_text(family = "Arial", size = 10),
                axis.text = element_text(family = 'Arial', size = 10, color="black"),
                panel.grid = element_blank())
my_palette <- c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")

#Load file source
ps.rds <- readRDS("ps.rarefied.test.dbrda.rds")

#take necessary info
genus <- subset_samples(ps.rds, Group %in% c("FFT","Control"))
# genus <- subset_samples(genus, Group %in% c("CON","FFT","FVT","FRT"))
# genus <- subset_samples(genus, Group %in% c("CON","FFT","iFFT"))

#combine the otu from same genus in phylose?
genus <- tax_glom(genus, "Genus")


genus_tab <- as.data.frame(otu_table(genus))
tax_tab <- phyloseq::tax_table(genus)
tax_tab <- data.table::as.data.table(tax_tab, keep.rownames = TRUE)
index <- match(rownames(genus_tab),tax_tab$rn)
#Because of vOTU duplicate: further reason need to be check
tax_tab$Genus <- paste(tax_tab$Genus,tax_tab$rn, sep = "_")

rownames(genus_tab) <- as.character(tax_tab$Genus)[index]

# write.csv(tax_tab,"glom_genus_Simone.csv")
#keep rows that rowsum>0
genus_tab <- genus_tab[rowSums(genus_tab)>0,]
genus_tab <- t(genus_tab)

#hellinger transformation:  concentrate dataset by kai fang
genus_tab_h <- decostand(genus_tab, method = "hellinger")

mapping <- data.frame(sample_data(genus))
mapping$Group


genus_tab_h <- as.matrix(genus_tab_h)
dis_bray <- vegdist(genus_tab_h, method = "bray")
#subset metadata
env <- subset(mapping, select = c(Group, Litter,Sex))
db_rda <- capscale(genus_tab_h~., env, distance = "bray", add = TRUE)
dbanov <- anova.cca(db_rda, by="term", permutaitons = 999)

dbanov

# db_rda <- capscale(genus_tab_h~Group*Litter, env, distance = "bray", add = TRUE)
# dbanov <- anova.cca(db_rda, by="term", permutaitons = 999)
# dbanov
db_rda <- capscale(genus_tab_h~Group+Litter, env, distance = "bray", add = TRUE)
dbanov <- anova.cca(db_rda, by="term", permutaitons = 999)
dbanov
# In this case, the results show that the interaction between Group and Litter is
#significant (Pr(>F) = 0.011), indicating that the effect of Group on the microbial 
#community composition depends on the level of Litter, and vice versa. 
# db_rda <- capscale(genus_tab_h~Group+Litter, env, distance = "bray", add = TRUE)
# dbanov <- anova.cca(db_rda, by="term", permutaitons = 999)
# dbanov



#extract dbrda using scaling 1
db_rda.scaling1 <- summary(db_rda, scaling=1)
#RsquareAdj() extract R2 ?RsquareAdj() 
r2 <- RsquareAdj(db_rda)
db_rda_noadj <- r2$r.squared #original R2
db_rda_adj <- r2$adj.r.squared       #adjusted R2
db_rda_noadj
db_rda_adj

# recalculate using adjusted r2
#eig is te zheng zhi
db_rda_exp_adj <- db_rda_adj * db_rda$CCA$eig/sum(db_rda$CCA$eig)
db_rda_eig_adj <- db_rda_exp_adj * db_rda$tot.chi

db_rda_exp_adj
sum(db_rda_exp_adj)
db_rda_eig_adj
##above is vector summerize based on 3 axis cap1 cap2 and cap3


# permutation test
#globel tests, on all axis 999 ? anova.cca
db_rda_test <- anova.cca(db_rda, permutations = 999)
db_rda_test
#separate test on each axises
db_rda_test_axis <- anova.cca(db_rda, by = 'axis', permutations = 999)
# adjust P（BH???
db_rda_test_axis$`Pr(>F)` <- p.adjust(db_rda_test_axis$`Pr(>F)`, method = 'BH')
db_rda_test_axis


# co-linear
#The vif.cca() function is used to compute the variance inflation factors (VIF)
#for the environmental variables in a constrained ordination analysis represented 
#by a cca object in the vegan package in R. The VIF measures the degree of
#collinearity among the environmental variables, and values greater than 1 indicate 
#that there is some degree of collinearity.
#When applied to a cca object, vif.cca() computes the VIF for each environmental 
#variable included in the analysis, based on the correlation between each variable 
#and the other variables included in the analysis.
#In the line of code you provided, vif.cca() is applied to the db_rda object,
#which represents the results of a constrained ordination analysis. The resulting 
#output will be a vector of VIF values, one for each environmental variable included in the analysis.
#By examining the VIF values, it is possible to identify whether there is significant 
#collinearity among the environmental variables included in the analysis. If the 
#VIF values are all relatively low (e.g., less than 5), then there is little 
#collinearity and the environmental variables can be considered to be independent. 
#If the VIF values are high (e.g., greater than 10), then there may be significant 
#collinearity, which can affect the interpretation and validity of subsequent
#statistical analyses. In such cases, it may be necessary to either remove some 
#of the highly correlated variables or to use alternative methods that are more robust to collinearity.  
vif.cca(db_rda)
#vPartition of Lingoes adjusted squared Bray distance in dbRDA 
db_rda_vp <- varpart(dis_bray, env['Group'], env['Litter'], scale = FALSE, add = TRUE)
db_rda_vp
##this make dbrda barplot
library(vegan)
tab <- data.frame(var=c("Group", "Litter","Residuals", "Intersect"), Adj_R2=round(c(db_rda_vp$part$indfract[c(1,3,4), 3], 0), digits = 3))
# tab <- data.frame(var=c("Group", "Litter", "Sex","Residuals", "Intersect"), Adj_R2=round(c(dis_bray[["part"]][["indfract"]][c(1,3,4), 3], 0), digits = 3))
head(tab)
tab$var <- factor(tab$var, levels = c("Residuals", "Intersect", "Litter", "FMT" ))
# make a barplot for visualization
p_bar <- ggplot(data = tab, aes(x= var, y= Adj_R2)) +
  geom_bar(stat = "identity", color="lightgray", fill="lightgray", position=position_dodge()) +
  labs(x = "", y = "", title = "Adjusted R-squared") +
  geom_text(aes(label=Adj_R2, y=Adj_R2 + 0.04), vjust=1, color="black",
            position = position_dodge(0.9), size=3.5)+
  coord_flip() + 
  theme_classic() +
  mytheme +
  theme(legend.position = "none")
p_bar
##make the p bar which is dbrda barplot

# Group explained merely
##Here we calculate a number
anova.cca(capscale(dis_bray~Group+Condition(Litter), env, add = TRUE), permutations = 999)
anova.cca(capscale(dis_bray~Litter+Condition(Group), env, add = TRUE), permutations = 999)
# no co-explained

# visualizatiom with ggplot2
#extract r2 and enviroment factors，firtst 2 axies，scaling1
db_rda.scaling1 <- summary(db_rda, scaling = 1)
#The zuo biao of db_rda.scaling , only keep cap1 and cap2
db_rda.site <- data.frame(db_rda.scaling1$sites)[1:2]
db_rda.env <- data.frame(db_rda.scaling1$biplot)[1:2]
db_rda.spe <- data.frame(db_rda.scaling1$species)[1:2]
# only take the co-microbiome mean abundance > 1%


#same with some part above
genus_tab.rel <- transform_sample_counts(genus, function(x) x/sum(x)*100)
genus_tab.rel
genus_tab.rel <- as.data.frame(otu_table(genus_tab.rel))
tax_tab.rel <- as.data.frame(phyloseq::tax_table(genus))
tax_tab.rel$Genus <- paste(tax_tab.rel$Genus,rownames(tax_tab.rel), sep = "_")
tax_tab.rel <- data.table::as.data.table(tax_tab.rel, keep.rownames=TRUE)
index <- match(rownames(genus_tab.rel), tax_tab.rel$rn)
rownames(genus_tab.rel) <- as.character(tax_tab.rel$Genus)[index]
#find the core data from the otu table
genus_core <- genus_tab.rel[rowMeans(genus_tab.rel)>1,]
#find the location of core data on cap1 and cap2 map
db_rda.spe.core <- db_rda.spe[rownames(genus_core),]


#add names and groups
rownames(db_rda.env)
db_rda.env$name <- c("FFT","FRT","FVT","Litter","GroupFFT:Litter","GroupFRT:Litter","GroupFVT:Litter")
rownames(db_rda.site)
db_rda.site$name <- rownames(db_rda.site)
index <- match(db_rda.site$name, rownames(mapping))
db_rda.site$group <- mapping$Group[index]
db_rda.site$Litter <- mapping$Litter[index]

db_rda.spe.core$name <- rownames(db_rda.spe.core)
#calculate the adujsted r squared
#calculate x cap1 and y cap2
exp_adj <- db_rda$CCA$eig/sum(db_rda$CCA$eig)
rda1_exp <- paste('db-RDA1:', round(exp_adj[1]*100, 2), '%')
rda2_exp <- paste('db-RDA2:', round(exp_adj[2]*100, 2), '%')

# get ggplto2 default colors
# get ggplto2 default colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(5)


db_rda.site$Litter <- paste("Litter", db_rda.site$Litter, sep = "_")

p_dbrda <- ggplot(db_rda.site, aes(CAP1, CAP2)) +
  geom_point(aes(color = group, shape = Litter), size = 4) +
  scale_color_manual(values = c(cols[1], cols[2], cols[3], cols[4])) +
  xlim(-1, 1.5) +
  ylim(-1, 1.2) +
  theme_classic() +
  mytheme +
  labs(x = rda1_exp, y = rda2_exp) +
  guides(color = guide_legend(title = "FVT")) +
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  scale_color_manual(values=c("#87c99d","#41b6c4","#A3dbff"),
                     breaks = c( "CON", "FFT",  "iFFT"),
                     labels = c("CON", "FFT",  "iFFT"))+
  # scale_color_manual(values=c("#80cdc1","#018571","#dfc27d","#05B9E2"),
  #                    breaks = c( "CON", "FFT", "FVT", "FRT"),
  #                    labels = c( "CON", "FFT", "FVT", "FRT"))+
  geom_segment(data = db_rda.spe.core, aes(x = 0, y = 0, xend = CAP1 * 0.4, yend = CAP2 * 0.4), 
               arrow = arrow(length = unit(0.2, 'cm')), size = 0.3, color = 'blue') +
  geom_text(data = db_rda.spe.core, aes(CAP1 * 0.43, CAP2 * 0.46, label = name), color = 'black', size = 3)

p_dbrda
 
