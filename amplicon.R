install.packages("ggprism")
BiocManager::install("microbiome")
library(ggprism)
library(devtools)
library(phyloseq)
library(ggplot2)# Load the devtools package
# install_github("microbiome/microbiome") # Install the package

ps.rarefied <- readRDS("Desktop/16s results/ps.rarefied.test.dbrda.rds")

ps.rarefied.nop <- phyloseq::subset_samples(ps.rarefied, Group %in% c("FFT","Control"))
otus <- as.data.frame(otu_table(ps.rarefied.nop))
metadata <- sample_data(ps.rarefied.nop)
metadata <- subset(metadata, select=c("Litter", "Group", "Sex"))
sample_data(ps.rarefied.nop) <-metadata

cap_ordistep <- function(otus, metadata){
  # remove single level metadata
  uniques <- apply(metadata, 2, function(x) length(unique(x)))
  metadata <- metadata[,uniques!=1]
  # ordistep
  mod0 <- capscale(t(otus) ~ 1, metadata, distance = "bray")
  mod1 <- capscale(t(otus) ~ ., metadata, distance = "bray")
  mod.bray <- ordistep(mod0, scope = formula(mod1), permutations = how(nperm = 999))
  return(mod.bray$anova)
}

plot_cap <- function(mod.c, color_c, shape_c){
  pcs <- mod.c$CA$u
  caps <- mod.c$CCA$u
  pc1.explained <- round(mod.c$CA$eig[1]/mod.c$tot.chi * 100, 1)
  cap1.explained <- round(mod.c$CCA$eig[1]/mod.c$tot.chi * 100, 1)
  
  p.matrix <- cbind(pcs, caps, metadata)
  var.color <- eval(substitute(color_c), p.matrix)
  var.shape <- eval(substitute(shape_c), p.matrix)
#  p1 <- ggplot(p.matrix,
#               aes(x= CAP1, y=MDS1, color=var.color, shape=var.shape)) +
#    geom_point(size=3) +
#    labs(x= paste0("CAP1 (",cap1.explained, " %)"), y= paste0("MDS1 (", pc1.explained, " %)")) +
#    ggtitle("Effect of FFT corrected for litter effect") +
#    theme_prism() +
#    scale_shape_prism()
  
#  return(p1)
#}

  theme_prism <- function() {
    theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  }
  
  p1 <- ggplot(p.matrix,
               aes(x= CAP1, y=MDS1, color=var.color, shape=var.shape)) +
    geom_point(size=3) +
    ggtitle("Effect of FFT corrected for litter effect") +
    labs(x= paste0("CAP1 (",cap1.explained, " %)"), y= paste0("MDS1 (", pc1.explained, " %)")) +
    theme_prism() +
    scale_shape_prism()
  
  return(p1)
}

#p1.0 <- print(p1 + ggtitle("Effect of FFT corrected for litter effect"))
  
#-----------------all samples
ps.rarefied.a <- prune_taxa(rowSums(otu_table(ps.rarefied.nop))!=0, ps.rarefied.nop)
otus <- as.data.frame(otu_table(ps.rarefied.a))
#hellinger transform
#otus.hell <- decostand (otus, 'hell')
metadata <- data.frame(sample_data(ps.rarefied.a))
mod.a <- cap_ordistep(otus, metadata) #Substrate Donor
library(vegan)
#------------------model
mod.f <- cap_ordistep(otus, metadata)
# f1 Donor as conditional factor
mod.f1 <- capscale(t(otus) ~ Litter + Condition(Group), data = metadata, distance = "bray")
#anova(mod.bray.f1, by = "term")
p.cap.f1 <- plot_cap(mod.f1, Litter, Group)

# f2 Substrate as conditional factor
mod.f2 <- capscale(t(otus) ~ Group + Condition(Litter), data = metadata, distance = "bray")
p.cap.f2 <- plot_cap(mod.f2, Group, Litter)


# save
dir.create("results/dbRDA", showWarnings = FALSE)
p.ggplots <- list(p.cap.f1, p.cap.f2)
p.names <- c("bray-Substrate-DonorCON", "bray-Donor-SubstrateCON")

ggsave.adj <- function(x,y){
  ggsave(paste0("db_rda.R", x,".pdf"), y, device = cairo_pdf, height = 4, width = 6)
}

mapply(ggsave.adj, p.names, p.ggplots)

#install.packages("Cairo")
#library(Cairo)

#dir.create("Desktop/16s results/db_rda.R", showWarnings = FALSE)
#p.ggplots <- list(p.cap.f1, p.cap.f2)
#p.names <- c("bray-Substrate-DonorCON", "bray-Donor-SubstrateCON")

#ggsave.adj <- function(x, y) {
#  ggsave(paste0("Desktop/16s results/db_rda.R/", x, ".pdf"), y, device = "pdf", height = 4, width = 6)
#}

#mapply(ggsave.adj, p.names, p.ggplots)

