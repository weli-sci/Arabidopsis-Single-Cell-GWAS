# Goal: create a plot where x axis is gene location and y axis in hit/SNP location.

# What we need:
# list of hits
# gene locations
# a way to plot continuous variable (base pair) nested within an ordered categorical vairalbe (chr)
# on the same axis. (potential code here: https://www.r-bloggers.com/2019/03/using-r-plotting-the-genome-on-a-line/)
### get list of hits, can use "celltype_specificity.txt"

library(dplyr)
library(tidyr)
library(data.table)

# Get data#####
# setwd("~/Documents/Research_projects/Single_cell/roots/")
# dat <- read.table("celltype_specificity.txt", header = T)
dat <- read.table("celltype_specificity_CellType.txt", header = T) #this file includes info in cell types
names(dat)
str(dat)

dat <- separate(dat, hit, c("Chr", "bp", "gene"), sep = "_")
dat$Freq <- NULL
dat$p <- NULL
dat$gene <- as.character(dat$gene)
dat$bp <- as.integer(dat$bp)
dat$cell_type <- as.character(dat$cell_type)

gID <- read.table("gene_index_metadata.txt", header = T)
colnames(gID) <- c("gene_name", "gene")
table(unique(gID$gene))
levels(gID$gene_name)
gPOS <- read.table("gene_start.txt", header = F)
colnames(gPOS) <- c("gChr", "gbp", "gene_name")
head(gPOS)
head(gID)

gdat <- left_join(gID, gPOS, "gene_name")
head(gdat)
str(gdat)
gdat$gChr <- as.character(gdat$gChr)
gdat$gene <- as.character(gdat$gene)

names(gdat)
table(unique(gdat$gene))
str(dat)
str(gdat)
d <- left_join(dat, gdat, "gene")
head(d)

library (ggplot2)
plot(d$bp,d$gbp)

names(d)
summary(d)
str(d)
d$bp<-as.integer(d$bp)

table(d$Chr,d$gChr)
table(d2$Chr,d2$gChr)
d2<-d[!(d$gChr=="C" | d$gChr=="M" | is.na(d$gChr)),]
d<-d2

table(d$cell_type)
sum(table(d$cell_type))


# Renaming cell types#####
d2$cell_type<-as.factor(d2$cell_type)
d<-d2
d$cell_typeCode<-d$cell_type
# table(d$cell_type)
# table(d$cell_typeCode)
table(d$cell_type, d$cell_typeCode)
levels(d$cell_type)[levels(d$cell_typeCode)=="1"] <- "Cortex"
levels(d$cell_type)[levels(d$cell_typeCode)=="2"] <- "Hair cells"
levels(d$cell_type)[levels(d$cell_typeCode)=="3"] <- "Columella"
levels(d$cell_type)[levels(d$cell_typeCode)=="4"] <- "Non-hair cells"
levels(d$cell_type)[levels(d$cell_typeCode)=="5"] <- "Endodermis"
levels(d$cell_type)[levels(d$cell_typeCode)=="6"] <- "Phloem pole pericycle"
levels(d$cell_type)[levels(d$cell_typeCode)=="7"] <- "Xylem"
levels(d$cell_type)[levels(d$cell_typeCode)=="8"] <- "Xylem pole pericycle"
levels(d$cell_type)[levels(d$cell_typeCode)=="9"] <- "Root cap tip"
levels(d$cell_type)[levels(d$cell_typeCode)=="10"] <- "Root cap early"
levels(d$cell_type)[levels(d$cell_typeCode)=="11"] <- "Root cap lateral"
levels(d$cell_type)[levels(d$cell_typeCode)=="12"] <- "Epidermis meristem"
levels(d$cell_type)[levels(d$cell_typeCode)=="13"] <- "Procambium"
levels(d$cell_type)[levels(d$cell_typeCode)=="14"] <- "Phloem and companion cells"
levels(d$cell_type)[levels(d$cell_typeCode)=="15"] <- "Protoxylem"




# Writing the locations flat along the axes####
Chr_lengths <- summarise(group_by(d, Chr), length = max(bp))
gChr_lengths <- summarise(group_by(d, gChr), length = max(gbp))

flatten_coordinatesChr <- function(Chr, coord, Chr_lengths) {
  coord_flat <- coord
  offset <- 0
  for (contig_ix in 1:nrow(Chr_lengths)) {
    on_contig <- Chr == Chr_lengths$Chr[contig_ix]
    coord_flat[on_contig] <- coord[on_contig] + offset
    offset <- offset + Chr_lengths$length[contig_ix]
  }
  coord_flat
}

flatten_coordinatesgChr <- function(gChr, coord, gChr_lengths) {
  coord_flat <- coord
  offset <- 0
  for (contig_ix in 1:nrow(gChr_lengths)) {
    on_contig <- gChr == gChr_lengths$gChr[contig_ix]
    coord_flat[on_contig] <- coord[on_contig] + offset
    offset <- offset + gChr_lengths$length[contig_ix]
  }
  coord_flat
}

d$bp_flat <- flatten_coordinatesChr(d$Chr, d$bp, Chr_lengths)
d$gbp_flat <- flatten_coordinatesgChr(d$gChr,d$gbp, gChr_lengths)
Chr_lengths$length_flat <- flatten_coordinatesChr(Chr_lengths$Chr, Chr_lengths$length, Chr_lengths)
gChr_lengths$length_flat <- flatten_coordinatesgChr(gChr_lengths$gChr, gChr_lengths$length, gChr_lengths)

axis_coordChr <- c(0, Chr_lengths$length_flat[-nrow(Chr_lengths)])
axis_coordgChr <- c(0, gChr_lengths$length_flat[-nrow(gChr_lengths)])

# Ploting cell types separately####
library(ggplot2)
library(viridis)
library(ggpointdensity)
library(gridExtra)

levels(d$cell_type)
table(d$cell_type)

summary(d)

write.csv(d, "d.csv", row.names = F)
d<-read.csv("d.csv", header = T)

plotsvec<-vector("list",length=length(unique(d$cell_type)))
for (i in unique(d$cell_type)) {   #could also use levels(d$cell_type) as it is a factor
  p<-ggplot(d[d$cell_type==i,], aes(x = bp_flat, y = gbp_flat))+
    geom_pointdensity(na.rm = T) +
    scale_color_viridis(direction = 1, alpha = 0.35,option = "H")+
    geom_vline(aes(xintercept = length_flat), data = gChr_lengths,na.rm = T) +
    scale_y_continuous(breaks = axis_coordgChr,
                       labels = gChr_lengths$gChr,
                       limits = c(0, max(gChr_lengths$length_flat))) +
    scale_x_continuous(breaks = axis_coordChr,
                       labels = Chr_lengths$Chr,
                       limits = c(0, max(Chr_lengths$length_flat))) +
    ylab("Expressed Gene Chromossomic Start Location") + xlab ("eQTL Chromossomic Location") + 
    ggtitle (i)+ theme_bw()
  
  plotsvec[[i]]<-p
}
plotsvec[1:15]<-NULL #because it was adding 15 extra NULL lists int he beginning for some reason
summary(plotsvec)
plotsvec

# arrange the single plots into a a file with a page per graph####
library(gridExtra)
ggsave(
  filename = "Cell-type plot in each page.pdf", 
  plot = marrangeGrob(plotsvec, nrow=1, ncol=1)) #if I change the nrow/ncol values, it starts putting graphs int he same page


# Generate all single plots into a single page/panel

# Show the actual data, the dots, with an indication of density - has a better idea of cis/trans within each cell type 
 p<-ggplot(d, aes(x = bp_flat, y = gbp_flat))+
    facet_wrap(cell_type ~ .,scales = "free", ncol = 3)+ #this separate plots by cell types, facet_wrap is used for a variable with many levels
    geom_pointdensity() +
    scale_color_viridis(direction = 1, alpha = 0.35,option = "H")+
    geom_vline(aes(xintercept = length_flat), data = gChr_lengths) +
    scale_y_continuous(breaks = axis_coordgChr,
                       labels = gChr_lengths$gChr,
                       limits = c(0, max(gChr_lengths$length_flat))) +
    scale_x_continuous(breaks = axis_coordChr,
                       labels = Chr_lengths$Chr,
                       limits = c(0, max(Chr_lengths$length_flat))) +
    ylab("Expressed Gene Chromossomic Start Location") + xlab ("eQTL Chromossomic Location") + 
    ggtitle ("Cis-trans regulatory locations by cell-type")+ theme_bw()
  
 # Density of associations across cell types (visualysing the count (density) of dots inside hexagonal areas)
ggplot(d, aes(x = bp_flat, y = gbp_flat))+
    facet_wrap(cell_type ~ .,scales = "free", ncol = 3)+ #this separate plots by cell types, facet_wrap is used for a variable with many levels
  geom_hex (bins=3, aes(alpha= ..count..))+
  # scale_color_viridis(direction = 1, alpha = 0.35,option = "H")+
    geom_vline(aes(xintercept = length_flat), data = gChr_lengths) +
    scale_y_continuous(breaks = axis_coordgChr,
                       labels = gChr_lengths$gChr,
                       limits = c(0, max(gChr_lengths$length_flat))) +
    scale_x_continuous(breaks = axis_coordChr,
                       labels = Chr_lengths$Chr,
                       limits = c(0, max(Chr_lengths$length_flat))) +
    ylab("Expressed Gene Chromossomic Start Location") + xlab ("eQTL Chromossomic Location") + 
    ggtitle ("Cis-trans regulatory locations by cell-type")+ theme_bw()
  
# save.image(p,file = "Cis trans plots per cell type.jpeg")


# All cell types together plot#####
plot_AllCellT <- ggplot(d, aes(y = gbp_flat,
                               x = bp_flat)) +
  geom_pointdensity() +
  scale_color_viridis(direction = 1, alpha = 0.25,option = "H")+
  geom_vline(aes(xintercept = length_flat),
             data = gChr_lengths) +
  scale_y_continuous(breaks = axis_coordgChr,
                     labels = gChr_lengths$gChr,
                     limits = c(0, max(gChr_lengths$length_flat))) +
  scale_x_continuous(breaks = axis_coordChr,
                     labels = Chr_lengths$Chr,
                     limits = c(0, max(Chr_lengths$length_flat))) +
  ylab("Expressed Gene Chromossomic Start Location") + xlab ("eQTL Chromossomic Location") + theme_bw()

names(d)




