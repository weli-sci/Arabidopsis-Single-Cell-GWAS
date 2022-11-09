
09/11/2022 #Weli
# Creating a table with the number of hits for each celltype/chromossome, 
# and the percentage of those hits that are in cis regulatory regions

# get files
ciswind<- read.table("gene_cis_window.txt", header = T) #genes start and end
head(ciswind)
tail(ciswind)
names(ciswind)
str(ciswind)

library(dplyr)
str(d) #df with the hits for the respective expressed gene (the expressed gene may apear multiple times if multiple locations were asscoaited)
str(ciswind) #df with the list of expressed genes, a gene appear only once.
ciswind$gChr<-as.integer(ciswind$gChr)
d$gbp_start<-d$gbp
d$gbp<-NULL
d2<-dplyr::inner_join(d, ciswind, by = NULL)
tail(d2)
table(is.na(d2$cis_window_start_bp))
table(is.na(ciswind$cis_window_start_bp))
table(is.na(d2$gbp_start))
table(is.na(ciswind$gbp_start))
table(is.na(d$gbp_start))

# Combine Cell_type and chr

d2$cell_gchr<-paste(d2$cell_typeCode, d2$gChr,sep = "_")

# Count hits for each cell_gchr

dhit<-as.data.frame(table(d2$cell_gchr))
colnames(dhit) <- c("cell_gchr","nhits")

# count how many cis and how many trans per cell_gchr

d2$Cis<-ifelse(d2$bp<d2$cis_window_end_bp & d2$bp>d2$cis_window_start_bp, 1,0)

cis_c<-as.data.frame(table (d2$cell_gchr, d2$Cis))
head(cis_c)
colnames(cis_c)<-c("cell_gchr","Cis","count")
trans_c<-cis_c[cis_c$Cis==0,]
cis_c<-cis_c[cis_c$Cis==1,]
trans_c$Cis<-NULL
cis_c$Cis<-NULL
colnames(trans_c)<-c("cell_gchr","count_trans")
colnames(cis_c)<-c("cell_gchr","count_cis")

# Create table merging files above

djoin<-dplyr::left_join(cis_c,trans_c, by=NULL)
djoin2<-dplyr::left_join(dhit, djoin, by=NULL)
djoin2$Per_cis<-(djoin2$count_cis*100/djoin2$nhits)
names(d)
library(tidyr)
d3<-tidyr::separate(djoin2, cell_gchr, sep = "_", into = c("cell_typeCode", "gChr"))
d3$cell_type<-d3$cell_typeCode
d3$cell_type<-as.factor(d3$cell_type)
d3$cell_typeCode<-as.factor(d3$cell_typeCode)

levels(d3$cell_type)
levels(d3$cell_type)[levels(d3$cell_type)=="1"] <- "Cortex"
levels(d3$cell_type)[levels(d3$cell_typeCode)=="2"] <- "Hair cells"
levels(d3$cell_type)[levels(d3$cell_typeCode)=="3"] <- "Columella"
levels(d3$cell_type)[levels(d3$cell_typeCode)=="4"] <- "Non-hair cells"
levels(d3$cell_type)[levels(d3$cell_typeCode)=="5"] <- "Endodermis"
levels(d3$cell_type)[levels(d3$cell_typeCode)=="6"] <- "Phloem pole pericycle"
levels(d3$cell_type)[levels(d3$cell_typeCode)=="7"] <- "Xylem"
levels(d3$cell_type)[levels(d3$cell_typeCode)=="8"] <- "Xylem pole pericycle"
levels(d3$cell_type)[levels(d3$cell_typeCode)=="9"] <- "Root cap tip"
levels(d3$cell_type)[levels(d3$cell_typeCode)=="10"] <- "Root cap early"
levels(d3$cell_type)[levels(d3$cell_typeCode)=="11"] <- "Root cap lateral"
levels(d3$cell_type)[levels(d3$cell_typeCode)=="12"] <- "Epidermis meristem"
levels(d3$cell_type)[levels(d3$cell_typeCode)=="13"] <- "Procambium"
levels(d3$cell_type)[levels(d3$cell_typeCode)=="14"] <- "Phloem and companion cells"
levels(d3$cell_type)[levels(d3$cell_typeCode)=="15"] <- "Protoxylem"

write.csv(d3,file="table with % of cis regions.csv", row.names = F)

# Getting middle of the gene, just if needed in future

ciswind<-ciswind[complete.cases(ciswind$gbp_start),]
ciswind$gbp_mid<-NA
for (i in 1:nrow(ciswind)) {
  median<-  (ciswind$gbp_start[i]+ ciswind$gbp_end[i])/2
  ciswind$gbp_mid[i]<-median
}

head(d)

