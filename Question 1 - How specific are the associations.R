
# Question 1: How specific are the eQTL associations

# Goal: Plot histogram of how many cell-types a hit is found
# (hit refers to a significant association between a genome position and RNA count of a particular gene)

# We are using the data stores in the RSB FlashHeart supercomputer (@130.56.33.253) throught the Terminal
#If you have the storage space, you can do this in your computer

# Go to Terminal
# Log in to Flashheart
cd /home/adamr/single_cell/eQTL_results/ # change directory
ls -ltrh #list files by date
R #go to R


###Step 1: get a list of SNPs that had at least one hit in each cell-type####

# A 'hit' is a strictly signifficant eQTL (SNP - RNA count association)
library(RSQLite)
library(data.table)
library(dplyr)
library(tidyr)

#make connection to new database
conn <- dbConnect(RSQLite::SQLite(), "results_columella.db")
hits <- dbGetQuery(conn, "SELECT * FROM GWAS_results WHERE p < 7.764071e-13")
#p is extremally conservative

snps <- dplyr::select(hits, Chr, bp, gene) #in hits, take columns Chr, bp, and gene
write.table(snps, "hit_list_columella.txt", quote = F, sep = "\t", row.names = F)
dbDisconnect(conn)

###  Cortex

#make connection to new database

conn <- dbConnect(RSQLite::SQLite(), "results_cortex.db")

hits <- dbGetQuery(conn, "SELECT * FROM GWAS_results WHERE p < 7.764071e-13")

snps <- select(hits, Chr, bp, gene)

write.table(snps, "hit_list_cortex.txt", quote = F, sep = "\t", row.names = F)

dbDisconnect(conn)


###Edodermins

#make connection to new database

conn <- dbConnect(RSQLite::SQLite(), "results_endodermis.db")

hits <- dbGetQuery(conn, "SELECT * FROM GWAS_results WHERE p < 7.764071e-13")

snps <- select(hits, Chr, bp, gene)

write.table(snps, "hit_list_endodermis.txt", quote = F, sep = "\t", row.names = F)

dbDisconnect(conn)


###Epidermis Meristem

#make connection to new database

conn <- dbConnect(RSQLite::SQLite(), "results_epidermis_meristem.db")

hits <- dbGetQuery(conn, "SELECT * FROM GWAS_results WHERE p < 7.764071e-13")

snps <- select(hits, Chr, bp, gene)

write.table(snps, "hit_list_epidermis_meristem.txt", quote = F, sep = "\t", row.names = F)

dbDisconnect(conn)


###Hair cells

#make connection to new database

conn <- dbConnect(RSQLite::SQLite(), "results_hair_cells.db")

hits <- dbGetQuery(conn, "SELECT * FROM GWAS_results WHERE p < 7.764071e-13")

snps <- select(hits, Chr, bp, gene)

write.table(snps, "hit_list_hair_cells.txt", quote = F, sep = "\t", row.names = F)

dbDisconnect(conn)


###Nonhair cells

#make connection to new database

conn <- dbConnect(RSQLite::SQLite(), "results_nonhair_cells.db")

hits <- dbGetQuery(conn, "SELECT * FROM GWAS_results WHERE p < 7.764071e-13")

snps <- select(hits, Chr, bp, gene)

write.table(snps, "hit_list_nonhair_cells.txt", quote = F, sep = "\t", row.names = F)

dbDisconnect(conn)


###Phloem and companion cells

#make connection to new database

conn <- dbConnect(RSQLite::SQLite(), "results_phloem_and_companion_cells.db")

hits <- dbGetQuery(conn, "SELECT * FROM GWAS_results WHERE p < 7.764071e-13")

snps <- select(hits, Chr, bp, gene)

write.table(snps, "hit_list_phloem_and_companion_cells.txt", quote = F, sep = "\t", row.names = F)

dbDisconnect(conn)


###Phloem Pole pericycle

#make connection to new database

conn <- dbConnect(RSQLite::SQLite(), "results_phloem_pole_pericycle.db")

hits <- dbGetQuery(conn, "SELECT * FROM GWAS_results WHERE p < 7.764071e-13")

snps <- select(hits, Chr, bp, gene)

write.table(snps, "hit_list_phloem_pole_pericycle.txt", quote = F, sep = "\t", row.names = F)

dbDisconnect(conn)


###Procabium

#make connection to new database

conn <- dbConnect(RSQLite::SQLite(), "results_procambium.db")

hits <- dbGetQuery(conn, "SELECT * FROM GWAS_results WHERE p < 7.764071e-13")

snps <- select(hits, Chr, bp, gene)

write.table(snps, "hit_list_procambium.txt", quote = F, sep = "\t", row.names = F)

dbDisconnect(conn)


###Protoxylem

#make connection to new database

conn <- dbConnect(RSQLite::SQLite(), "results_protoxylem.db")

hits <- dbGetQuery(conn, "SELECT * FROM GWAS_results WHERE p < 7.764071e-13")

snps <- select(hits, Chr, bp, gene)

write.table(snps, "hit_list_protoxylem.txt", quote = F, sep = "\t", row.names = F)

dbDisconnect(conn)


###Root cap early

#make connection to new database

conn <- dbConnect(RSQLite::SQLite(), "results_root_cap_early.db")

hits <- dbGetQuery(conn, "SELECT * FROM GWAS_results WHERE p < 7.764071e-13")

snps <- select(hits, Chr, bp, gene)

write.table(snps, "hit_list_root_cap_early.txt", quote = F, sep = "\t", row.names = F)

dbDisconnect(conn)


###Root cap lateral

#make connection to new database

conn <- dbConnect(RSQLite::SQLite(), "results_root_cap_lateral.db")

hits <- dbGetQuery(conn, "SELECT * FROM GWAS_results WHERE p < 7.764071e-13")

snps <- select(hits, Chr, bp, gene)

write.table(snps, "hit_list_root_cap_lateral.txt", quote = F, sep = "\t", row.names = F)

dbDisconnect(conn)


###Root cap tip

#make connection to new database

conn <- dbConnect(RSQLite::SQLite(), "results_root_cap_tip.db")

hits <- dbGetQuery(conn, "SELECT * FROM GWAS_results WHERE p < 7.764071e-13")

snps <- select(hits, Chr, bp, gene)

write.table(snps, "hit_list_root_cap_tip.txt", quote = F, sep = "\t", row.names = F)

dbDisconnect(conn)


###Xylem

#make connection to new database

conn <- dbConnect(RSQLite::SQLite(), "results_xylem.db")

hits <- dbGetQuery(conn, "SELECT * FROM GWAS_results WHERE p < 7.764071e-13")

snps <- select(hits, Chr, bp, gene)

write.table(snps, "hit_list_xylem.txt", quote = F, sep = "\t", row.names = F)

dbDisconnect(conn)


###Xylem pole pericycle

#make connection to new database

conn <- dbConnect(RSQLite::SQLite(), "results_xylem_pole_pericycle.db")

hits <- dbGetQuery(conn, "SELECT * FROM GWAS_results WHERE p < 7.764071e-13")

snps <- select(hits, Chr, bp, gene)

write.table(snps, "hit_list_xylem_pole_pericycle.txt", quote = F, sep = "\t", row.names = F)

dbDisconnect(conn)


###Step 2: create a list of hits study wide#####
# This puts all unique hits across all cell types together, produced this data file:
# hit_list_all.txt

library(dplyr)
library(tidyr)
library(data.table)

#read in hits list for each cell-type

h1 <- read.table("hit_list_columella.txt", header =T)
h2 <- read.table("hit_list_cortex.txt", header =T)
h3 <- read.table("hit_list_endodermis.txt", header =T)
h4 <- read.table("hit_list_epidermis_meristem.txt", header =T)
h5 <- read.table("hit_list_hair_cells.txt", header =T)
h6 <- read.table("hit_list_nonhair_cells.txt", header =T)
h7 <- read.table("hit_list_phloem_and_companion_cells.txt", header =T)
h8 <- read.table("hit_list_phloem_pole_pericycle.txt", header =T)
h9 <- read.table("hit_list_procambium.txt", header =T)
h10 <- read.table("hit_list_protoxylem.txt", header =T)
h11 <- read.table("hit_list_root_cap_early.txt", header =T)
h12 <- read.table("hit_list_root_cap_lateral.txt", header =T)
h13 <- read.table("hit_list_root_cap_tip.txt", header =T)
h14 <- read.table("hit_list_xylem.txt", header =T)
h15 <- read.table("hit_list_xylem_pole_pericycle.txt", header =T)

hits <- rbind(h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12, h13, h14, h15)
hits2 <- unique(hits)
write.table(hits2, "hit_list_all.txt", quote = F, sep = "\t", row.names = F)

###


###Step 3: extract hits from the different databases####

#This basically adds again the eQTLresults to the new file with all genome wide hits
#So after this, I need to filter the data for the desired p-value again
# Produced  file 'c.txt: full list of hits, and the frequency of cell types its occurs

install.packages("RSQLite")
library(RSQLite)
library(data.table)
library(dplyr)
library(tidyr)

#make connection to new database

conn <- dbConnect(RSQLite::SQLite(), "results_columella.db")

results <- dbReadTable(conn, "GWAS_results")

hits <- read.table("hit_list_all.txt", header = T)

filtered <- left_join(hits, results, by = c("Chr", "bp", "gene"))

write.table(filtered, "columella_all_hits.txt", quote = F, sep = "\t", row.names = F)

dbDisconnect(conn)

###

#make connection to new database

conn <- dbConnect(RSQLite::SQLite(), "results_cortex.db")

results <- dbReadTable(conn, "GWAS_results")

hits <- read.table("hit_list_all.txt", header = T)

filtered <- left_join(hits, results, by = c("Chr", "bp", "gene"))

write.table(filtered, "cortex_all_hits.txt", quote = F, sep = "\t", row.names = F)

dbDisconnect(conn)

###

#make connection to new database

conn <- dbConnect(RSQLite::SQLite(), "results_endodermis.db")

results <- dbReadTable(conn, "GWAS_results")

hits <- read.table("hit_list_all.txt", header = T)

filtered <- left_join(hits, results, by = c("Chr", "bp", "gene"))

write.table(filtered, "endodermis_all_hits.txt", quote = F, sep = "\t", row.names = F)

dbDisconnect(conn)

###

#make connection to new database

conn <- dbConnect(RSQLite::SQLite(), "results_epidermis_meristem.db")

results <- dbReadTable(conn, "GWAS_results")

hits <- read.table("hit_list_all.txt", header = T)

filtered <- left_join(hits, results, by = c("Chr", "bp", "gene"))

write.table(filtered, "epidermis_meristem_all_hits.txt", quote = F, sep = "\t", row.names = F)

dbDisconnect(conn)

###

#make connection to new database

conn <- dbConnect(RSQLite::SQLite(), "results_hair_cells.db")

results <- dbReadTable(conn, "GWAS_results")

hits <- read.table("hit_list_all.txt", header = T)

filtered <- left_join(hits, results, by = c("Chr", "bp", "gene"))

write.table(filtered, "hair_cells_all_hits.txt", quote = F, sep = "\t", row.names = F)

dbDisconnect(conn)

###

#make connection to new database

conn <- dbConnect(RSQLite::SQLite(), "results_nonhair_cells.db")

results <- dbReadTable(conn, "GWAS_results")

hits <- read.table("hit_list_all.txt", header = T)

filtered <- left_join(hits, results, by = c("Chr", "bp", "gene"))

write.table(filtered, "nonhair_cells_all_hits.txt", quote = F, sep = "\t", row.names = F)

dbDisconnect(conn)

###

#make connection to new database

conn <- dbConnect(RSQLite::SQLite(), "results_phloem_and_companion_cells.db")

results <- dbReadTable(conn, "GWAS_results")

hits <- read.table("hit_list_all.txt", header = T)

filtered <- left_join(hits, results, by = c("Chr", "bp", "gene"))

write.table(filtered, "phloem_and_companion_cells_all_hits.txt", quote = F, sep = "\t", row.names = F)

dbDisconnect(conn)

###

#make connection to new database

conn <- dbConnect(RSQLite::SQLite(), "results_phloem_pole_pericycle.db")

results <- dbReadTable(conn, "GWAS_results")

hits <- read.table("hit_list_all.txt", header = T)

filtered <- left_join(hits, results, by = c("Chr", "bp", "gene"))

write.table(filtered, "phloem_pole_pericycle_all_hits.txt", quote = F, sep = "\t", row.names = F)

dbDisconnect(conn)

###

#make connection to new database

conn <- dbConnect(RSQLite::SQLite(), "results_procambium.db")

results <- dbReadTable(conn, "GWAS_results")

hits <- read.table("hit_list_all.txt", header = T)

filtered <- left_join(hits, results, by = c("Chr", "bp", "gene"))

write.table(filtered, "procambium_all_hits.txt", quote = F, sep = "\t", row.names = F)

dbDisconnect(conn)

###

#make connection to new database

conn <- dbConnect(RSQLite::SQLite(), "results_protoxylem.db")

results <- dbReadTable(conn, "GWAS_results")

hits <- read.table("hit_list_all.txt", header = T)

filtered <- left_join(hits, results, by = c("Chr", "bp", "gene"))

write.table(filtered, "protoxylem_all_hits.txt", quote = F, sep = "\t", row.names = F)

dbDisconnect(conn)

###

#make connection to new database

conn <- dbConnect(RSQLite::SQLite(), "results_root_cap_early.db")

results <- dbReadTable(conn, "GWAS_results")

hits <- read.table("hit_list_all.txt", header = T)

filtered <- left_join(hits, results, by = c("Chr", "bp", "gene"))

write.table(filtered, "root_cap_early_all_hits.txt", quote = F, sep = "\t", row.names = F)

dbDisconnect(conn)

###

#make connection to new database

conn <- dbConnect(RSQLite::SQLite(), "results_root_cap_lateral.db")

results <- dbReadTable(conn, "GWAS_results")

hits <- read.table("hit_list_all.txt", header = T)

filtered <- left_join(hits, results, by = c("Chr", "bp", "gene"))

write.table(filtered, "root_cap_lateral_all_hits.txt", quote = F, sep = "\t", row.names = F)

dbDisconnect(conn)

###

#make connection to new database

conn <- dbConnect(RSQLite::SQLite(), "results_root_cap_tip.db")

results <- dbReadTable(conn, "GWAS_results")

hits <- read.table("hit_list_all.txt", header = T)

filtered <- left_join(hits, results, by = c("Chr", "bp", "gene"))

write.table(filtered, "root_cap_tip_all_hits.txt", quote = F, sep = "\t", row.names = F)

dbDisconnect(conn)

###

#make connection to new database

conn <- dbConnect(RSQLite::SQLite(), "results_xylem.db")

results <- dbReadTable(conn, "GWAS_results")

hits <- read.table("hit_list_all.txt", header = T)

filtered <- left_join(hits, results, by = c("Chr", "bp", "gene"))

write.table(filtered, "xylem_all_hits.txt", quote = F, sep = "\t", row.names = F)

dbDisconnect(conn)

###

#make connection to new database

conn <- dbConnect(RSQLite::SQLite(), "results_xylem_pole_pericycle.db")

results <- dbReadTable(conn, "GWAS_results")

hits <- read.table("hit_list_all.txt", header = T)

filtered <- left_join(hits, results, by = c("Chr", "bp", "gene"))

write.table(filtered, "xylem_pole_pericycle_all_hits.txt", quote = F, sep = "\t", row.names = F)

dbDisconnect(conn)

###

# Bind it all together

h1 <- read.table("columella_all_hits.txt", header =T)
h2 <- read.table("cortex_all_hits.txt", header =T)
h3 <- read.table("endodermis_all_hits.txt", header =T)
h4 <- read.table("epidermis_meristem_all_hits.txt", header =T)
h5 <- read.table("hair_cells_all_hits.txt", header =T)
h6 <- read.table("nonhair_cells_all_hits.txt", header =T)
h7 <- read.table("phloem_and_companion_cells_all_hits.txt", header =T)
h8 <- read.table("phloem_pole_pericycle_all_hits.txt", header =T)
h9 <- read.table("procambium_all_hits.txt", header =T)
h10 <- read.table("protoxylem_all_hits.txt", header =T)
h11 <- read.table("root_cap_early_all_hits.txt", header =T)
h12 <- read.table("root_cap_lateral_all_hits.txt", header =T)
h13 <- read.table("root_cap_tip_all_hits.txt", header =T)
h14 <- read.table("xylem_all_hits.txt", header =T)
h15 <- read.table("xylem_pole_pericycle_all_hits.txt", header =T)

dat <- rbind(h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12, h13, h14, h15)

dat2 <- filter(dat, p < 7.764071e-13) #keep only signifficant hits for each cell type

dat2 <- tidyr::unite(dat2, hit, 1,2,3, remove = TRUE) #creates a
# new column call 'hit' that paste together the info in columns 1,2, and 3
# remove=TRUE remove the input columns
# sep = "_" default


dat3 <- as.data.frame(table(dat2$hit)) #get the frequency of each hit
colnames(dat3)[1] <- "hit"

write.table(dat3, "c.txt", quote = F, row.names = F, sep = "\t")

# file 'c.txt" has a full list of hits, and the frequency of cell types its occurs

dt <- read.table("c.txt", header = T)

library(ggplot2)
ggplot(dt, aes(x=Freq)) +
  geom_histogram(position="identity", alpha=0.5,bins = 14,color = 'darkgrey', fill = NA)+ #,color = 'black', fill = NA
  geom_density(alpha=0.2)+
  geom_vline(data=dt, aes(xintercept=median(Freq)),color="purple",
             linetype="dashed")+
  geom_vline(data=dt, aes(xintercept=mean(Freq),color="grey"),
             linetype="dashed", show.legend=TRUE)+
  scale_color_manual(name = "Statistics", values = c(Median = "purple", Mean = "grey"))+ #adds legend
  stat_bin(binwidth=1, geom="text", aes(label=after_stat(count)), vjust=0)+ #adds txt/freq on each bar
  scale_x_continuous(breaks = seq(1,14, by = 1)) + #adds all values of x
  labs(title="Histogram of how many cell-types a hit is found",x="Count of cell-types", y = "Frequency")+
  theme_classic() 
#The logic here is: 14,843 hits (very significant eQTLs) were only found in one type of cells..

