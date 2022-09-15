



dt <- read.table("c.txt", header = T)
hist(dt$Freq, xlab="Number of Cells", main="eQTL present in number of cell types")



# Question 1: Were the associations common across different cell-types or unique for specific cell-types
FreqHits_Cell <- read.table("FreqHits_Cell_TableAsscoaitionsPerCellType1.txt", header = T) #file with all associations againts all cell types, the one with that association or not
FreqHits_Cell3<-FreqHits_Cell[FreqHits_Cell$Freq==1,] # keep only rows that had an association
FreqHits_Cell4<-as.data.frame(table(FreqHits_Cell3$Var1)) #get for each unique association, the number of times it appears in the data file, hence, 
# a count of how many cell-types had that association

FreqHits_Cell4$FreqSum<-FreqHits_Cell4$Freq 
FreqHits_Cell4$Freq<-NULL
hist(FreqHits_Cell4$FreqSum) #plot to see the frequency of associations from the ones specific to a single cell to the ones shared across all cell cell types
# Answer: many of the associations were shared across 4-5 types of cells, only a few were unique to 
# a cell type, and there were several associations shared across all cell-types


#Question 2: of all assocaitions, which cell-types had more associations?
WhichCell_TypeAll <- left_join(FreqHits_Cell3,FreqHits_Cell4, by = c("Var1")) #join the file with the count of cell types an association has to a file where you know the cell type IDs
WhichCell_TypeAllTb<-table(WhichCell_TypeAll$Var1,WhichCell_TypeAll$Var2) #Var1 means association ID and Var2 is the cell type ID
write.table(WhichCell_TypeAllTb, 'WhichCell_TypeAllTb.txt', quote = F, sep = "\t", row.names = F)
WhichCell_TypeAllTb2 <- read.table("WhichCell_TypeAllTb.txt", header = T)


table(WhichCell_TypeAll$Var2)
order(table(WhichCell_TypeAll$Var2), decreasing =T)

plot(table(WhichCell_TypeAll$Var2))
# Answer: cell-type 2, 8, 4, 12, 10, and 5 respectively (range 11040 (cell-type 5) to 15412 (cell-type 2) associations) had more associations compared 
#to the others (range 4699 to 9854 associations)
# The order of importance is: 2  8  4 12 10  5  1  6 13  7 11  9 14  3 15


# Question 3: Were the 4-5 freq cell types usually the same ones?####

library(dplyr)
library(tidyr)

FreqHits_Cell5<-FreqHits_Cell4[FreqHits_Cell4$Freq==4 | FreqHits_Cell4$Freq==5,] #filter only observations with 4 or 5 counts of cell type

names(FreqHits_Cell3)
names(FreqHits_Cell5)
FreqHits_Cell5$FreqSum<-FreqHits_Cell5$Freq
FreqHits_Cell5$Freq<-NULL
WhichCell_Type <- left_join(FreqHits_Cell3,FreqHits_Cell5, by = c("Var1"))

WhichCell_Type<-WhichCell_Type[complete.cases(WhichCell_Type$FreqSum),]

WhichCell_Type4_5FreqTb<-table(WhichCell_Type$Var1,WhichCell_Type$Var2)
write.csv(WhichCell_Type4_5FreqTb, 'WhichCell_Type4_5FreqTb.csv', row.names = F)
WhichCell_Type4_5FreqTb2 <- read.csv("WhichCell_Type4_5FreqTb.csv", header = T)
print(unique(WhichCell_Type$Var2))
table(WhichCell_Type$Var2)
order(table(WhichCell_Type$Var2), decreasing =T)

# Answer: all cell types feature in those association with frequency of 4-5
#  However, the most frequent cell-types are 2,  8, 4, 12, 10, and 5
# The whole order of importance is: 2  8  4 12 10  5  1  6 13 11  9  7 15  3 14 - very similar to the pattern in
# the whole data file.

# 2: hair cells
# 8: xylem pole pericycle
# 4: non hair cells
# 12: epidermis meristem
# 10: root cap early
# 5: endodermis







