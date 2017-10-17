### Installer les packages necessaire

library("reshape2") ### pour fonction melt dans la correlation

### Ouvrir les 3 fichiers
goi = read.table("GOI.txt", header=F)
goi_fur1 = read.table("GOI_FUR_1.txt", header=T)
goi_fur2 = read.table("GOI_FUR_2.txt", header=T)
expression_table_all = read.table("P_EXPR_antiobios.txt", header=T)

### Modification pour rendre les tables plus jolies
colnames(goi) = "gene"

### Triage des tables fur pour garder les genes significatifs

goi_fur1_sorted = subset(goi_fur1, abs(log2.fold_change._Dfur.wt.) >= 2)
goi_fur2_sorted = subset(goi_fur2, abs(log2.fold_change._.Dfur.wt.) >= 2)

### Extraction de la liste des gènes
gene_list1 = as.vector(goi_fur1_sorted[,"Gene"])
gene_list2 = as.vector(goi_fur2_sorted[,"Gene"])
gene_list3 = as.vector(goi[,"gene"])
gene_list = c(gene_list1, gene_list2, gene_list3)

expression_table = subset(expression_table_all, expression_table_all$IDENTIFIER %in% gene_list)
rownames(expression_table) = expression_table[,"IDENTIFIER"]
expression_table = expression_table[,3:20]
expression_table = t(expression_table)

### Corrélation

correlation =cor(expression_table)
list_correlation = melt(correlation)

#On filtre les résultats de corrélation :
list_correlation=list_correlation[abs(list_correlation["value"])>0.95,]

#On élimine les doublons
list_correlation["Alphabétique"]<-as.character(list_correlation[,"Var1"])<as.character(list_correlation[,"Var2"])
list_correlation=list_correlation[list_correlation[,4]==TRUE,]
list_correlation=list_correlation[,0:2]
colnames(list_correlation) = c("gene1", "gene2")

#On enregistre dans un fichier
write.csv(list_correlation, "liste.csv", row.names=F)

#sARN 
#lecture de la table
interac_table = read.table("sRNA_interaction.txt",header = T)
#filtrer pour garder seulement les genes d'interet
interac_table_filter = subset(interac_table, interac_table$CIBLES %in% unlist(list_correlation)) 
write.csv(interac_table_filter, "liste_sRNA.csv")