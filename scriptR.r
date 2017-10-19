### Installer les packages necessaire

library("reshape2") ### pour fonction melt dans la correlation

### Ouvrir les fichiers
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
write.table(gene_list, "gene_list.txt", row.names=F, col.names=F, quote=F)

expression_table = subset(expression_table_all, expression_table_all$IDENTIFIER %in% gene_list)
rownames(expression_table) = expression_table[,"IDENTIFIER"]
expression_table = expression_table[,3:20]
expression_table = t(expression_table)

### Corrélation

correlation =cor(expression_table)
list_correlation = melt(correlation)

#On filtre les résultats de corrélation :
list_correlation=list_correlation[abs(list_correlation["value"])>0.99,]

#On élimine les doublons
list_correlation["Alphabétique"]<-as.character(list_correlation[,"Var1"])<as.character(list_correlation[,"Var2"])
list_correlation=list_correlation[list_correlation[,4]==TRUE,]
list_correlation=list_correlation[,0:3]

#Préparation du .sif
colnames(list_correlation) = c("node1", "node2", "lien")
list_correlation[,3] = "correlation"
list_correlation = list_correlation[, c(1,3,2)]

#Ajout des sRNA

#lecture de la table
interac_table = read.table("sRNA_interaction.txt",header = T)
#filtrer pour garder seulement les genes d'interet
interac_table_filter = subset(interac_table, interac_table$CIBLES %in% unlist(list_correlation[,1])) 
interac_table_filter = interac_table_filter[, c(1,3,2)]
#Fusion des tables
colnames(interac_table_filter) = c("node1", "lien", "node2")
liste_arn_gene= merge(list_correlation, interac_table_filter, by = intersect(names(list_correlation), names(interac_table_filter)), all.x = T, all.y = T)

#On enregistre dans un fichier

write.table(liste_arn_gene, "liste99.sif", row.names=F, col.names=F, quote=F)
