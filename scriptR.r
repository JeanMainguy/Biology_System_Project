### Ouvrir les 3 fichiers
goi = read.table("GOI.txt", header=F)
goi_fur1 = read.table("GOI_FUR_1.txt", header=T)
goi_fur2 = read.table("GOI_FUR_2.txt", header=T)

### Modification pour rendre les tables plus jolies
colnames(goi) = "gene"

### Triage des tables fur pour garder les genes significatifs

goi_fur1_sorted = subset(goi_fur1, abs(log2.fold_change._Dfur.wt.) >= 2)
goi_fur2_sorted = subset(goi_fur2, abs(log2.fold_change._.Dfur.wt.) >= 2)

### Extraction de la liste des gènes (TODO)
