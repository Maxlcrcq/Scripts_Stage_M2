{
  library(readxl)
  library(writexl)
  library(papeR)
  library(microseq)
  library(phangorn)
  library(muscle)
  library(msa)
  library(Biostrings)
  library(seqinr)
  library(rBLAST)
  library(jsonlite)
  library(Biostrings)
  library(biomartr)
  library(readtext)
}

####recherche sequence proteines conservees####

setwd("/chemin/WD")

liste_genomes <- list.files("/chemin/Genomes_200pb_faa", full.names = F, recursive = F)
liste_genomes <- gsub(".faa","",liste_genomes)

table_all_genes <- read.csv("gene_presence_absence.csv", sep=";")
table_conserved_genes <- table_all_genes[table_all_genes$No..isolates == 61,]

colnames(table_conserved_genes) <- gsub("_fna","",colnames(table_conserved_genes))

genes_conserve <- table_conserved_genes$Gene

table_gene <- data.frame(genes_conserve)

####lecture fichier faa####

genes_conserve <- table_conserved_genes$Gene
table_gene <- data.frame(genes_conserve)

for (j in 1:length(liste_genomes)) {
# for (j in 1:4) {
  # j <- 2
  genome <- liste_genomes[j]
  fichier_faa <- list.files(path ="/chemin/Genomes_200pb_faa",
                              pattern = genome,
                              full.names = F)
  fastaFile <- readAAStringSet(paste0("/chemin/Genomes_200pb_faa/",fichier_faa))
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  df <- data.frame(seq_name, sequence)
  liste_gene <- c()
  
  for (i in 1:length(genes_conserve)) {
    # i <- 1
    item <- table_conserved_genes[i,grep(genome, colnames(table_conserved_genes))]
    seq_prot <- df$sequence[grep(item, df$seq_name)]  
    liste_gene <- c(liste_gene, seq_prot)
  }
  
  col <- c(colnames(table_gene),genome)
  table_gene <- cbind(table_gene, liste_gene)
  colnames(table_gene) <- col
}

write_xlsx(table_gene, "/chemin/table_conserved_genes_S_uberis.xlsx")

table_gene_transpose <- data.frame(t(table_gene))
colnames(table_gene_transpose) <- table_gene_transpose[1,]
table_gene_transpose <- table_gene_transpose[-1,]

for (i in 1:70) {
  # i <- 2
  deb <- 10 * (i -1) + 1
  fin <- deb + 9
  table_prot_fusion[[i]] <- tidyr::unite(table_gene_transpose[,deb : fin], col = "union", sep = "")
  table_prot_fusion[[i]]$size <- nchar(table_prot_fusion[[i]][,1])
  table_prot_fusion[[i]]$genome <- rownames(table_prot_fusion[[i]])
}

# write_xlsx(table_prot_fusion_1, "table_prot_fusion_S_uberis_1.xlsx")
for (j in 1:70) {
  # j <- 1
  prot_fusion <- AAStringSet(c(table_prot_fusion[[j]]$union))
  names <- table_prot_fusion[[j]]$genome
  names(prot_fusion) <- table_prot_fusion[[j]]$genome
  fusion_fasta <- as.list(prot_fusion)
  names(fusion_fasta)
  write.fasta(fusion_fasta, names, paste0("/chemin/proteine_fusion_uberis_",j,".faa"))
  aln <- muscle::muscle(prot_fusion)
  AAstr = as(aln, "AAStringSet")
  noms_table <- paste0("/chemin/alignement_strepto_uberis_w_s_porcinus_",j,".faa")
  writeXStringSet(AAstr, file = noms_table)
}

prot_fusion_final <- data.frame(Souche = rownames(table_prot_fusion[[1]]), sequence = rep("XXXX",61))
for (j in 1:70) {
  aln <- readAAStringSet(paste0("/chemin/alignement_strepto_uberis_w_porcinus_",j,".faa"))
  prot_fusion_final$sequence <- paste0(prot_fusion_final$sequence, paste(aln))
}

list_prot_fusion <- as.list(prot_fusion_final$sequence)
names(list_prot_fusion) <- prot_fusion_final$Souche
write.fasta(list_prot_fusion, prot_fusion_final$Souche, paste0("/chemin/alignement_strepto_uberis_w_porcinus.faa"))

text <- readtext("/chemin/alignement_edited_no_gaps_nj.txt")

table_souches <- read_xlsx("/chemin/Sub_Casdar_subset_2011_DNAseq.xlsx", sheet = "data_Florence_with_ST")

for (i in (1:nrow(table_souches))) {
  text$text <- gsub(table_souches$Tree_node_label_1[i], table_souches$ID_souche[i], text$text)
}

write.table(text$text, "/chemin/alignement_edited_no_gaps_nj.txt", quote = F)


####statistiques ####

setwd("/chemin/Genomes_200pb_faa")

liste_genomes <- list.files("/chemin/Genomes_200pb_faa", full.names = F, recursive = F)
liste_genomes <- gsub(".faa","",liste_genomes)
db_nb_CDS <- data.frame()

for (j in 1:length(liste_genomes)) {
  # for (j in 1:4) {
  # j <- 1
  genome <- liste_genomes[j]
  fichier_faa <- list.files(path ="/chemin/Genomes_200pb_faa",
                            pattern = genome,
                            full.names = F)
  fastaFile <- readAAStringSet(paste0("/chemin/Genomes_200pb_faa/",fichier_faa))
  seq_name = c(genome,length(fastaFile))
  db_nb_CDS <- rbind(db_nb_CDS, seq_name)
  
}

colnames(db_nb_CDS) <- c("Strain", "Nb_CDS")
write_xlsx(db_nb_CDS, "//chemin/table_conserved_genes_uberis_CDS.xlsx")