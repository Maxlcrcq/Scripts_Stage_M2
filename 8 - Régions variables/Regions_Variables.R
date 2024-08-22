# BiocManager::install(c("bios2mds", "msa", "seqrflp", "spgs", "msa", "phangorn", "microseq", "rgff", "rtracklayer", "spliceR","rBLAST"))

###Chargement des librairies####
{
  library(readxl)
  library(writexl)
  library(papeR)
  library(microseq)
  library(phangorn)
  library(muscle)
  library(msa)
  library(Biostrings)
  library(rBLAST)
  library(jsonlite)
  library(biomartr)
  library(bios2mds)
  library(dplyr)
  library(rtracklayer)
  library(tidyr)
  library(spgs) # to reverse complement DNA sequence
  library(forstringr) # to use str_split
  library(DECIPHER)
  library(seqinr) # to manipulate DNA and AA sequences
  library(ggplot2)
  library(phytools) # to use force.ultrametric
  library(ggpmisc)
  library(stringr)
  library(patchwork)
}

setwd("/home/maximeleclercq/Bureau/Clean/7 - Analyses_Regions_variables")

####recherche sequence proteines conservees####

liste_genomes <- list.files("/home/maximeleclercq/Bureau/Clean/3 - Prokka/FNA",
                            pattern = ".fna$", full.names = FALSE, recursive = FALSE)
liste_genomes <- gsub(".fna", "", liste_genomes)

#suppression du génome de S. porcinus de la liste des gènomes
liste_genomes <- liste_genomes[liste_genomes != "Porcinus"]

#chargement de la table des ID des gènes conservés après Roary
table_all_genes <- read.csv("/home/maximeleclercq/Bureau/Clean/4- Roary/Roary_sans_porcinus/gene_presence_absence.csv",sep = ";")

# table_all_genes <- table_all_genes[,-grep("901553735", colnames(table_all_genes))]
table_conserved_genes <- table_all_genes[table_all_genes$No..isolates == 88,]
colnames(table_conserved_genes) <- gsub("_fna", "", colnames(table_conserved_genes))
table_conserved_genes <- table_conserved_genes[table_conserved_genes$Gene != "group_3",]

genes_conserve <- table_conserved_genes$Gene

#créer une table avec toutes les séquences de tous les gènes conservés
table_gene <- data.frame(matrix(ncol = length(genes_conserve) + 1, nrow = 0))
colnames(table_gene) <- c("genome", genes_conserve)

table_gene <- data.frame()

for (j in 1:length(liste_genomes)) {
  # for (j in 1:4) {
  # j <- 85
  genome <- liste_genomes[j]
  genome_gff <- data.frame(readGFF(paste0("/home/maximeleclercq/Bureau/Clean/3 - Prokka/GFF/",genome,".gff")))
  fichier_fna <- list.files(path ="/home/maximeleclercq/Bureau/Clean/3 - Prokka/FNA/",
                            pattern = paste0(genome,".fna$"),
                            full.names = F)
  fastaFile <- readDNAStringSet(paste0("/home/maximeleclercq/Bureau/Clean/3 - Prokka/FNA/",fichier_fna))
  names(fastaFile)
  length(fastaFile)
  fastaFile@ranges@width
  gcf <- as.character(grep("GCF", genome))
  ja <- as.character(grep("JA", genome))
  
  length(gcf)
  if (length(gcf) + length(ja) > 0 ) {
    gffFile1 <- paste0("/home/maximeleclercq/Bureau/Clean/3 - Prokka/GFF/",genome,".gff")
    table_gff_test <- read.delim(gffFile1, header=F)
    col_titre <- grep("##", table_gff_test[,1])
    names_contigs <- table_gff_test[col_titre,1]
    names_contigs <- names_contigs[-c(1, length(names_contigs))]
    names_contigs <- substr(names_contigs, 19, nchar(names_contigs))
    # unlist(str_split(names_contigs," ", 3))
    
    names_test <- unlist(gregexpr(" ", names_contigs))[seq(1,80, by = 2)]
    names_contigs <- substr(names_contigs, 1,names_test -1)
    seq_name = names_contigs
  } else {
    seq_name = names(fastaFile)
  }
  
  # seq_name = names_contigs
  # seq_name = names(fastaFile)
  seq_name
  sequence = paste(fastaFile)
  df <- data.frame(seq_name, sequence)
  
  genome <- paste0("^", genome)
  
  liste_gene <- c()
  
  for (i in 1:length(genes_conserve)) {
    # i <- 999
    item <- table_conserved_genes[i,grep(genome, colnames(table_conserved_genes))]
    if (unlist(gregexpr("\\t", item)) >1) {
      item <- substr(item, 1,unlist(gregexpr("\\t", item))[1] -1)
    }
    item <- paste0(item, "$")
    contig <- genome_gff$seqid[grep(item, genome_gff$ID)]
    start <- genome_gff$start[grep(item, genome_gff$ID)]
    end <- genome_gff$end[grep(item, genome_gff$ID)]
    strand <- genome_gff$strand[grep(item, genome_gff$ID)]
    sequence <- subset(df, df$seq_name == contig)
    if (strand == "+") {
      sequence <- substr(sequence[1,2], start, end)
    } else {
      sequence <- substr(sequence[1,2], start, end)
      sequence <- toupper(spgs::reverseComplement(sequence, content=c("dna")))
    }
    liste_gene <- c(liste_gene, sequence)
  }
  
  
  table_gene <- rbind(table_gene, liste_gene)
}

colnames(table_gene) <- genes_conserve
rownames(table_gene) <- liste_genomes


write_xlsx(table_gene, "table_gene.xlsx")

table_gene <- data.frame(read_xlsx("table_gene.xlsx"))

#créer un table avec les gènes conservés alignés et les gaps éliminés
grep("axe1", colnames(table_gene))
table_genes_trimmed <- data.frame(Header = liste_genomes)

for (i in 1:length(genes_conserve)) {
  # i <- 1
  sequences_input <- DNAStringSet(table_gene[,i])
  names(sequences_input) <- liste_genomes
  aln <- msa(inputSeqs = sequences_input,
             method=c("Muscle"))

  writeXStringSet(unmasked(aln), file="test.fasta")

  aln_fasta <- readFasta("test.fasta")
  aln_trimmed <- data.frame(msaTrim(aln_fasta, gap.end = 0, gap.mid = 0))

  col <- colnames(table_genes_trimmed)
  table_genes_trimmed <- merge(table_genes_trimmed, aln_trimmed, by="Header")
  colnames(table_genes_trimmed) <- c(col, genes_conserve[i])
}

#nom des colonnes avec des NA
num_col_NA <- colnames(table_genes_trimmed[colSums(is.na(table_genes_trimmed)) != 0])
# suppression de cette colonne
table_genes_trimmed_wo_NA <- table_genes_trimmed[, !is.element(colnames(table_genes_trimmed), num_col_NA)]

write_xlsx(table_genes_trimmed_wo_NA, "table_genes_trimmed_wo_S_porcinus.xlsx")

table_genes_trimmed <- read_xlsx("table_genes_trimmed_wo_S_porcinus.xlsx")

####travail sur les alignements générés####
table_genes_trimmed <- read_xlsx("table_genes_trimmed_wo_S_porcinus.xlsx")
grep("axe1-6A", colnames(table_genes_trimmed))
grep("group_3$", colnames(table_genes_trimmed))

table_all_genes <- read.csv("/home/maximeleclercq/Bureau/Clean/4- Roary/Roary_sans_porcinus/gene_presence_absence.csv", sep=";")

table_conserved_genes <- table_all_genes[table_all_genes$No..isolates == 88,]
table_conserved_genes <- table_conserved_genes[table_conserved_genes$Gene!="group_3",]

grep("axe1-6A", table_conserved_genes$Gene)
grep("group_3$", table_conserved_genes$Gene)
colnames(table_conserved_genes) <- gsub("_fna","",colnames(table_conserved_genes))

####récupérer l'ordre gènes de 0140J####

#GCF_000009545.1
genome_0140J_gff <- data.frame(readGFF(paste0("/home/maximeleclercq/Bureau/Clean/3 - Prokka/GFF/GCF_000009545.gff")))
genome_0140J_gff_gene <- genome_0140J_gff[genome_0140J_gff$type == "gene",]
liste_gene_0140J_ordonne <- data.frame(Key=seq(1:nrow(genome_0140J_gff_gene)), GCF_000009545_1_ASM954v1 = genome_0140J_gff_gene$locus_tag, start = genome_0140J_gff_gene$start, end = genome_0140J_gff_gene$end)

#récupérer les identifiants des gènes conservés de O140J
colonne_0140J <- grep("09545", colnames(table_conserved_genes))
ordre_genes_conserves_0140J <- data.frame(gene_ID=table_conserved_genes$Gene, GCF_000009545_1_ASM954v1 = table_conserved_genes[, colonne_0140J])
ordre_genes_conserves_0140J <- merge(liste_gene_0140J_ordonne, ordre_genes_conserves_0140J, by="GCF_000009545_1_ASM954v1")


#réduire la table des séquences des gènes à ceux présents dans O140J
table_genes_trimmed_ordonne <- table_genes_trimmed[,c("Header", ordre_genes_conserves_0140J$gene_ID)]
taille_genes <- nchar(table_genes_trimmed_ordonne[1, -1])


ordre_genes_conserves_0140J <- cbind(ordre_genes_conserves_0140J, taille_genes)

ordre_genes_conserves_0140J$end_concatenated_genome <- cumsum(ordre_genes_conserves_0140J$taille_genes)
ordre_genes_conserves_0140J$start_concatenated_genome <- ordre_genes_conserves_0140J$end_concatenated_genome - ordre_genes_conserves_0140J$taille_genes + 1



table_genes_fusionne <- table_genes_trimmed_ordonne %>% 
  unite(genome, c(2:1407), sep = "", remove = TRUE)

write.csv2(table_genes_fusionne,"table_genes_fusionne.csv")

write.fasta(sequences = as.list(table_genes_fusionne$genome), names = table_genes_fusionne$Header, file.out = "core_genome_ordonne_wo_porcinus.fna")


#Creation table pour stockage des données de likelihood
table_likelihood <- data.frame()

nchar(table_genes_fusionne[2,2])
(nchar(table_genes_fusionne[1,2]) - 5000) / 250
# 2125 * 30 / 3600

{
debut <- Sys.time() 
for (i in (0:5072)) {
  # i <- 0
  deb <- 250 * i
  fin <- deb + 5000
  subtable1 <- data.frame(Header = table_genes_fusionne$Header, Sequence = substr(table_genes_fusionne$genome, deb ,fin))

  #calcul du nombre de sites polymorphes
  subtable2 <- data.frame(do.call(rbind, strsplit(subtable1$Sequence, "")))
  test <- apply(subtable2 , 2, unique)
  polymorphic <- as.vector(unlist(lapply(test, length)))
  number_polymorphism <- length(polymorphic[ polymorphic > 1])
    
  
  #
  
  write.fasta(sequences = as.list(subtable1$Sequence), names = subtable1$Header, file.out = "subtable.fna")
  # subset_fasta <- readDNAStringSet("subtable.fna")
  
  aln <- read.alignment("subtable.fna", format = "fasta")
  data <- as.phyDat(aln, type = "DNA")
  mt <- modelTest(data, model=c("HKY"), 
                  control = pml.control(trace = 0))
  
  pml_results <- as.pml(mt, "HKY")
  
  table_likelihood <- rbind(table_likelihood, c(number_polymorphism,pml_results$logLik ))
  }
fin_script <- Sys.time()
duree <- fin_script - debut
}

colnames(table_likelihood) <- c("polymorphism","likelihood")

#verification que write.phyDat(data, file = "reverse.fna", format = "fasta")

write_xlsx(table_likelihood, "table_likelihood_sans_porcinus.xlsx")

table_likelihood <- read_xlsx("table_likelihood_sans_porcinus.xlsx")

#### Graphs ####

# Charger les données
table_likelihood <- read_xlsx("table_likelihood_sans_porcinus.xlsx")
table_genes_fusionne <- read.csv("table_genes_fusionne.csv", sep = ";")

# Extraire la séquence de la colonne genome
sequence <- as.character(table_genes_fusionne$genome[1])

# Fonction pour calculer le pourcentage de GC
calculate_gc_content <- function(sequence) {
  gc_count <- sum(str_count(sequence, "[GCgc]"))
  total_count <- nchar(sequence)
  return((gc_count / total_count) * 100)
}

# Initialiser les paramètres de la fenêtre glissante
window_size <- 5000
step_size <- 250
positions <- seq(1, nchar(sequence) - window_size + 1, by = step_size)

# Initialiser la colonne GC%
gc_values <- sapply(positions, function(pos) {
  window_seq <- substr(sequence, pos, pos + window_size - 1)
  calculate_gc_content(window_seq)
})

# Ajouter les valeurs de GC% et les positions à la table_likelihood
table_likelihood$gc_content <- gc_values
table_likelihood$position <- positions

# Ajuster les modèles et calculer les z-scores
model <- lm(table_likelihood$likelihood ~ table_likelihood$polymorphism)
table_likelihood$estimated <- model$coefficients[1] + model$coefficients[2] * table_likelihood$polymorphism
table_likelihood$deviation <- table_likelihood$likelihood - table_likelihood$estimated
valeur_reduite <- sd(table_likelihood$likelihood)
table_likelihood$z_score <- table_likelihood$deviation / valeur_reduite

# Définir le seuil
seuil <- -1.75

# Trouver les segments où le z-score est inférieur au seuil
sub_positions <- which(table_likelihood$z_score < seuil)

# Identifier les débuts et fins des segments continus
segments <- split(sub_positions, cumsum(c(1, diff(sub_positions)) != 1))
positions_interets <- sapply(segments, function(segment) mean(table_likelihood$position[segment]))

# Créer le graphe des z-scores avec des lignes verticales
G1 <- ggplot(table_likelihood, aes(x = position, y = z_score)) +
  geom_point(color = "blue") +
  geom_vline(xintercept = positions_interets, linetype = "dashed", color = "black") +
  theme(
    axis.text.x = element_text(angle = 90),
    panel.grid = element_line(colour = "lightgrey")
  ) +
  scale_x_continuous(breaks = seq(0, max(table_likelihood$position), 25000)) +
  labs(y = "Z-Score")

# Créer le graphe des GC%
G2 <- ggplot(table_likelihood, aes(x = position, y = gc_content)) +
  geom_line(color = "red") +
  geom_hline(yintercept = 36.45, linetype = "dashed", color = "darkorange4") +
  geom_vline(xintercept = positions_interets, linetype = "dashed", color = "black") +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_line(colour = "lightgrey"),
    legend.position = "none"
  ) +
  scale_x_continuous(breaks = seq(0, max(table_likelihood$position), 25000)) +
  labs(y = "GC%", x = "Position")

# Fusionner les deux graphes
Graph_combine <- G2 / G1

# Afficher le graphe combiné
print(Graph_combine)

# Analyse des minima
table_minima <- data.frame(num = c(1:4), position = c(87000, 150000, 715000, 1137000))

for (i in 1:4) {
  mini <- table_minima$position[i] - 10000
  maxi <- table_minima$position[i] + 10000
  table_likelihood_subset <- subset(table_likelihood, table_likelihood$position > mini & table_likelihood$position < maxi)
  
  g <- ggplot(table_likelihood_subset, aes(x = position)) +
    geom_point(aes(y = deviation)) +
    theme(axis.text.x = element_text(angle = 90),
          panel.grid = element_line(colour = "lightgrey")) +
    scale_x_continuous(breaks = seq(mini, maxi, 1000))
  
  print(g)
}

mini <- 650000
maxi <- 770000
table_likelihood_subset <- subset(table_likelihood, table_likelihood$position > mini & table_likelihood$position < maxi)

g <- ggplot(table_likelihood_subset, aes(x = position)) +
  geom_point(aes(y = deviation)) +
  theme(axis.text.x = element_text(angle = 90),
        panel.grid = element_line(colour = "lightgrey")) +
  scale_x_continuous(breaks = seq(mini, maxi, 5000))

print(g)

print(positions_interets)

#### Analyse Région Spécifique ####
#definition region####

gene1 <- c("dnaG")
gene2 <- c("infC")

df_region <- data.frame(Strain=c(""), Sequence = c(""))

for (j in 1:length(liste_genomes)) {
  # for (j in 1:4) {
  # j <- 2
  genome <- liste_genomes[j]
  
  ID_gene1_souche_etudiee <- table_conserved_genes[grep(gene1, table_conserved_genes$Gene), grep(genome, colnames(table_conserved_genes))]

  ID_gene2_souche_etudiee <- table_conserved_genes[grep(gene2, table_conserved_genes$Gene), grep(genome, colnames(table_conserved_genes))]
  
  gfffile3 <- paste0("/home/maximeleclercq/Bureau/Clean/3 - Prokka/GFF/",genome,".gff")
  genome_gff <- data.frame(readGFF(gfffile3))
  genome_gff <- genome_gff[genome_gff$type=="CDS",]
  
  contig <- as.character(genome_gff$seqid[grep(ID_gene1_souche_etudiee, genome_gff$ID)])
  
  contig <- gsub("gnl\\|.+\\|","",contig)
  
  region_gene1_start <- genome_gff$start[grep(ID_gene1_souche_etudiee, genome_gff$ID)]
  region_gene1_end <- genome_gff$end[grep(ID_gene1_souche_etudiee, genome_gff$ID)]
  strand_gene1 <- genome_gff$strand[grep(ID_gene1_souche_etudiee, genome_gff$ID)]
  
  region_gene2_start <- genome_gff$start[grep(ID_gene2_souche_etudiee, genome_gff$ID)]
  region_gene2_end <- genome_gff$end[grep(ID_gene2_souche_etudiee, genome_gff$ID)]
  
  region_start <- min(region_gene1_start,
                      region_gene1_end,
                      region_gene2_start,
                      region_gene2_end)
  
  region_end <- max(region_gene1_start,
                    region_gene1_end,
                    region_gene2_start,
                    region_gene2_end)
  
  ###
  
  table_gff_test_2 <- read.delim(gfffile3, header=F)
  col_titre <- grep("##", table_gff_test_2[,1])
  names_contigs <- table_gff_test_2[col_titre,1]
  names_contigs <- names_contigs[-c(1, length(names_contigs))]
  names_contigs <- substr(names_contigs, 19, nchar(names_contigs))
  names_test <- unlist(gregexpr(" ", names_contigs))[seq(1,80, by = 2)]
  names_contigs <- substr(names_contigs, 1,names_test -1)
  
  contig <- paste0(contig, "$")
  num_contig <- grep(contig, names_contigs)
  
  ###
  fichier_fna <- list.files(path ="/home/maximeleclercq/Bureau/Clean/3 - Prokka/FNA/",
                            pattern = paste0(genome,".fna$"),
                            full.names = F)
  fastaFile <- readDNAStringSet(paste0("/home/maximeleclercq/Bureau/Clean/3 - Prokka/FNA/",fichier_fna))
  
  

  DNA_seq <- substr(as.character(fastaFile[[num_contig]]), region_start, region_end)
  
  if(strand_gene1 == "-") {
    DNA_seq <- toupper(spgs::reverseComplement(DNA_seq, content=c("dna")))
  }
  
  df_region <- rbind(df_region, c(genome, DNA_seq))
}

df_region <- df_region[-1,]

write.csv2(df_region, "df_region.csv")

table_df <- read.csv("df_region.csv", sep=";")

#_####
#DECIPHER####
# 

table_df <- read.csv("df_region.csv", sep=";")

sequences_input <- DNAStringSet(table_df[,3])
names(sequences_input) <- table_df[,2]

table_df$longueur_region <- nchar(table_df[,3])

library(ggplot2)
ggplot(table_df, aes(x=longueur_region)) +
  geom_histogram()


DNA <- AlignSeqs(sequences_input)
names(DNA)
table_DNA_aln_DECIPHER <- data.frame(sequence = as.character(DNA))

write.fasta(sequences = as.list(table_DNA_aln_DECIPHER$sequence), names = names(DNA), file.out = "DECIPHER_aln_DNA.fna")