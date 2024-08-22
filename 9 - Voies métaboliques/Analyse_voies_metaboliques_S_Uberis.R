library(mixOmics)
library(WallomicsData)
library(cowplot)
library(pheatmap)
library(viridis)

setwd("/home/maximeleclercq/Bureau/Clean/9 - Voies Métaboliques")

### Agréger le fichier d'origine ###
# Lire le CSV
df <- read.csv2("Rapport_Voies_Metaboliques_V2.csv", sep = "\t", dec = ".")

# Copier les valeurs de Pathway.Organism dans Super.Pathway si Super.Pathway est vide
df$Super.Pathway <- ifelse(df$Super.Pathway == "", df$Pathway.Organism, df$Super.Pathway)

# Aggréger
df_sans_sp <- df[, !grepl("Pathway.Organism", colnames(df))]
df_agg <- aggregate(df_sans_sp[,-1], by=list(Super.Pathway = df_sans_sp$Super.Pathway),FUN = mean)

# Enregistrer le nouveau CSV agrégé
write.csv(df_agg, "Rapport_Voies_Metaboliques_V2_agrege.csv", row.names = FALSE)

### Création des heatmaps ###
# Lire les CSV agrégés et originaux
heatmap1 <- read.csv2("Rapport_Voies_Metaboliques_V2_agrege.csv", sep = ",", dec = ".")
heatmap2 <- read.csv2("Rapport_Voies_Metaboliques_V2.csv", sep = "\t", dec = ".")

# Créer les matrices de données pour les heatmaps
matrix1 <- as.matrix(heatmap1[,-1])
matrix2 <- as.matrix(heatmap2[,-c(1,2)])

# Définir les noms des lignes
rownames(matrix1) <- heatmap1$Super.Pathway
rownames(matrix2) <- heatmap2$Pathway.Organism

# Définir les annotations
data_annotations <- read.csv2("Mammite-Pathogènes_V2.csv", sep = "\t")

# Créer le data frame des annotations
annotation_col <- data.frame(
  Mammite = factor(data_annotations$Mammite),
  Pathogène_Homme = factor(data_annotations$Homme),
  Espèce = factor(data_annotations$Espèce)
)

rownames(annotation_col) <- data_annotations$Souches

# Définir les couleurs des annotations
annotation_colors <- list(
  Mammite = c("Oui" = "blue", "Non" = "lightgrey"),
  Pathogène_Homme = c("Oui" = "darkorange", "Non" = "#B80F0A"),
  Espèce = c("Streptococcus agalactiae" = "green",
             "Streptococcus canis" = "#fbfeed",
             "Streptococcus dysgalactiae" = "grey",
             "Streptococcus parauberis" = "red",
             "Streptococcus pyogenes" = "turquoise",
             "Streptococcus uberis" = "darkorange")
)

# Heatmap pour matrix2 et récupération de l'ordre des colonnes
p  <- pheatmap(matrix2,
               annotation_col = annotation_col,
               show_colnames = F,
               annotation_colors = annotation_colors,
               color = viridis_pal(option = "B", direction = 1)(256),
               main = "Heatmap non agrégée des voies métaboliques des différents streptocoques"
)

# Extraction de l'ordre des colonnes à partir de la heatmap de matrix2
column_order <- colnames(matrix2)[p$tree_col$order]

# Réorganiser les colonnes de matrix1 selon l'ordre de matrix2
matrix1 <- matrix1[, column_order, drop=FALSE]

# Heatmap pour matrix1 avec l'ordre des colonnes de matrix2
pheatmap(matrix1,
         annotation_col = annotation_col,
         show_colnames = FALSE,
         cluster_cols = FALSE,
         annotation_colors = annotation_colors,
         color = viridis_pal(option = "B", direction = 1)(256),
         main = "Heatmap agrégée des voies métaboliques des différents streptocoques",
         fontsize = 14,
         fontsize_row = 14,
)

table(annotation_col$Espèce)


# Filtrer les données pour ne garder que les colonnes Streptococcus.uberis
matrix2_uberis <- matrix2[, grep("Streptococcus.uberis", colnames(matrix2))]
matrix1_uberis <- matrix1[, grep("Streptococcus.uberis", colnames(matrix1))]

# Créer les heatmaps pour les Streptococcus.uberis
pheatmap(matrix2_uberis,
         show_colnames = FALSE,
         cluster_cols = T,
         color = viridis_pal(option = "B", direction = 1)(256),
         main = "Heatmap non agrégée des voies métaboliques de Streptococcus uberis"
)

pheatmap(matrix1_uberis,
         show_colnames = FALSE,
         cluster_cols = T,
         color = viridis_pal(option = "B", direction = 1)(256),
         main = "Heatmap agrégée des voies métaboliques de Streptococcus uberis",
         fontsize = 14,
         fontsize_row = 14,
)


table_complete <- read.csv2("Rapport_Voies_Metaboliques_V2.csv", sep = "\t", dec = ".")
Table_type <- read.csv2("Clinical type Uberis.csv", sep = "\t", dec = ".")

# Filtrer les éléments contenant des NA dans les types cliniques
Table_type <- na.omit(Table_type)

Table_type$Type_Clinique2 <- Table_type$Type_Clinique
Table_type$Type_Clinique2[Table_type$Type_Clinique == "Clinique" |Table_type$Type_Clinique == "Mixte"] <- "Autre"

Table_type$Type_Clinique3 <- Table_type$Type_Clinique
Table_type$Type_Clinique3[Table_type$Type_Clinique == "Subclinique" |Table_type$Type_Clinique == "Mixte"] <- "Autre"

Table_no_mixt <- Table_type[Table_type$Type_Clinique != "Mixte", ]

# Sélectionner la première colonne et les colonnes contenant "Streptococcus uberis"
selected_columns <- c(names(table_complete)[1], grep("Streptococcus.uberis", names(table_complete), value = TRUE))
new_table <- table_complete[, selected_columns]
write.csv2(new_table, file = "Matrix_Uberis_Only.csv")

# Transposer la table sans la première colonne
table_transpose <- t(new_table[,-1])
colnames(table_transpose) <- new_table$Pathway.Organism
table_complete <- new_table

# Fonction pour trouver la position du deuxième point
position_deuxieme_point <- function(x) {
  positions <- gregexpr("\\.", x)[[1]]  # Trouver toutes les positions des points
  if (length(positions) >= 2) {
    return(positions[2])  # Retourner la position du deuxième point
  } else {
    return(NA)  # Retourner NA si la chaîne n'a pas deux points
  }
}

# Appliquer la fonction à chaque élément du vecteur
result <- sapply(colnames(table_complete[,-c(1,2)]), position_deuxieme_point)
Type <- substr(colnames(table_complete[,-c(1,2)]), 1, result - 1)
names(Type) <- colnames(table_complete[,-c(1,2)])

# Harmoniser les noms de colonnes pour correspondre à "Table_type"
Table_type <- Table_type[,-c(3,4)]

# Filtrer les souches dans la table_complete pour ne conserver que celles présentes dans Table_type
valid_souches <- intersect(colnames(table_complete[,-c(1,2)]), Table_type$Souches)
table_complete <- table_complete[, c("Pathway.Organism", valid_souches)]

# Créer la matrice de données
matrix <- as.matrix(t(table_complete[,-1]))
colnames(matrix) <- table_complete$Pathway.Organism
var <- apply(matrix, 2, var)
matrix <- matrix[, var != 0]

# Associer les types cliniques aux souches de la matrice filtrée
Type_clinique <- Table_type$Type_Clinique2[match(rownames(matrix), Table_type$Souches)]
names(Type_clinique) <- rownames(matrix)

# Effectuer l'analyse PCA
pca_meta <- mixOmics::pca(matrix, scale = TRUE)
mixOmics::plotIndiv(pca_meta, group = Type_clinique, legend = TRUE, ind.names = FALSE)
plotVar(pca_meta, var.names = FALSE, legend = TRUE)

# PLS-DA
# Associer les types cliniques aux souches de la matrice filtrée
Type_clinique <- Table_type$Type_Clinique[match(rownames(matrix), Table_type$Souches)]
names(Type_clinique) <- rownames(matrix)
Result_PLSDA_metabolic_1 <- plsda(X = matrix, Y = Type_clinique, ncomp = 2)

plotIndiv(Result_PLSDA_metabolic_1, ind.names = FALSE,
          title = "PLSDA Subclinique, Clinique, Mixte",
          group = factor(Table_type$Type_Clinique, levels = c("Subclinique","Clinique","Mixte")),
          # fill = factor(Table_type$Type_Clinique, levels = c("Subclinique","Clinique","Mixte")),
          pch = factor(Table_type$Type_Clinique, levels = c("Subclinique","Clinique","Mixte")),
          legend = TRUE, comp = c(1, 2))
plotVar(Result_PLSDA_metabolic_1, var.names = T, cex = 3,
        title = "Correlation Circle Plot PLSDA Subclinique, Clinique, Mixte")
plotLoadings(Result_PLSDA_metabolic_1, contrib = 'max', method = 'mean', comp = 1)


Type_clinique <- Table_type$Type_Clinique2[match(rownames(matrix), Table_type$Souches)]
names(Type_clinique) <- rownames(matrix)
Result_PLSDA_metabolic_2 <- plsda(X = matrix, Y = Type_clinique, ncomp = 2)
plotIndiv(Result_PLSDA_metabolic_2,
          title = "PLSDA Subclinique, Autre",
          ind.names = FALSE,
          group = factor(Table_type$Type_Clinique2, levels = c("Subclinique","Autre")),
          legend = TRUE,
          comp = c(1, 2),
          pch = factor(Table_type$Type_Clinique2, levels = c("Subclinique","Autre")))
plotVar(Result_PLSDA_metabolic_2, var.names = T, cex = 3,
        title = "Correlation Circle Plot PLSDA Subclinique, Autre")
plotLoadings(Result_PLSDA_metabolic_2, contrib = 'max', method = 'mean', comp = 1)

Type_clinique <- Table_type$Type_Clinique3[match(rownames(matrix), Table_type$Souches)]
names(Type_clinique) <- rownames(matrix)
Result_PLSDA_metabolic_3 <- plsda(X = matrix, Y = Type_clinique, ncomp = 2)
plotIndiv(Result_PLSDA_metabolic_3,
          title = "PLSDA Clinique, Autre",
          ind.names = FALSE,
          group = factor(Table_type$Type_Clinique3, levels = c("Clinique","Autre")),
          legend = TRUE,
          comp = c(1, 2),
          pch = factor(Table_type$Type_Clinique3, levels = c("Clinique","Autre")))

plotVar(Result_PLSDA_metabolic_3, var.names = T, cex = 3,
        title = "Correlation Circle Plot PLSDA Clinique, Autre")
plotLoadings(Result_PLSDA_metabolic_3, contrib = 'max', method = 'mean', comp = 1)

# PLSDA non mixte
# Filtrer les souches pour ne conserver que celles présentes dans Table_no_mixt
valid_souches_no_mixt <- intersect(colnames(table_complete[,-c(1)]), Table_no_mixt$Souches)
matrix_no_mixt <- matrix[valid_souches_no_mixt, ]
var_no_mixt <- apply(matrix_no_mixt, 2, var)
matrix_no_mixt <- matrix_no_mixt[, var_no_mixt != 0]  # Remove columns with zero variance
Type_clinique_no_mixt <- Table_no_mixt$Type_Clinique[match(rownames(matrix_no_mixt), Table_no_mixt$Souches)]
names(Type_clinique_no_mixt) <- rownames(matrix_no_mixt)

Result_PLSDA_no_mixt <- plsda(X = matrix_no_mixt, Y = Type_clinique_no_mixt, ncomp = 2)

plotIndiv(Result_PLSDA_no_mixt, ind.names = FALSE, group = Type_clinique_no_mixt, legend = TRUE, title = "PLSDA non mixte")

plotVar(Result_PLSDA_no_mixt, var.names = T, cex = 3,
        title = "Correlation Circle Plot PLSDA Non Mixte")
plotLoadings(Result_PLSDA_no_mixt, contrib = 'max', method = 'mean', comp = 1)


# sPLS-DA
Result_sPLSDA <- splsda(X = matrix, Y = Type_clinique, ncomp = 2, keepX = c(50,50))
Selected_sPLSDA_metabolic <- selectVar(Result_sPLSDA, comp = 1)

Type_clinique <- Table_type$Type_Clinique[match(rownames(matrix), Table_type$Souches)]
names(Type_clinique) <- rownames(matrix)
Result_sPLSDA_metabolic_1 <- splsda(X = matrix, Y = Type_clinique, ncomp = 2, keepX = c(50,50))
plotIndiv(Result_sPLSDA_metabolic_1,
          title = "sPLSDA souches Streptococcus uberis",
          ind.names = FALSE,
          group = factor(Table_type$Type_Clinique, levels = c("Subclinique","Clinique","Mixte")),
          pch = factor(Table_type$Type_Clinique, levels = c("Subclinique","Clinique","Mixte")),
          legend = TRUE, 
          comp = c(1, 2),
          legend.text = list(cex = 20),
          size.xlabel = 17,
          size.ylabel = 17,
          size.legend = 17,
          legend.title = "")

plotVar(Result_sPLSDA_metabolic_1, var.names = T, pch = 16, cex = 3, col = 2,
        title = "Correlation Circle Plot sPLSDA souches Streptococcus uberis")
plotLoadings(Result_sPLSDA_metabolic_1, contrib = 'max', method = 'mean', comp = 1,
        title = "Contribution sPLSDA souches Streptococcus uberis",
        size.xlabel = 17,
        size.ylabel = 17,
        size.legend = 1.5,
        ndisplay = 20,
        legend = T,
        size.name = 0.8
        )

Type_clinique <- Table_type$Type_Clinique2[match(rownames(matrix), Table_type$Souches)]
names(Type_clinique) <- rownames(matrix)
Result_sPLSDA_metabolic_2 <- splsda(X = matrix, Y = Type_clinique, ncomp = 2, keepX = c(50,50))
plotIndiv(Result_sPLSDA_metabolic_2,
          title = "sPLSDA Subclinique, Autre",
          ind.names = FALSE,
          group = factor(Table_type$Type_Clinique2, levels = c("Subclinique","Autre")),
          # fill = factor(Table_type$Type_Clinique2, levels = c("Subclinique","Autre")),
          pch = factor(Table_type$Type_Clinique2, levels = c("Subclinique","Autre")),
          legend = TRUE, comp = c(1, 2))
plotVar(Result_sPLSDA_metabolic_2, var.names = T, pch = 16, cex = 3, col = 2,
        title = "Correlation Circle Plot sPLSDA Subclinique, Autre")
plotLoadings(Result_sPLSDA_metabolic_2, contrib = 'max', method = 'mean', comp = 1)

Type_clinique <- Table_type$Type_Clinique3[match(rownames(matrix), Table_type$Souches)]
names(Type_clinique) <- rownames(matrix)
Result_sPLSDA_metabolic_3 <- splsda(X = matrix, Y = Type_clinique, ncomp = 2, keepX = c(50,50))
plotIndiv(Result_sPLSDA_metabolic_3,
          title = "sPLSDA Clinique, Autre",
          ind.names = FALSE,
          group = factor(Table_type$Type_Clinique3, levels = c("Clinique","Autre")),
          # fill = factor(Table_type$Type_Clinique3, levels = c("Clinique","Autre")),
          pch = factor(Table_type$Type_Clinique3, levels = c("Clinique","Autre")),
          legend = TRUE, comp = c(1, 2))
plotVar(Result_sPLSDA_metabolic_3, var.names = T, pch = 16, cex = 3, col = 2,
        title = "Correlation Circle Plot sPLSDA Clinique, Autre")
plotLoadings(Result_sPLSDA_metabolic_3, contrib = 'max', method = 'mean', comp = 1)

# sPLSDA table non mixte
Result_sPLSDA_no_mixt <- splsda(X = matrix_no_mixt, Y = Type_clinique_no_mixt, ncomp = 2, keepX = c(50,50))

plotIndiv(Result_sPLSDA_no_mixt, title = "sPLSDA non mixte", ind.names = FALSE, legend = TRUE, comp = c(1,2))

plotVar(Result_sPLSDA_no_mixt, var.names = T, pch = 16, cex = 3, col = 2,
        title = "Correlation Circle Plot sPLSDA Non Mixte")
plotLoadings(Result_sPLSDA_no_mixt, contrib = 'max', method = 'mean', comp = 1)

Selected_sPLSDA_metabolic_no_mixt <- selectVar(Result_sPLSDA_no_mixt, comp = 1)

# Heatmap tous les types cliniques
annotation_col <- data.frame(Type_Clinique = factor(Table_type$Type_Clinique))
rownames(annotation_col) <- Table_type$Souches

# Définir les couleurs des annotations
annotation_colors <- list(
  Type_Clinique = c("Clinique" = "red", "Subclinique" = "darkgreen", "Mixte" = "yellow")
)

# Filtrer la matrice pour ne garder que les variables sélectionnées
selected_variables <- Selected_sPLSDA_metabolic$name
matrix_selected <- matrix[, selected_variables]

# Créer la heatmap avec les annotations
pheatmap(t(matrix_selected),
         annotation_col = annotation_col,
         clustering_method = "ward.D",
         cutree_cols = 1,
         show_colnames = FALSE,
         color = viridis_pal(option = "B", direction = 1)(256),
         annotation_colors = annotation_colors,
         fontsize = 14)

# Heatmap types cliniques non-mixte
annotation_col_no_mixt <- data.frame(Type_Clinique = factor(Table_no_mixt$Type_Clinique))
rownames(annotation_col_no_mixt) <- Table_no_mixt$Souches

# Définir les couleurs des annotations
annotation_colors_no_mixt <- list(
  Type_Clinique = c("Clinique" = "red", "Subclinique" = "darkgreen")
)

# Filtrer la matrice pour ne garder que les variables sélectionnées
selected_variables_no_mixt <- Selected_sPLSDA_metabolic_no_mixt$name
matrix_selected_no_mixt <- matrix_no_mixt[, selected_variables_no_mixt]

# Créer le heatmap avec les annotations pour le dataset non mixte
pheatmap(t(matrix_selected_no_mixt),
         annotation_col = annotation_col_no_mixt,
         clustering_method = "ward.D",
         cutree_cols = 1,
         show_colnames = FALSE,
         color = viridis_pal(option = "B", direction = 1)(256),
         annotation_colors = annotation_colors_no_mixt)


# Représentations voies spécifiques
table_merged <- merge(Table_type, table_transpose, by = "Souches")
col_graph <- grep("pyridoxal 5_-phosphate salvage pathway", colnames(table_merged))
col_graph <- grep("lactose degradation III", colnames(table_merged))
col_graph <- grep("acyl-ACP thioesterase", colnames(table_merged))

col_graph <- grep("NAD salvage pathway I$", colnames(table_merged))

ggplot(table_merged, aes(x = Type_Clinique, y = as.numeric(table_merged[,col_graph]))) +
  geom_jitter(width =0.2,
              height = 0) +
  ylim(0,1)


ggplot(table_merged, aes(x = as.factor(Type_Clinique), y = as.numeric(table_merged[,col_graph]))) +
  geom_violin(trim=FALSE)

#stat
matrix <- as.matrix(t(table_complete[,-c(1,2)]))
colnames(matrix) <- table_complete$Pathway.Organism
data_frame_pres_abs <- as.data.frame(matrix)
data_frame_pres_abs[data_frame_pres_abs == 0] <- "absent"
data_frame_pres_abs[data_frame_pres_abs == 1] <- "complete"
data_frame_pres_abs[data_frame_pres_abs < 1 & data_frame_pres_abs >0] <- "incomplete"

data_frame_pres_abs <- cbind(Type, Mammites, Human, data_frame_pres_abs)

chi_square_table <- data.frame(Pathway = colnames(data_frame_pres_abs))
for (pathway in 2:ncol(data_frame_pres_abs)) {
  # pathway <- 2
  colnames(data_frame_pres_abs)[pathway]
  chi <- chisq.test(table(data_frame_pres_abs$Type, data_frame_pres_abs[, pathway]), simulate.p.value = TRUE, B = 2000)
  chi_square_table[pathway,2] <- chi$p.value
}

chi_square_table$adj_pvalue <- p.adjust(chi_square_table$V2, method = "fdr")

summary(chi_square_table$adj_pvalue < 0.05)  


#mammites
chi_square_table_mastitis <- data.frame(Pathway = colnames(data_frame_pres_abs))
for (pathway in 2:ncol(data_frame_pres_abs)) {
  # pathway <- 2
  colnames(data_frame_pres_abs)[pathway]
  chi <- chisq.test(table(data_frame_pres_abs$Mammites, data_frame_pres_abs[, pathway]), simulate.p.value = TRUE, B = 2000)
  chi_square_table_mastitis[pathway,2] <- chi$p.value
}

chi_square_table_mastitis$adj_pvalue <- p.adjust(chi_square_table_mastitis$V2, method = "fdr")

summary(chi_square_table_mastitis$adj_pvalue < 0.01)


#human pathogen
chi_square_table_human <- data.frame(Pathway = colnames(data_frame_pres_abs))
for (pathway in 2:ncol(data_frame_pres_abs)) {
  # pathway <- 2
  colnames(data_frame_pres_abs)[pathway]
  chi <- chisq.test(table(data_frame_pres_abs$Human, data_frame_pres_abs[, pathway]), simulate.p.value = TRUE, B = 2000)
  chi_square_table_human[pathway,2] <- chi$p.value
}

chi_square_table_human$adj_pvalue <- p.adjust(chi_square_table_human$V2, method = "fdr")

summary(chi_square_table_human$adj_pvalue < 0.01)  

table(data_frame_pres_abs$Type)


#table_8_uberis
table_sans_uberis <- data_frame_pres_abs[-grep("s.uberis", data_frame_pres_abs$Type),]
rownum_uberis <- grep("s.uberis", data_frame_pres_abs$Type)
chi_square_table_8_uberis <- data.frame(Pathway = colnames(table_sans_uberis))

for (i in 1:100) {
  random <- sample(rownum_uberis)[1:8]
  table_8_uberis <- rbind(table_sans_uberis, data_frame_pres_abs[random,])
  table(table_8_uberis$Type)
  for (pathway in 2:ncol(table_8_uberis)) {
   # pathway <- 2
   colnames(table_8_uberis)[pathway]
    chi <- chisq.test(table(table_8_uberis$Type, table_8_uberis[, pathway]))
    chi_square_table_8_uberis[pathway, i+1] <- chi$p.value
  }
  chi_square_table_8_uberis[,i+1] <- p.adjust(chi_square_table_8_uberis[,i+1], method = "BH")
}

chi_square_table_8_uberis_mean <- apply(chi_square_table_8_uberis[,-1],1, median)
chi_square_table_8_uberis_mean <- data.frame(Pathway = colnames(table_8_uberis), pvalue = chi_square_table_8_uberis_mean)
chi_square_table_8_uberis_mean <- chi_square_table_8_uberis_mean[-c(1:3),]
summary(chi_square_table_8_uberis_mean$pvalue < 0.05)  