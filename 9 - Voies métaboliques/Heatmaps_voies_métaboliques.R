#Chargement des packages
library(pheatmap)
library(dplyr)
library(viridis)

setwd("/home/maximeleclercq/Bureau/Clean/9 - Voies Métaboliques/V2")
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

