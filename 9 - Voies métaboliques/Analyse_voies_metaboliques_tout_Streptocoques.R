library(mixOmics)
library(WallomicsData)
library(cowplot)
library(pheatmap)
library(viridis)

setwd("/home/maximeleclercq/Bureau/Clean/9 - Voies Métaboliques")
table_complete <- read.csv2("Rapport_Voies_Metaboliques_V2.csv", sep = "\t", dec = ".")

table_transpose <- t(table_complete[,-2])

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

Species <-substr(colnames(table_complete[,-c(1,2)]), 1, result -1 )
names(Species) <- colnames(table_complete[,-c(1,2)])
{
        Mammites <- Species
        Mammites[grep("agalactiae", Mammites)] <- "mastitis"
        Mammites[grep("uberis", Mammites)] <- "mastitis"
        Mammites[grep("dysga", Mammites)] <- "mastitis"
        Mammites[grep("paraub", Mammites)] <- "mastitis"
        Mammites[grep("pyoge", Mammites)] <- "non-mastitis"
        Mammites[grep("canis", Mammites)] <- "non-mastitis"
        
        Human <- Species
        Human[grep("agalactiae", Human)] <- "pathogen"
        Human[grep("uberis", Human)] <- "non-pathogen"
        Human[grep("dysga", Human)] <- "pathogen"
        Human[grep("paraub", Human)] <- "non-pathogen"
        Human[grep("pyoge", Human)] <- "pathogen"
        Human[grep("canis", Human)] <- "pathogen"
}

matrix <- as.matrix(t(table_complete[,-c(1,2)]))
colnames(matrix) <- table_complete$Pathway.Organism
var <- apply(matrix, 2, var)

matrix <- matrix[,var!=0]
pca_meta <- pca(matrix, scale = T)

# Represent the individuals according to their ecotype
plotIndiv(pca_meta,
          group = Species,
          legend = TRUE,
          ind.names = F
          )

plotVar(pca_meta,
        var.names = F,
        legend = TRUE
        )

#PLS-DA
Result_PLSDA_metabolic <- 
  plsda(X = matrix, Y = Species,
        ncomp = 2)

plotIndiv(Result_PLSDA_metabolic,
          ind.names = F,
          legend = T,
          comp = c(1,2))

plotIndiv(Result_PLSDA_metabolic,
          ind.names = FALSE,
          group = Species, 
          legend = TRUE,
          title = "PLSDA tout Strepto")

plotVar(Result_PLSDA_metabolic,
        var.names = T, pch = 16, cex = 3, col = 2,
        title = "Correlation Circle Plot PLSDA tout Strepto")

plotLoadings(Result_PLSDA_metabolic,
             contrib = 'max', method = 'mean', comp = 1,
             title = "Contribution PLSDA tout Strepto")



# sPLS-DA
Result_sPLSDA_metabolic <-
  splsda(X = matrix,
         Y = Species,
         ncomp = 2,
         keepX = c(50,50))



plotIndiv(Result_sPLSDA_metabolic,
          title = "sPLSDA tout Strepto",
          ind.names = F,
          legend = T)

plotVar(Result_sPLSDA_metabolic,
        var.names = T, pch = 16, cex = 3, col = 2,
        title = "Correlation Circle Plot sPLSDA tout Strepto")

plotLoadings(Result_sPLSDA_metabolic,
             contrib = 'max', method = 'mean', comp = 1,
             title = "Contribution sPLSDA tout Strepto")

Selected_sPLSDA_metabolic <-
  selectVar(Result_sPLSDA_metabolic, comp = 1)

# Heatmap
data_annotations <- read.csv2("Mammite-Pathogènes_V2.csv", sep = "\t")

annotation_col <- data.frame(
  Mammite = factor(data_annotations$Mammite),
  Pathogène_Homme = factor(data_annotations$Homme),
  Espèce = factor(data_annotations$Espèce)
)
rownames(annotation_col) <- data_annotations$Souches

annotation_colors <- list(
  Mammite = c("Oui" = "blue", "Non" = "lightgrey"),
  Pathogène_Homme = c("Oui" = "purple", "Non" = "#B80F0A"),
  Espèce = c("Streptococcus agalactiae" = "green",
             "Streptococcus canis" = "#fbfeed",
             "Streptococcus dysgalactiae" = "grey",
             "Streptococcus parauberis" = "red",
             "Streptococcus pyogenes" = "turquoise",
             "Streptococcus uberis" = "darkorange")
)

pheatmap(t(matrix[, is.element(colnames(matrix), Selected_sPLSDA_metabolic$name)]),
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         clustering_method = "ward.D2",
         cutree_cols = 6,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         color = viridis_pal(option = "B", direction = 1)(256),
         show_colnames = F,
         fontsize = 14,
         main = "Heatmap des voies métaboliques des différentes espèces de Streptocoques"
         )

#stat
matrix <- as.matrix(t(table_complete[,-c(1,2)]))
colnames(matrix) <- table_complete$Pathway.Organism
data_frame_pres_abs <- as.data.frame(matrix)
data_frame_pres_abs[data_frame_pres_abs == 0] <- "absent"
data_frame_pres_abs[data_frame_pres_abs == 1] <- "complete"
data_frame_pres_abs[data_frame_pres_abs < 1 & data_frame_pres_abs >0] <- "incomplete"

data_frame_pres_abs <- cbind(Species, Mammites, Human, data_frame_pres_abs)

chi_square_table <- data.frame(Pathway = colnames(data_frame_pres_abs))
for (pathway in 2:ncol(data_frame_pres_abs)) {
  # pathway <- 2
  colnames(data_frame_pres_abs)[pathway]
  chi <- chisq.test(table(data_frame_pres_abs$Species, data_frame_pres_abs[, pathway]), simulate.p.value = TRUE, B = 2000)
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

table(data_frame_pres_abs$Species)


#table_8_uberis
table_sans_uberis <- data_frame_pres_abs[-grep("s.uberis", data_frame_pres_abs$Species),]
rownum_uberis <- grep("s.uberis", data_frame_pres_abs$Species)
chi_square_table_8_uberis <- data.frame(Pathway = colnames(table_sans_uberis))

for (i in 1:100) {
  random <- sample(rownum_uberis)[1:15]
  table_8_uberis <- rbind(table_sans_uberis, data_frame_pres_abs[random,])
  table(table_8_uberis$Species)
  for (pathway in 2:ncol(table_8_uberis)) {
   # pathway <- 2
   colnames(table_8_uberis)[pathway]
    chi <- chisq.test(table(table_8_uberis$Species, table_8_uberis[, pathway]))
    chi_square_table_8_uberis[pathway, i+1] <- chi$p.value
  }
  chi_square_table_8_uberis[,i+1] <- p.adjust(chi_square_table_8_uberis[,i+1], method = "bonferroni")
}

chi_square_table_8_uberis_mean <- apply(chi_square_table_8_uberis[,-1],1, median)
chi_square_table_8_uberis_mean <- data.frame(Pathway = colnames(table_8_uberis), pvalue = chi_square_table_8_uberis_mean)
chi_square_table_8_uberis_mean <- chi_square_table_8_uberis_mean[-c(1:3),]
summary(chi_square_table_8_uberis_mean$pvalue < 0.05)  



