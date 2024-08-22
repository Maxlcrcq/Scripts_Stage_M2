library(seqinr)
library(writexl)
library(readxl)

# Charger le fichier Excel
data <- readxl::read_excel("C:/Users/lecle/Desktop/table_VF_streptococcus_vezina.xlsx")

for (i in 1:nrow(data)) {
  id <- data[i, "ID"]
  sequence <- unlist(strsplit(as.character(data[i, "AA_sequence"]), ""))
  pdf_name <- paste0(id, ".pdf")
  pdf(pdf_name, 12, 12)
  plot_title <- paste("Dotplot", id)
  seqinr::dotPlot(sequence, sequence, wsize = 20, nmatch = 5, main = plot_title)
  dev.off()
}