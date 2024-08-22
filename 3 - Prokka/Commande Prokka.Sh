#!/bin/bash

# Définir le dossier contenant les fichiers .fasta
dossier_fasta=/save_home/mxleclercq/Raw_data_reads/

# Définir le dossier de sortie pour les résultats de Prokka
dossier_output=/save_home/mxleclercq/Raw_data_reads/Resultat_Prokka/

# Activer l'environnement conda pour Prokka
conda activate prokka-1.14.6

# Boucle sur les fichiers .fasta dans le dossier
for fichier in "$dossier_fasta"/*.fasta; do
  # Extraire le nom du fichier sans l'extension
  prefix=$(basename "$fichier" .fasta)
  
  # Exécuter la commande Prokka
  prokka --outdir "$dossier_output/$prefix" --prefix "$prefix" "$fichier" --centre X --compliant
done

# Désactiver l'environnement conda
conda deactivate