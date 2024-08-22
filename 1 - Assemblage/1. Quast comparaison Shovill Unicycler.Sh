#!/bin/bash
conda activate quast-5.0.2 &&

# Permet de réaliser l'analyse quast afin de vérifier les données d'assemblages des différentes souches de S. uberis

# Définir les répertoires d'entrée et de sortie
input_base="/save_projet/mastisub/resultats_shovill"
output_base="/save_projet/mastisub/Quast"

# Trouver tous les fichiers fasta dans le répertoire d'entrée et ses sous-répertoires
find "$input_base" -type f -name "*.fasta" | while read file; do
    # Get the relative path of the file
    relative_path="${file#$input_base/}"
    
    # Extraire le nom du répertoire contenant le fichier
    relative_dir="$(dirname "$relative_path")"
    
    # Créer le répertoire de sortie correspondant
    mkdir -p "$output_base/$relative_dir"
    
    # Commande Quast et sauvegarde des résultats dans le répertoire de sortie correspondant
    quast.py "$file" -o "$output_base/$relative_dir"

done

&& conda deactivate