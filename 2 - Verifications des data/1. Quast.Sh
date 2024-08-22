#!/bin/bash
# Permet de réaliser l'analyse quast afin de vérifier les données d'assemblages des différentes souches de S. uberis
conda activate quast-5.0.2 &&

# Répertoires d'entrée et de sortie
dossier_input="/save_projet/mastisub/resultats_shovill"
dossier_output="/save_projet/mastisub/Quast"

find "$dossier_input" -type f -name "*.fasta" | while read file; do
    relative_path="${file#$dossier_input/}"
    relative_dir="$(dirname "$relative_path")"
    mkdir -p "$dossier_output/$relative_dir"
    quast.py "$file" -o "$dossier_output/$relative_dir"

done

&& conda deactivate