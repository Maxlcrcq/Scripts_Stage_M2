#!/bin/bash

# Comparaison quast des assembleurs Shovill et Unicycler

for file in /save_projet/mastisub/resultats_shovill/*/*.fasta; do
  base_name=$(basename "$file" .fasta)
  dir_name=$(dirname "$file")
  quast_out_dir=/save_projet/mastisub/resultats_quast/${dir_name##*/}
  mkdir -p "$quast_out_dir"
  conda activate quast-5.0.2 && quast.py "$file" "/save_projet/mastisub/resultats_assemblages/${dir_name##*/}/${base_name}.fasta" -o "$quast_out_dir" --circos && conda deactivate
done