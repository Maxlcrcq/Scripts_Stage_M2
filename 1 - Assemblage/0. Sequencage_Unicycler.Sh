#!/bin/bash

# Permet l'assemblage des séquençages de S. uberis en utilisant Unicycler

# Répertoire de travail
WORKDIR=/chemin/vers/les/fichiers

# Répertoire de sortie
OUTDIR=/chemin/vers/le/dossier/de/sortie

# Boucle sur les fichiers _1.fq.gz
for file in ${WORKDIR}/*_1.fq.gz; do
  # Extraire le nom de l'échantillon
  sample=${file%_1.fq.gz}

  # Trouver le fichier _2.fq.gz correspondant
  R2=${sample}_2.fq.gz

  # Vérifier si le fichier _2.fq.gz existe
  if [ -f "${WORKDIR}/${R2}" ]; then
    # Exécuter la commande Unicycler
    conda activate unicycler-0.5.0 && \
    unicycler -1 ${file} -2 ${WORKDIR}/${R2} -o ${OUTDIR}/${sample} && \
    conda deactivate
  else
    echo "Erreur : fichier ${R2} introuvable"
  fi
done