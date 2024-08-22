#!/bin/bash

# Permet l'assemblage des séquençages de S. uberis en utilisant Shovill

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
    # Exécuter la commande Shovill
    conda activate shovill-1.1.0 && \
    shovill --R1 ${file} --R2 ${WORKDIR}/${R2} --outdir ${OUTDIR}/${sample} && \
    conda deactivate
  else
    echo "Erreur : fichier ${R2} introuvable"
  fi
done