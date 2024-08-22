#!/bin/bash
# Réalise les roary sur l'ensemble des gff du dossier input

conda activate roary-3.13.0 && roary -e --mafft -f /Chemin/dossier/output -p 8 /Chemin/dossier.input/*.gff -i 80 -cd 99 -g 100000 && conda deactivate