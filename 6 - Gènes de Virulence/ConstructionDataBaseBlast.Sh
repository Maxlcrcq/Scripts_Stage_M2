#!/bin/bash

# Permet de créer une base de données pour chaque souoche .faa

for file in *.faa; do
  base_name=$(basename "$file" .faa)
  makeblastdb -in "$file" -dbtype prot -parse_seqids -out "$base_name"
done