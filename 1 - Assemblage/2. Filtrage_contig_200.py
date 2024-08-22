# Permet de filtrer les séquences fasta de moins de 200 nucléotides pour l'envoi sur la plateforme du genoscope
import os

from Bio import SeqIO

dossier_entrée = 'Chemin/vers/le/dossier/entrée'
dossier_sortie = 'Chemin/vers/le/dossier/sortie'
if not os.path.exists(dossier_sortie):
    os.makedirs(dossier_sortie)

for file_name in os.listdir(dossier_entrée):
    if file_name.endswith('.fasta') or file_name.endswith('.fna'):
        file_path = os.path.join(dossier_entrée, file_name)
        output_file_path = os.path.join(dossier_sortie, file_name)
        with open(file_path, 'r') as input_file, open(output_file_path, 'w') as output_file:
            for record in SeqIO.parse(input_file, 'fasta'):
                if len(record.seq) >= 200:
                    SeqIO.write(record, output_file, 'fasta')
