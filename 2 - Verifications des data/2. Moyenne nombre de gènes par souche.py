#Permet de calculer la moyenne du nombre de gènes présents dans les différentes souches de Streptococcus uberis

import os

def moyenne_genes(dossier):
    total_genes = 0
    total_fichiers = 0
    for fichier in os.listdir(dossier):
        if fichier.endswith(".txt"):
            chemin_fichier = os.path.join(dossier, fichier)
            with open(chemin_fichier, 'r') as f:
                lignes = f.readlines()
                for ligne in lignes:
                    if ligne.startswith("gene:"):
                        nb_genes = int(ligne.split(":")[1])
                        total_genes += nb_genes
                        total_fichiers += 1
                        break
    if total_fichiers > 0
        moyenne = total_genes / total_fichiers
        return moyenne
    else:
        return "Erreur : aucun fichier valide trouvé"

# Utilisation du script
dossier = "chemin/vers/le/dossier"
moyenne = moyenne_genes(dossier)
print("La moyenne de gènes présents est :", moyenne)
