# Permet de lancer l'alignement des séquences protéiques de S. uberis, obtenues à l'aide du script R, avec MUSCLE

conda activate muscle-3.8.1551

for ((i=1; i<=70; i++))
do
    muscle -in alignement_strepto_uberis_w_s_porcinus_$i.faa -out Aligned_S_uberis_$i.fasta
done

conda deactivate
