# Calcul des OTU 

Vous trouverez la description complète du TP [ici](https://docs.google.com/document/d/1qWNqPZ9Ecd-yZ5Hpl6n2zd7ZGtHPjf3yaW1ulKRdWnk/edit?usp=sharing).

## Résultats

Pytest : \
![](/results/pytest.png)

Exécution : \
![](/results/execution.png)

Nous identifions de l'ordre de 100 OTUs dans l'échantillon proposé, contenant 8 espèces bactériennes en quantité d'ADN équivalente et 2 espèces de levure 6 fois moins représentées.

vsearch : \
![](/results/vsearch.png)


L'échantillon à l'origine des *reads* séquencés contenant seulement 10 espèces, nous nous sommes intéressées à réduire le pourcentage d'identité pour le regroupement des OTUS en cours d'algorithme.  
La version modifiée du script (ajout d'un argument optionnel en plus) est consultable sur la branche `test-treshold_bene`, les fichiers de résultats sont présents dans la branche `master`. Exécution du programme en 8 minutes 26 secondes environ.

vsearch identity treshold 95 : \
![](/results/vsearch_treshold95.PNG)

Le changement du pourcentage d'identité (de 97 à 95%, nomnbre d'OTUs du même ordre que la référence dans l'expérience de Kunin *et al.*, 2010) réduit le nombre d'OTUs identifiées à seulement 45 au lieu de 114.  
Toutefois, une séquence demeure non reconnue dans le lot d'OTUs identifiées. Il peut s'agir d'une séquence contaminante lors de la préparation du mélange (autre espèce) ou à celle de la librairie de séquençage (préparateur).   


## Introduction

L’objectif de ce TP sera de calculer les OTU obtenues à partir d’un séquençage “mock”. Nous n’avons amplifié que les bactéries (et non les champignons). 8 espèces sont ainsi attendues.

Vous devrez développer un programme effectuant une dé-duplication en séquence complète (“dereplication full length”), une recherche des séquences chimériques et un regroupement basé sur un algorithme glouton (“Abundance Greedy Clustering”).  


## Installation des dépendances

Vous utiliserez les librairies nwalign3, pytest et pylint de Python:
```
pip3 install --user nwalign3 pytest pylint pytest-cov
```

## Utilisation

Vous devrez développer un programme python3 effectuant une dé-duplication en séquence complète (“dereplication full length”), une recherche des séquences chimériques et un regroupement basé sur un algorithme glouton (“Abundance Greedy Clustering”). Il prendra pour arguments:

 -i, -amplicon_file fichier contenant des séquences au format FASTA
 -s, -minseqlen Longueur minimum des séquences (optionnel - valeur par défaut 400)
 -m, -mincount Comptage minimum des séquences (optionnel - valeur par défaut 10)
 -c, -chunk_size Taille des partitions de séquence (optionnel - valeur par défaut 100)
 -k, -kmer_size Longueur des “kmer” (optionnel - valeur par défaut 8)
 -o, -output_file fichier de sortie avec les OTU au format FASTA

 ## Tests

Vous testerez vos fonctions à l’aide de la commande pytest --cov=agc à exécuter dans le dossier agc-tp/. En raison de cette contrainte, les noms des fonctions ne seront pas libre. Il sera donc impératif de respecter le nom des fonctions “imposées”, de même que leur caractéristique et paramètres. 
Vous vérifierez également la qualité syntaxique de votre programme en exécutant la commande: pylint agc.py

## Contact

En cas de questions, vous pouvez me contacter par email: amine.ghozlane[at]pasteur.fr

