# README

## Objectifs

Définir une configuration des paramètres d'entrée du modèle PROSAIL:
1. D'abord une configuration optimale par site.
2. Voir s'il est possible d'en trouver une généralisable à tout site.

## Méthodes
- La configuration ATBD (https://hal.inrae.fr/hal-03584016) est la configuration de référence. Pour chaque configuration testée, nous partons de la configuration ATBD (co-distributions désactivées) et modifiont une ou plusieurs variables. 
- Les paramètres à modifier sont définis.
- Echantillonner un certain nombre de points pour ne pas avoir à inverser sur tous les pixels (gain de temps).
- Pour chaque configuration, calculer la corrélation de Pearson entre le LAI Sentinel-2 obtenu par inversion de PROSAIL avec le LAI LiDAR.
- Récupérer les distributions des paramètres associés à la meilleure corrélation.

### Paramètres à modifier
- LAI (autre possibilité: mettre la distribution du LAI LiDAR)
- LIDFa
- EWT
- LMA
- CHL
- ...
Chacun peut avoir une distribution uniforme ou gaussienne.

### Echantillonnage
- Echantillonnage aléatoire ou stratifié entre 2,000 et 10,000 points par site.

### Parallélisation
Le processus pourrait être parallélisé avec une librairie comme future.