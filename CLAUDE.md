# CLAUDE.md — Projet LiDAR × Sentinel-2 LAI (Corroyez et al., RSE révision)

## Contexte scientifique

Ce projet compare des estimations de Leaf Area Index (LAI) dérivées de deux
sources : airborne LiDAR (ALS, méthode gap fraction avec loi de Beer-Lambert,
Bouvier et al. 2015) et Sentinel-2 (inversion hybride PROSAIL via le package R
`prosail` de Féret & de Boissieu). Trois forêts tempérées décidues françaises :
Aigoual (hêtraie montagnarde, Cévennes), Blois (chênaie sessile de plaine),
Mormal (chênaie-hêtraie de plaine, Nord). Acquisitions ALS et S2 en été 2021,
résolution de travail 10 m (grille S2).

Le papier (RSE-D-25-04417) identifie trois facteurs expliquant les écarts
LAI_ALS vs LAI_S2 : (1) profondeur effective de canopée vue par S2 (d_opt),
(2) paramétrisation du LUT PROSAIL pour forêts vs croplands ATBD,
(3) hétérogénéité horizontale de la canopée (CHM_std). Le repositionnement
stratégique de la révision passe d'une logique "improve accuracy" à
"improve inter-sensor consistency" car il n'y a pas de validation terrain.

Le projet est actuellement en phase de révision majeure suite aux reviews.
Les analyses additionnelles à produire concernent principalement :
- Sensibilité de LAI_ALS au coefficient d'extinction k (test {0.4, 0.5, 0.6})
- Sensibilité aux seuils h_min (2/3/5 m) et fCover (80/90/95 %)
- Comparaison CHM_std vs DSM_std comme métrique d'hétérogénéité
- Sélection multi-critères de d_opt (Pareto au lieu de Pearson seul)
- Analyse "Mormal phase 1 only" pour quantifier l'impact du décalage temporel

## Structure du dépôt

Le dossier racine est `~/Documents/NC_Full/`. Il est historiquement en
désordre. Claude Code ne doit PAS explorer ou modifier librement tout le dépôt.

### Dossiers en lecture seule — NE JAMAIS MODIFIER NI SUPPRIMER
- `01_DATA/` : données brutes (LiDAR, S2, BDForêt, GEDI). En lecture seule.
  Ne jamais écrire dedans. Ne jamais lancer de téléchargement qui y écrirait.
- `03_RESULTS/` : sorties historiques des pipelines. Lecture seule par défaut.
  Si une analyse doit réécrire, me demander confirmation explicite avant.
- `04_FIGURES/` : figures historiques. Lecture seule.
- `paper1_Figures/` : figures finales du manuscrit soumis. Lecture seule
  absolue, c'est la version de référence.
- `PROSAIL-Optimization/` : sous-projet dédié à l'optimisation PROSAIL, a sa
  propre logique. Ne pas toucher sans instruction explicite.
- `02_CODES/archive/`, `02_CODES/other/`, `02_CODES/other_finalize_p1/` :
  scripts legacy, souvent obsolètes ou redondants. Lire pour comprendre
  l'historique mais ne jamais exécuter ni modifier.

### Dossiers de travail pour la révision — zone de travail principale
- `revision/` : à créer. C'est là que va tout le nouveau code refactorisé et
  les analyses de révision. Toute nouvelle fonction, tout nouveau script,
  toute nouvelle figure va là. Structure attendue :
  - `revision/R/` : fonctions (une fonction publique par fichier quand c'est
    raisonnable, groupées par thème sinon)
  - `revision/scripts/` : scripts d'orchestration numérotés
    (`01_prepare_data.R`, `02_lai_als_sensitivity.R`, etc.)
  - `revision/output/figures/` : nouvelles figures
  - `revision/output/tables/` : nouvelles tables
  - `revision/output/intermediate/` : résultats intermédiaires (rds, parquet)
  - `revision/tests/` : tests unitaires basiques
  - `revision/vignettes/` : documentation exécutable
  - `revision/README.md` : point d'entrée pour les encadrants

### Scripts existants à refactoriser (pipeline du papier soumis)
Les scripts qui produisent les figures du papier soumis et qui sont la cible
prioritaire du refactoring se trouvent principalement dans :
- `02_CODES/LiDAR/` : calcul des métriques LiDAR, PAD profiles, masks
- `02_CODES/Sentinel_2/` : téléchargement, preprocessing, inversion PROSAIL
- `02_CODES/libraries/` : fonctions utilitaires (functions_*.R)
- `PROSAIL-Optimization/02_CODES/Main_*.R` : scripts d'analyse PROSAIL

Avant de refactoriser, toujours me demander quel script correspond à quelle
figure du papier. Ne jamais supposer.

## Conventions de code

### Langage et style
- R principalement. Python uniquement si un script Python existe déjà et doit
  être modifié (ex : téléchargements GEDI).
- Style : préférer du R fonctionnel simple et lisible. Pas d'OOP S4/R6 sauf
  si déjà présent. Pas de métaprogrammation inutile.
- Utiliser `data.table` et base R de préférence. `dplyr`/`tidyr` acceptés si
  le script existant l'utilise déjà, ne pas mélanger les deux dans un même
  fichier.
- Toujours utiliser `here::here()` pour les chemins, jamais de chemins
  absolus hardcodés, jamais de `setwd()`.
- Encodage UTF-8, fins de ligne LF.
- Indentation : 2 espaces, pas de tabs.
- Longueur de ligne : 90 caractères max.

### Packages autorisés
- Géospatial : `terra`, `sf`, `lidR`, `stars`
- Statistiques : `mgcv`, `stats`, base R
- PROSAIL : `prosail`, `prospect` (versions CRAN de Féret)
- Manipulation : `data.table`, optionnellement `dplyr` si cohérence locale
- Graphiques : `ggplot2`, `patchwork`, `scales`
- Utilitaires : `here`, `fs`, `cli`, `glue`, `yaml`
- Parallélisation : `future`, `future.apply`

Ne PAS introduire de nouveaux packages sans me demander. En particulier,
éviter `tidymodels`, `targets`, `drake`, `renv` sauf instruction explicite.

### Documentation des fonctions
Toute nouvelle fonction dans `revision/R/` doit avoir un bloc roxygen2
complet :
- `@title` en une phrase
- `@description` en 2-4 phrases expliquant le rôle scientifique
- `@param` pour chaque argument avec type et unité (ex : "hauteur en mètres")
- `@return` avec type et structure
- `@examples` avec un exemple minimal reproductible si possible
- `@references` quand une fonction implémente une méthode publiée
  (ex : Bouvier et al. 2015 pour LAD via gap fraction)

### Nommage
- Fonctions : `snake_case`, verbes (`compute_lad_profile`, `estimate_lai_als`)
- Variables : `snake_case`, noms explicites avec unités quand pertinent
  (`h_min_m`, `lai_als_dopt`)
- Fichiers R : `snake_case.R`
- Scripts numérotés : `NN_short_description.R`
- Pas d'abréviations cryptiques. `lidar` pas `ldr`, `sentinel` pas `s2` dans
  les noms de fonctions (ok dans les noms de variables courtes).

### Gestion des chemins et I/O
- Tous les chemins passent par `here::here()` relatif au root du projet.
- Les lectures depuis `01_DATA/` sont OK. Les écritures vers `01_DATA/`
  sont INTERDITES.
- Toutes les sorties vont dans `revision/output/`.
- Fichiers intermédiaires : utiliser `qs` ou `rds` pour du R, `parquet`
  pour du tabulaire volumineux.
- Les fichiers intermédiaires lourds (>100 MB) doivent être dans
  `revision/output/intermediate/` et ajoutés au `.gitignore`.

## Règles de sécurité et de prudence

### Interdictions absolues
- Ne jamais exécuter `rm`, `file.remove`, `unlink` sur un fichier ou dossier
  sans me le demander explicitement et attendre ma confirmation.
- Ne jamais modifier ou écrire dans `01_DATA/`, `03_RESULTS/`, `04_FIGURES/`,
  `paper1_Figures/`.
- Ne jamais lancer de téléchargement S2, GEDI, Copernicus sans confirmation.
  Les credentials Copernicus sont dans `01_DATA/copernicus-credentials*.yml`
  et ne doivent jamais être lus, copiés, ou affichés.
- Ne jamais faire de `git push` vers un remote. Les commits locaux sont OK
  après ma validation.
- Ne jamais modifier `CLAUDE.md` lui-même sans me demander.

### Calculs lourds — règle du dry run obligatoire
Pour toute analyse qui va tourner sur :
- plus de 10 000 pixels S2, ou
- plus de 10 configurations PROSAIL, ou
- un traitement LiDAR sur un point cloud complet (pas un tuile),

commencer impérativement par un dry run sur un sous-échantillon (100 pixels,
2 configs, une tuile). Me montrer le résultat, les temps de calcul estimés,
et attendre ma validation avant le full run.

### Git
- Travailler sur des branches dédiées par thème :
  `revision/refactor-lidar`, `revision/k-sensitivity`,
  `revision/dopt-pareto`, etc.
- Commits atomiques, messages en anglais, format conventionnel :
  `feat(lai_als): add k sensitivity test`, `refactor(prosail): extract LUT`.
- Ne jamais merger vers `main` sans mon accord explicite.
- Avant toute modification d'un script existant dans `02_CODES/`, vérifier
  qu'on est sur une branche de révision et pas sur `main`.

## Priorités de travail (ordre suggéré)

Cette liste reflète l'ordre dans lequel on va probablement attaquer les
tâches. Ne pas sauter d'étapes sans raison.

1. **Inventaire et cartographie** : lire les scripts pertinents dans
   `02_CODES/LiDAR/`, `02_CODES/Sentinel_2/`, `02_CODES/libraries/` et
   produire un rapport (`revision/INVENTORY.md`) qui mappe chaque script
   aux figures du papier qu'il produit. Ne rien modifier à ce stade.

2. **Mise en place de l'infrastructure** : créer l'arborescence
   `revision/`, initialiser un `DESCRIPTION` minimal si utile, configurer
   `.gitignore`, créer le `README.md` pour les encadrants.

3. **Refactoring module par module**, dans cet ordre :
   a. Fonctions LiDAR (normalisation, CHM, DTM, LAD profiles, LAI_ALS)
   b. Fonctions Sentinel-2 (preprocessing, masks)
   c. Fonctions PROSAIL (LUT generation, inversion, metrics)
   d. Fonctions d'analyse (d_opt, heterogeneity, Pareto)
   e. Scripts d'orchestration qui chaînent le tout

   Un module à la fois. Review humain entre chaque. Commit git entre chaque.

4. **Analyses additionnelles pour la révision**, dans cet ordre :
   a. Sensibilité à k (rapide, post-processing des LAD existantes)
   b. Sensibilité à h_min et fCover (un site d'abord)
   c. CHM_std vs DSM_std (un site d'abord, puis extension)
   d. Sélection multi-critères de d_opt (Pareto)
   e. Analyse "Mormal phase 1 only"

5. **Documentation finale** : vignette Quarto reproductible sur un
   mini-dataset, README encadrants finalisé, éventuellement un DOI Zenodo.

## Spécificités scientifiques à respecter

- **LAD computation** : la formule est celle de Bouvier et al. 2015
  (Eq. 1 du papier), avec k = 0.5 par défaut. Toute modification de k doit
  être paramétrable, pas hardcodée.
- **Normalisation CHM vs DTM** : les deux approches doivent rester
  disponibles. Le papier utilise CHM-based en principal, DTM-based en
  comparaison. Ne pas supprimer l'une au profit de l'autre.
- **Seuil de végétation h < 2 m mis à 0** : c'est une convention Bouvier,
  elle doit être paramétrable (`h_min`) pour les tests de sensibilité.
- **fCover > 90%** : filtre actuel pour exclure gaps et lisières. Doit être
  paramétrable pour les tests de sensibilité.
- **PROSAIL inversion** : utiliser `prosail::Invert_PROSAIL` ou l'équivalent
  courant du package. Ne pas réimplémenter l'inversion. Vérifier la version
  du package utilisée et la documenter.
- **Coordinate system** : tout en UTM 31N (EPSG:32631). Les données Lambert-93
  brutes doivent être reprojetées. Ne jamais mélanger les deux.

## Communication

- Je suis francophone, Claude Code peut me répondre en français ou en
  anglais selon ce qui est le plus clair pour le sujet technique.
- Pour le code lui-même (commentaires, docstrings, noms), tout en anglais.
- Pour les messages de commit, en anglais.
- Pour les rapports d'analyse (`INVENTORY.md`, résumés), français ok.
- Quand Claude Code n'est pas sûr d'un choix méthodologique scientifique,
  il demande avant d'agir plutôt que d'inventer. Exemples : quelle métrique
  Pareto retenir, quel seuil pour un test, quel package si ambiguïté.
- Quand Claude Code rencontre un script existant dont la logique est peu
  claire, il me demande son rôle avant de le refactoriser.

## Checklist avant toute session longue

Avant une session de travail >1h, Claude Code doit :
1. Confirmer quelle branche git est active.
2. Confirmer quel module ou analyse est la cible de la session.
3. Lister les fichiers qu'il prévoit de créer ou modifier.
4. Attendre mon go avant de commencer.
