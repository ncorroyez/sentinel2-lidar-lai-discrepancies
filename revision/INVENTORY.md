# INVENTORY.md — Cartographie du code legacy

> Ce fichier est rempli par Claude Code à l'étape 1 du workflow de révision.
> Objectif : mapper chaque script existant aux figures et tables du papier
> soumis (RSE-D-25-04417), identifier les dépendances, et repérer les
> redondances / scripts obsolètes avant le refactoring.
>
> **Règle** : ne rien modifier pendant cette phase. Uniquement lire et
> documenter. Toute ambiguïté → demander à Nathan avant de trancher.

---

## 1. Résumé exécutif

*À remplir en dernier, après l'inventaire complet.*

- Nombre total de scripts R inspectés :
- Nombre de scripts actifs (utilisés dans le pipeline du papier) :
- Nombre de scripts legacy / obsolètes / redondants :
- Principales sources de fragilité identifiées :
- Ordre de refactoring recommandé :

---

## 2. Mapping figures du papier → scripts sources

Pour chaque figure et table du manuscrit soumis, identifier le(s) script(s)
qui la produit. Si incertain, mettre `?` et lister les candidats.

### Figures principales

| Figure | Description courte | Script(s) source | Fonctions clés appelées | Données d'entrée | Statut |
|--------|--------------------|------------------|-------------------------|------------------|--------|
| Fig. 1 | Maps LAI_ALS et LAI_S2 sur les 3 forêts | | | | |
| Fig. 2 | Workflow méthodologique | | | | |
| Fig. 3 | PAD profiles CHM vs DTM (exemple Aigoual) | | | | |
| Fig. 4 | Scatter LAI_S2_ATBD vs LAI_ALS par site | | | | |
| Fig. 5 | Corrélation vs depth (CHM/DTM, 4 panels) | | | | |
| Fig. 6 | Proportion configs PROSAIL dans Pareto front 1 | | | | |
| Fig. 7 | Évolution r/RMSE/bias/slope vs hétérogénéité | | | | |

### Tables principales

| Table | Description courte | Script source | Statut |
|-------|--------------------|---------------|--------|
| Table 1 | Dates d'acquisition ALS et S2 | (manuel) | |
| Table 2 | Distributions paramètres PROSAIL testées | (manuel) | |
| Table 3 | d_opt et r_CHM par site | | |
| Table 4 | d_opt et r_DTM par site | | |
| Table 5 | r_uniform, r_random, r_all par site | | |
| Table 6 | Performances Pareto front 1 vs ATBD | | |

### Figures supplémentaires (Appendix A.1 à A.14)

| Figure | Description | Script source | Statut |
|--------|-------------|---------------|--------|
| A.1 | Distribution des max heights par site | | |
| A.2 | PAD profiles CHM vs DTM tous sites | | |
| A.3 | Vérification ATBD vs SNAP | | |
| A.4 | Représentativité des classes d'hétérogénéité | | |
| A.5 | Sensibilité dmin (10/15/20 m) avec r_random | | |
| A.6 | (à identifier dans le manuscrit) | | |
| A.7 | Pareto fronts par site | | |
| A.8 | Scatter hétérogénéité Aigoual | | |
| A.9 | Scatter hétérogénéité Blois | | |
| A.10 | Scatter hétérogénéité Mormal | | |
| A.11 | Scatter hétérogénéité sites combined | | |
| A.12 | Distributions paramètres PROSAIL | | |
| A.13 | Relations LAI_ALS / CHM_SD / Max Height | | |
| A.14 | Histogrammes réflectances B03 B04 B08 | | |

---

## 3. Inventaire script par script

Pour chaque script R pertinent, remplir une fiche. Priorité aux dossiers :
`02_CODES/LiDAR/`, `02_CODES/Sentinel_2/`, `02_CODES/libraries/`,
`PROSAIL-Optimization/02_CODES/`. Ignorer `02_CODES/archive/` et
`02_CODES/other/` sauf si un script actif y fait référence.

### Template de fiche

```
#### `chemin/relatif/script.R`

- **Rôle** : (1-2 phrases)
- **Lu par** : (autres scripts qui le sourcent)
- **Source** : (scripts ou fichiers de libraries qu'il appelle)
- **Entrées** : (fichiers de données lus, depuis quel dossier)
- **Sorties** : (fichiers écrits, vers quel dossier)
- **Figures/tables produites** : (renvoi à la table section 2)
- **Dépendances packages** :
- **État du code** : clean / fonctionnel mais brouillon / cassé / obsolète
- **Candidat refactoring** : oui / non / à discuter
- **Notes** : (redondance avec autre script, TODO dans le code, questions)
```

### 3.1 Scripts LiDAR — `02_CODES/LiDAR/`

---

#### `02_CODES/LiDAR/0.change_filenames_L93_UTM31N.R`

- **Rôle** : Renomme les fichiers LAS en copiant leurs coordonnées L93 vers des coordonnées UTM31N dans le nom de fichier (sans reprojection des points, juste renommage).
- **Lu par** : aucun autre script
- **Source** : aucune library interne
- **Entrées** : `/media/corroyez/MyPassport/01_DATA/Aigoual/LiDAR/1-las_l93/` (chemin absolu)
- **Sorties** : `/media/corroyez/MyPassport/01_DATA/Aigoual/LiDAR/3-las_normalized_utm/`
- **Figures/tables produites** : aucune
- **Dépendances packages** : `sf`, `dplyr`, `stringr`
- **État du code** : fonctionnel mais brouillon (chemin absolu hardcodé, Aigoual seulement)
- **Candidat refactoring** : non (utilitaire de préparation one-shot)
- **Notes** : Complémentaire à `0_convert_l93_into_utm.R`. Fait la reprojection spatiale du *nom* de fichier (coordonnées L93→UTM31N dans le basename), pas du nuage de points. À distinguer de `0_convert_l93_into_utm.R` qui fait la vraie reprojection des points via LASTools. Dernière mise à jour 2025-12-02 (récent).

---

#### `02_CODES/LiDAR/0_convert_l93_into_utm.R`

- **Rôle** : Reprojection des fichiers LAS de Lambert-93 vers UTM 31N via LASTools (`las2las64`), avec renommage automatique en coordonnées UTM.
- **Lu par** : aucun autre script
- **Source** : aucune library interne
- **Entrées** : `/media/corroyez/MyPassport/01_DATA/{site}/LiDAR/leaf_off/1-las_l93/` (chemin absolu)
- **Sorties** : `/media/corroyez/MyPassport/01_DATA/{site}/LiDAR/leaf_off/2-las_utm/`
- **Figures/tables produites** : aucune
- **Dépendances packages** : `tools`, `lidR`
- **État du code** : fonctionnel mais brouillon (chemin LASTools hardcodé `/home/corroyez/Downloads/LAStools/bin/`, site = "Reine" en dur)
- **Candidat refactoring** : non (utilitaire one-shot, dépend de LASTools externe)
- **Notes** : Nécessite LASTools installé localement. La version commentée en bas utilise `sf` + lidR sans LASTools (approche alternative abandonnée). Site actif = "Reine" (pas un site du papier).

---

#### `02_CODES/LiDAR/0_convert_l93_into_utm_mask.R`

- **Rôle** : Reprojection des shapefiles BDForêt de Lambert-93 (EPSG:2154) vers UTM 31N (EPSG:32631), sauvegarde en GeoPackage.
- **Lu par** : aucun autre script identifié
- **Source** : aucune library interne
- **Entrées** : `/media/corroyez/MyPassport/01_DATA/{site}/Geo_Files/bdforet_2.shp` (chemin absolu)
- **Sorties** : `/media/corroyez/MyPassport/01_DATA/{site}/Geo_Files/bdforet_2_utm.gpkg`
- **Figures/tables produites** : aucune
- **Dépendances packages** : `sf`
- **État du code** : clean (simple, court, lisible)
- **Candidat refactoring** : non (one-shot, déjà exécuté)
- **Notes** : Traite les 3 sites du papier (Aigoual, Blois, Mormal). Étape préalable indispensable à `2_create_vegetation_forest_masks.R`.

---

#### `02_CODES/LiDAR/0_normalize_heights.R`

- **Rôle** : Normalisation DSM des nuages de points LAS (via `end_dsm_norm` en catalog_map) pour les 3 sites du papier. Produit des LAS normalisés dans `4-las_normalized_utm_dsm_chunked/`.
- **Lu par** : aucun autre script
- **Source** : `functions_lidar.R`, `functions_plots.R`, `functions_chm.R`
- **Entrées** : `/media/corroyez/MyPassport/01_DATA/{site}/LiDAR/2-las_utm/` + grilles S2 depuis `03_RESULTS/`
- **Sorties** : `/media/corroyez/MyPassport/01_DATA/{site}/LiDAR/4-las_normalized_utm_dsm_chunked/`
- **Figures/tables produites** : aucune
- **Dépendances packages** : `lidR`, `terra`, `future`, `dplyr`, `stringr`, `raster`, `rgl`, `viridis`
- **État du code** : fonctionnel mais brouillon (chemins absolus hardcodés, `sites <- c("Blois")` en dur, très nombreux blocs commentés d'exploration)
- **Candidat refactoring** : oui (contient la logique de normalisation DSM — `end_dsm_norm` — critique pour le pipeline)
- **Notes** : Appelle `end_dsm_norm` défini dans `functions_lidar.R`. Le `setwd()` via `rstudioapi` est une fragilité. Les sites Aigoual et Mormal sont commentés. À comparer avec `0tmp_normalize_heights.R`.

---

#### `02_CODES/LiDAR/0tmp_normalize_heights.R`

- **Rôle** : Version antérieure de `0_normalize_heights.R`, même date (2024-03-15). Contient un `stop()` à la ligne 126 qui empêche son exécution au-delà de la configuration catalog.
- **Lu par** : aucun
- **Source** : `functions_lidar.R`, `functions_plots.R`, `functions_chm.R`
- **Entrées** : `/media/corroyez/My Passport/01_DATA/` (chemin avec espace — ancien disque)
- **Sorties** : `/media/corroyez/My Passport/01_DATA/{site}/LiDAR/4-las_normalized_utm_dsm_chunked/`
- **Figures/tables produites** : aucune (bloqué par `stop()`)
- **Dépendances packages** : idem `0_normalize_heights.R` + `conflicted`, `rlang`
- **État du code** : cassé (`stop()` intentionnel, brouillon d'exploration)
- **Candidat refactoring** : non (doublon de `0_normalize_heights.R`, clairement obsolète)
- **Notes** : Seule différence notable : chemin disque avec espace ("My Passport" vs "MyPassport"), site = Aigoual. **Candidat à archiver/supprimer.**

---

#### `02_CODES/LiDAR/0tmp3_calculate_lidar_metrics.R`

- **Rôle** : Version "tmp" de `3_calculate_lidar_metrics.R`, même header (2024-07-22). Utilise `lidRmetrics` (commenté dans la version principale). Structure très similaire.
- **Lu par** : aucun
- **Source** : `functions_lidar.R`, `functions_create_masks.R`, `functions_plots.R`
- **Entrées** : grilles S2 + nuages LAS normalisés depuis `/media/corroyez/My Passport/01_DATA/`
- **Sorties** : `03_RESULTS/{site}/Metrics/`
- **Figures/tables produites** : ?
- **Dépendances packages** : idem `3_calculate_lidar_metrics.R` + `lidRmetrics`
- **État du code** : brouillon (préfixe "tmp", `lidRmetrics` actif)
- **Candidat refactoring** : non (doublon, explorer si des fonctions de `lidRmetrics` ont été intégrées ailleurs)
- **Notes** : **Question 1 (section 8).** Le préfixe `0tmp3` est inhabituel. À clarifier avec Nathan.

---

#### `02_CODES/LiDAR/0_good_parcel_viz.R`

- **Rôle** : Script de visualisation exploratoire : identifie les pixels Mormal à forte std et faible mean height, puis affiche le nuage de points LAS autour de ces coordonnées.
- **Lu par** : aucun
- **Source** : aucune library interne
- **Entrées** : `03_RESULTS/Mormal/LiDAR/Heterogeneity_Masks/` + `/media/corroyez/MyPassport/01_DATA/Mormal/LiDAR/2-las_utm/`
- **Sorties** : aucune (visualisation écran)
- **Figures/tables produites** : aucune dans le pipeline
- **Dépendances packages** : `terra`, `lidR`, `sf`, `rgl`, `viridis`
- **État du code** : brouillon (header dit "3.calculate_lidar_metrics.R" — copier-coller non mis à jour)
- **Candidat refactoring** : non (exploration one-shot)
- **Notes** : Header trompeur. Fichier trop large pour un simple script de viz (>10k tokens) — contient probablement beaucoup de blocs commentés.

---

#### `02_CODES/LiDAR/0_viz_las.R`

- **Rôle** : Visualise des nuages de points LAS centrés autour de placettes de terrain (coordonnées issues d'un GeoJSON), en affichant les métriques associées (LAI, LCV, etc.) dans le nom de fichier PNG.
- **Lu par** : aucun
- **Source** : `functions_lidar.R`, `functions_general_tools.R`, `functions_plots.R`
- **Entrées** : `03_RESULTS/{site}/data_utm31n.geojson` + résultats métriques + `/media/corroyez/My Passport/01_DATA/{site}/LiDAR/3-las_normalized_utm/`
- **Sorties** : `04_FIGURES/{site}/clouds/*.png`
- **Figures/tables produites** : figures d'exploration (pas dans le papier)
- **Dépendances packages** : `lidR`, `terra`, `data.table`, `geojsonio`, `lmom`, `rgl`
- **État du code** : brouillon (header dit "3.calculate_lidar_metrics.R" — copier-coller non mis à jour ; utilise `rgl.snapshot`)
- **Candidat refactoring** : non
- **Notes** : Utilise `load_metrics_with_pad()` définie dans `functions_general_tools.R` (hors scope de cette session). Header trompeur identique à `0_good_parcel_viz.R`.

---

#### `02_CODES/LiDAR/1_create_chm.R`

- **Rôle** : Crée le DTM, DSM et CHM à 1 m de résolution pour chaque site, en appelant `create_chm()` de `functions_chm.R`.
- **Lu par** : aucun autre script
- **Source** : `functions_chm.R`
- **Entrées** : `/media/corroyez/MyPassport/01_DATA/{site}/LiDAR/{leaf_state}/2-las_utm/`
- **Sorties** : `03_RESULTS/{site}/LiDAR/{leaf_state}/dtm/`, `dsm/`, `chm/`
- **Figures/tables produites** : aucune directement (produit des rasters intermédiaires critiques)
- **Dépendances packages** : `lidR`, `terra`, `future`, `future.apply`
- **État du code** : fonctionnel mais brouillon (sites actifs = `c("Reine", "Hayes")` — PAS les sites du papier ; `leaf_state = "leaf_off"`)
- **Candidat refactoring** : oui (wrappeur du pipeline, doit gérer les 3 sites paper)
- **Notes** : **Question 2 (section 8).** Sites actuellement configurés pour Hayes et Reine (hors papier). Les résultats pour Aigoual/Blois/Mormal sont probablement déjà dans `03_RESULTS/`. Voir `1_create_chm (copy).R` pour la version antérieure avec ancienne signature de `create_chm`.

---

#### `02_CODES/LiDAR/1_create_chm (copy).R`

- **Rôle** : Version antérieure de `1_create_chm.R`, utilisant l'ancienne signature de `create_chm(data_path, site, results_path, resolution)` (sans `s2_rast` ni `leaf_state`).
- **Lu par** : aucun
- **Source** : `functions_chm.R`
- **Entrées** : `../../01_DATA/` (chemin relatif — ancienne convention)
- **Sorties** : `../../03_RESULTS/{site}/LiDAR/`
- **Figures/tables produites** : aucune
- **Dépendances packages** : `lidR`, `raster`, `terra`
- **État du code** : obsolète (signature incompatible avec `functions_chm.R` actuel)
- **Candidat refactoring** : non (**candidat à supprimer**)
- **Notes** : Header dit "0_create_chm.R" (copier-coller). Remplacé par `1_create_chm.R`.

---

#### `02_CODES/LiDAR/2_create_vegetation_forest_masks.R`

- **Rôle** : Crée les masques de composition forestière (feuillu/résineux, strict/flex) à 10 m à partir des shapefiles BDForêt en appelant `create_vegetation_forest_mask()`.
- **Lu par** : aucun autre script identifié
- **Source** : `functions_sentinel_2.R`, `functions_plots.R`, `functions_lidar.R`, `functions_create_masks.R`
- **Entrées** : `/media/corroyez/My Passport/01_DATA/{site}/Shapefiles/` + `03_RESULTS/{site}/`
- **Sorties** : `03_RESULTS/{site}/LiDAR/Heterogeneity_Masks/`
- **Figures/tables produites** : aucune directement (masques utilisés par scripts 3–6)
- **Dépendances packages** : `lidR`, `raster`, `terra`, `viridis`, `future`
- **État du code** : fonctionnel (traite les 3 sites)
- **Candidat refactoring** : oui (simplifié par refactoring de `functions_create_masks.R`)
- **Notes** : Chemin avec espace "My Passport". Source `functions_sentinel_2.R` qui semble hors scope LiDAR pur.

---

#### `02_CODES/LiDAR/2.calculate_25m_metrics.R`

- **Rôle** : Calcule les métriques LiDAR (LAI, LAD profiles, CHM std, etc.) à résolution 25 m sur les 3 sites, variante de `3_calculate_lidar_metrics.R`.
- **Lu par** : aucun
- **Source** : `functions_lidar.R`, `functions_create_masks.R`, `functions_plots.R`
- **Entrées** : grilles S2 + nuages LAS normalisés depuis `/media/corroyez/MyPassport/01_DATA/`
- **Sorties** : `03_RESULTS/{site}/Metrics/`
- **Figures/tables produites** : ?
- **Dépendances packages** : `lidR`, `lidRmetrics`, `terra`, `lmom`, `data.table`, `solaR`, `raster`, `rgl`
- **État du code** : fonctionnel mais brouillon (résolution 25 m, chemin absolu)
- **Candidat refactoring** : à discuter
- **Notes** : Dernière mise à jour 2025-11-13 — **très récent**, potentiellement lié à la révision. Actif sur les 3 sites paper. Utilise `lidRmetrics` (package hors liste CLAUDE.md — **Question 3**). Résolution 25 m ≠ résolution S2 standard (10 m).

---

#### `02_CODES/LiDAR/3_calculate_lidar_metrics.R`

- **Rôle** : Script principal de calcul des métriques LiDAR : LAI (gap fraction), LAD profiles (DTM et DSM), CHM, std, CV, VCI, VDR, lcv, lskew, rumple, gap fraction, shade, slope — à 10 m sur la grille S2.
- **Lu par** : aucun autre script (script d'orchestration)
- **Source** : `functions_lidar.R`, `functions_create_masks.R`, `functions_plots.R`
- **Entrées** : grilles S2 + `{site}/LiDAR/2-las_utm/` + `{site}/LiDAR/3-las_normalized_utm/` depuis `/media/corroyez/MyPassport/01_DATA/`
- **Sorties** : `03_RESULTS/{site}/Metrics/Raw/*.tif` + masques appliqués
- **Figures/tables produites** : produit les métriques sources pour Figs. 3–7 et Tables 3–6
- **Dépendances packages** : `lidR`, `terra`, `lmom`, `data.table`, `solaR`, `raster`, `rgl`, `viridis`
- **État du code** : fonctionnel mais brouillon (sites actifs = `c("Hayes", "Reine")` — PAS les sites du papier, voir Question 2)
- **Candidat refactoring** : **oui, prioritaire** (cœur du pipeline)
- **Notes** : Fichier volumineux (>13k tokens). La version active traite Hayes et Reine, PAS les 3 sites du papier. Les résultats dans `03_RESULTS/` ont été produits par une version antérieure du script. À comparer avec `0tmp3_calculate_lidar_metrics.R`. `solaR` est utilisé pour l'analyse d'ombrage.

---

#### `02_CODES/LiDAR/3.calculate_lidar_metrics2.R`

- **Rôle** : Variante de `3_calculate_lidar_metrics.R` à résolution 20 m, pour Blois uniquement, utilisant une API plus propre via `functions_lidar2.R` (`myPAI`, `myPAD_z1` avec z0=1 m).
- **Lu par** : aucun
- **Source** : `functions_lidar2.R`, `functions_create_masks.R`
- **Entrées** : grilles S2 Blois + `Blois/LiDAR/leaf_on/3-las_normalized_utm/`
- **Sorties** : `03_RESULTS/Blois/Metrics/Raw/` (hmax_p95, hmax_p98, fCover_dynamic, lai_z1, lad_profiles_z1)
- **Figures/tables produites** : ? (expérimental)
- **Dépendances packages** : `lidR`, `terra`
- **État du code** : fonctionnel mais limité (Blois seulement, 20 m)
- **Candidat refactoring** : à discuter — l'API `functions_lidar2.R` est plus propre que `functions_lidar.R` et mérite d'être étendue
- **Notes** : Utilise z0=1 m au lieu de z0=2 m (convention Bouvier). **Question 4.** Pas de `solaR`, pas de `lmom`. Uniquement `lidR` + `terra` — dépendances minimales, approche plus claire.

---

#### `02_CODES/LiDAR/3_arrange_ladstack.R`

- **Rôle** : Construit des profils PAD cumulatifs (intégration couche par couche du LADstack) pour les normalisations DSM et DTM, sur les 3 sites. Produit des rasters de PAD intégré à différentes profondeurs.
- **Lu par** : aucun autre script
- **Source** : `functions_lidar.R`, `functions_create_masks.R`, `functions_plots.R`
- **Entrées** : `03_RESULTS/{site}/Metrics/Deciduous_Only/ladstack_{dsm|dtm}.tif` + `lidarlai_res_10_m.tif` + `max_res_10_m.tif`
- **Sorties** : `03_RESULTS/{site}/Metrics/Deciduous_Only/testPADs/PAD_Profiles_{DSM|DTM}_{keepTrees|Above20}/`
- **Figures/tables produites** : intermédiaires pour Fig. 5 (corrélation vs depth), Tables 3–4
- **Dépendances packages** : `lidR`, `lidRmetrics`, `terra`, `data.table`, `solaR`, `lmom`
- **État du code** : fonctionnel (traite 3 sites), dernière MAJ 2025-03-10
- **Candidat refactoring** : oui (logique de cumsum par couche à nettoyer)
- **Notes** : Header dit "3b.arrange_ladstack.R" (copier-coller incorrect depuis `3b`). Traite les 3 sites. Distinguer de `3b_arrange_ladstack.R`. Utilise `lidRmetrics` et `solaR` (hors liste CLAUDE.md — **Question 3**).

---

#### `02_CODES/LiDAR/3b_arrange_ladstack.R`

- **Rôle** : Variante de `3_arrange_ladstack.R` pour Aigoual uniquement, utilisant des fichiers ladstack avec suffixe `_0408` (vraisemblablement une acquisition spécifique).
- **Lu par** : aucun
- **Source** : `functions_lidar.R`, `functions_create_masks.R`, `functions_plots.R`
- **Entrées** : `03_RESULTS/Aigoual/Metrics/Deciduous_Only/ladstack_{dsm|dtm}_0408.tif`
- **Sorties** : `03_RESULTS/Aigoual/Metrics/Deciduous_Only/testPADs/PAD_Profiles_*/3107/`
- **Figures/tables produites** : intermédiaires (Aigoual seulement)
- **Dépendances packages** : idem `3_arrange_ladstack.R`
- **État du code** : brouillon (Aigoual only, sous-dossier "3107" mystérieux)
- **Candidat refactoring** : non (cas particulier à clarifier avant)
- **Notes** : Header identique à `3_arrange_ladstack.R` ("3b.arrange_ladstack.R"). Le suffixe `_0408` et le sous-dossier `3107` sont cryptiques. **Question 5.**

---

#### `02_CODES/LiDAR/3_shadows_analysis.R`

- **Rôle** : Calcule l'analyse d'ombres portées sur le DSM pour estimer l'impact de l'ombrage sur les réflectances S2.
- **Lu par** : aucun
- **Source** : `functions_shadows_analysis.R`
- **Entrées** : `../../01_DATA/{site}/` (chemin relatif)
- **Sorties** : non identifiées (dans `perform_shadows_analysis`)
- **Figures/tables produites** : ? (non dans les figures principales identifiées)
- **Dépendances packages** : `lidR`, `data.table`, `terra`, `viridis`, `rayshader`
- **État du code** : fonctionnel mais vieux (2023-05-28)
- **Candidat refactoring** : non (script exploratoire, probablement obsolète)
- **Notes** : Utilise `rayshader` (non dans la liste CLAUDE.md). Site = Mormal uniquement. Antérieur à la version finale du pipeline. Probablement du code d'exploration.

---

#### `02_CODES/LiDAR/4_prepare_heterogeneity_depth_analysis.R`

- **Rôle** : Prépare les masques de quantiles d'hétérogénéité (par déciles ou intervalles égaux) sur les métriques CHM (std, cv, mean, rumple, etc.) pour les 3 sites.
- **Lu par** : aucun (orchestre `create_heterogeneity_quantiles`)
- **Source** : `functions_plots.R`, `functions_create_masks.R`
- **Entrées** : `03_RESULTS/{site}/Metrics/{composition_mask}/*.tif` (std, cv, mean, max, lcv, lskew, vci, vdr, shade, slope, gap_fraction, dtm)
- **Sorties** : `03_RESULTS/{site}/LiDAR/Heterogeneity_Masks/Quantiles/{mask}/{metric}/`
- **Figures/tables produites** : intermédiaires pour Figs. 7, A.4–A.11
- **Dépendances packages** : `lidR`, `raster`, `terra`, `viridis`, `future`
- **État du code** : fonctionnel (3 sites), masques actifs : `"Not_Masked"` only (les autres commentés)
- **Candidat refactoring** : oui (très long if/else pour identifier `heterogeneity_metric_name` — à remplacer)
- **Notes** : La métrique active est `shade` (toutes les autres commentées). Le masque `composition_mask = "Not_Masked"` uniquement. La variable `shade` est assignée deux fois par erreur : `shade <- gap_fraction <- terra::rast(...)`. Bug potentiel.

---

#### `02_CODES/LiDAR/5_heterogeneity_depth_analysis.R`

- **Rôle** : Analyse combinée hétérogénéité × profondeur optique : calcule les corrélations LAI_ALS vs LAI_S2 par classe d'hétérogénéité et par profondeur de canopée. Produit les figures clés du papier.
- **Lu par** : aucun
- **Source** : `functions_plots.R`, `functions_heterogeneity_depth_analysis.R`
- **Entrées** : `03_RESULTS/{site}/LiDAR/Heterogeneity_Masks/Quantiles/` + `Metrics/{mask}/PAD_Profiles_updated/`
- **Sorties** : `04_FIGURES/{site}/Vegetation_Heterogeneity_{Depth}_Analysis/`
- **Figures/tables produites** : **Fig. 5, Fig. 7** (corrélation vs depth, évolution r/RMSE vs hétérogénéité)
- **Dépendances packages** : `lidR`, `terra`, `viridis`, `future`
- **État du code** : fonctionnel mais `sites <- "Blois"` en dur (les 3 sites sont commentés)
- **Candidat refactoring** : **oui, prioritaire** (produit les figures centrales du papier)
- **Notes** : Trois boucles principales : (1) hétérogénéité seule, (2) profondeur seule, (3) combinée. La boucle CV_LAD est séparée. Délègue tout le travail à `functions_heterogeneity_depth_analysis.R` (hors scope de cette session).

---

#### `02_CODES/LiDAR/5_correlation_per_PAI_profiles.R`

- **Rôle** : Script exploratoire précoce de corrélation LAI_LiDAR vs LAI_S2 par classes de mean height. Probablement précurseur de `5_heterogeneity_depth_analysis.R`.
- **Lu par** : aucun
- **Source** : `../functions_plots.R`, `functions_lidar.R` (chemins *incorrects* — manque le sous-dossier `libraries/`)
- **Entrées** : `03_RESULTS/Blois/LiDAR/Heterogeneity_Masks/` (Blois seulement, format ENVI)
- **Sorties** : fichiers de plot
- **Figures/tables produites** : aucune dans le papier final
- **Dépendances packages** : `lidR`, `terra`, `viridis`, `future`
- **État du code** : **cassé** (chemins `source()` incorrects)
- **Candidat refactoring** : non (**candidat à archiver**)
- **Notes** : Header dit "3.correlation_per_PAI_profiles.R". `source("../functions_plots.R")` cherche au mauvais niveau. Lit des fichiers ENVI (`*.envi`) — format plus utilisé dans les scripts récents. Script très ancien.

---

#### `02_CODES/LiDAR/6a_explain_heterogeneity_unique_site.R`

- **Rôle** : Analyse Random Forest pour identifier les variables LiDAR expliquant l'hétérogénéité LAI_ALS vs LAI_S2, par site (approche site-spécifique). Co-écrit avec J.-B. Féret.
- **Lu par** : aucun
- **Source** : `functions_JBF.R`, `functions_plots.R`
- **Entrées** : `03_RESULTS/{site}/` (métriques chargées via fonctions JBF)
- **Sorties** : `03_RESULTS/{site}/` + figures
- **Figures/tables produites** : Fig. 7 (facteurs explicatifs), Tables 5–6 (performances par site)
- **Dépendances packages** : `foreach`, `doFuture`, `future`, `parallel`, `doParallel`, `stringr`, `randomForest`
- **État du code** : fonctionnel (actif sur `c('Aigoual')` — les 3 sites commentés)
- **Candidat refactoring** : oui (dépend de `functions_JBF.R` hors scope)
- **Notes** : Utilise `randomForest` et `doParallel` (hors liste CLAUDE.md). Co-auteur JBF. La plupart de la logique est dans `functions_JBF.R`.

---

#### `02_CODES/LiDAR/6b_explain_heterogeneity_mix_sites.R`

- **Rôle** : Analyse RF multi-sites (train sur 2 sites, test sur le 3e) pour évaluer la généralisation du modèle hétérogénéité → LAI discordance.
- **Lu par** : aucun
- **Source** : `functions_JBF.R`, `functions_plots.R`
- **Entrées** : `03_RESULTS/{site}/`
- **Sorties** : `03_RESULTS/{site}/` + figures
- **Figures/tables produites** : Fig. 7 (r_uniform, r_random, r_all), Table 5
- **Dépendances packages** : `foreach`, `doFuture`, `future`, `parallel`, `stringr`, `randomForest`
- **État du code** : fonctionnel (3 sites, combinaisons train/test)
- **Candidat refactoring** : oui
- **Notes** : Même header que `6b2` — les deux fichiers ont `title: "6b_explain_heterogeneity_mix_sites.R"`. Pas de `pls`, `test_sample_size = 50000`.

---

#### `02_CODES/LiDAR/6b2_explain_heterogeneity_mix_sites.R`

- **Rôle** : Variante de `6b` avec combinaisons train/test différentes (train = tous les 3 sites, test = chaque site individuellement). Ajoute le package `pls`. `test_sample_size = 3000` (beaucoup plus petit).
- **Lu par** : aucun
- **Source** : `functions_JBF.R`, `functions_plots.R`
- **Entrées** : `03_RESULTS/{site}/`
- **Sorties** : `03_RESULTS/{site}/` + figures
- **Figures/tables produites** : **Question 6** — laquelle de 6b ou 6b2 a produit les figures finales ?
- **Dépendances packages** : idem `6b` + `pls`
- **État du code** : brouillon (header identique à `6b` — doublon apparent)
- **Candidat refactoring** : à discuter après clarification
- **Notes** : Header strictement identique à `6b` (titre, auteurs, date). `test_sample_size` très différent (3000 vs 50000). `combinations` liste redéfinie avec train = all 3 sites.

---

#### `02_CODES/LiDAR/6c_explain_heterogeneity_stability.R`

- **Rôle** : Analyse de stabilité du modèle RF hétérogénéité : évalue la robustesse des variables sélectionnées à travers différentes initialisations (`set.seed`).
- **Lu par** : aucun
- **Source** : `functions_JBF.R`, `functions_plots.R`
- **Entrées** : `03_RESULTS/{site}/`
- **Sorties** : `03_RESULTS/{site}/` + figures
- **Figures/tables produites** : ? (figures de stabilité — probablement supplémentaires)
- **Dépendances packages** : `foreach`, `doFuture`, `future`, `parallel`, `doParallel`, `stringr`, `randomForest`
- **État du code** : fonctionnel (3 sites)
- **Candidat refactoring** : oui
- **Notes** : Même structure que `6a` mais avec `StabilityIndex`. Co-auteur JBF.

---

#### `02_CODES/LiDAR/max_hists.R`

- **Rôle** : Génère des histogrammes de distribution des hauteurs maximales (H_max) par site. Probablement Fig. A.1 du manuscrit.
- **Lu par** : aucun
- **Source** : aucune library interne
- **Entrées** : `03_RESULTS/{site}/Metrics/Deciduous_Only/max_res_10_m.tif`
- **Sorties** : figures (chemin non vu dans les 60 lignes lues)
- **Figures/tables produites** : vraisemblablement **Fig. A.1** (Distribution des max heights)
- **Dépendances packages** : `terra`, `ggplot2`, `dplyr`, `stringr`, `purrr`
- **État du code** : clean (court, lisible, pas de chemins absolus hardcodés)
- **Candidat refactoring** : non (script court autonome)
- **Notes** : Dernière MAJ 2025-08-25 — **très récent**. Utilise `here`-compatible `results_dir <- "../../03_RESULTS"`. Package `purrr` (hors liste CLAUDE.md).

---

#### `02_CODES/LiDAR/pad_boxplots.R`

- **Rôle** : Génère des boxplots de la distribution des valeurs PAD par couche de profondeur (PAD value vs depth).
- **Lu par** : aucun
- **Source** : aucune library interne
- **Entrées** : `03_RESULTS/{site}/Metrics/Deciduous_Only/PAD_Profiles_dsm_keepTrees/*.tif`
- **Sorties** : figures (non lues)
- **Figures/tables produites** : probablement **Fig. A.2** (PAD profiles tous sites)
- **Dépendances packages** : `terra`, `dplyr`, `stringr`, `ggplot2`, `purrr`
- **État du code** : clean (bien structuré, fonction encapsulée)
- **Candidat refactoring** : non
- **Notes** : Dernière MAJ 2025-07-23 — récent. Définit une fonction `make_pad_value_vs_depth_boxplot` inline (candidat à migrer dans `revision/R/`).

---

#### `02_CODES/LiDAR/main_viz_norms.R`

- **Rôle** : Visualise les nuages de points LAS normalisés (DTM et DSM) pour validation visuelle.
- **Lu par** : aucun
- **Source** : `functions_lidar.R`, `functions_plots.R`, `functions_chm.R`
- **Entrées** : `/media/corroyez/My Passport/01_DATA/{site}/LiDAR/2-las_utm/` et `3-las_normalized_utm/`
- **Sorties** : visualisations écran (non sauvegardées dans les 60 lignes lues)
- **Figures/tables produites** : aucune dans le papier
- **Dépendances packages** : `lidR`, `terra`, `rgl`, `viridis`, `raster`
- **État du code** : brouillon (exploration, site = Mormal seulement, chemin avec espace)
- **Candidat refactoring** : non
- **Notes** : Chemin disque ancien ("My Passport" avec espace).

---

#### `02_CODES/LiDAR/grid_metrics_final.Rmd`

- **Rôle** : Document R Markdown d'exploration LiDAR (calcul de métriques sur grille), par Marianne Laslier (août 2022).
- **Lu par** : aucun
- **Source** : aucune library interne (bibliothèques directement chargées dans le Rmd)
- **Entrées** : `D:\\LAS_Hauteur\\BLOIS` (chemin Windows absolu)
- **Sorties** : aucune (exploration)
- **Figures/tables produites** : aucune dans le papier
- **Dépendances packages** : `raster`, `rgdal`, `sp`, `lidR`, `tidyverse`
- **État du code** : **obsolète** (utilise `rgdal` déprécié, chemin Windows absolu, auteur externe)
- **Candidat refactoring** : non (**candidat à supprimer ou archiver**)
- **Notes** : Antérieur au pipeline. `rgdal` est retiré de CRAN depuis 2023.

### 3.2 Scripts Sentinel-2 — `02_CODES/Sentinel_2/`

*À remplir. Scripts prioritaires :*
- `Main_01_download_S2_sites.R`
- `Main_02_produceLAI.R` / `Main_02_produceLAI_v2.R`
- `0_verif_atbd_biophysical_toolbox.R` / `0.verif_atbd_biophysical_toolbox_v2.R`
- `1_mask_s2_ts.R`
- `3_train_predict_prosail.R` / `3.train_predict_prosail_ts.R`

> Note : idem, identifier les versions actives vs legacy.

### 3.3 Libraries — `02_CODES/libraries/`

*Inventaire des 5 fichiers demandés dans cette session. Les autres fichiers
(`functions_heterogeneity_depth_analysis.R`, `functions_JBF.R`, etc.)
sont à inventorier lors de sessions ultérieures.*

---

#### `02_CODES/libraries/functions_lidar.R`

- **Lu par** : `0_normalize_heights.R`, `0tmp_normalize_heights.R`, `0_viz_las.R`, `2_create_vegetation_forest_masks.R`, `3_calculate_lidar_metrics.R`, `0tmp3_calculate_lidar_metrics.R`, `3_arrange_ladstack.R`, `3b_arrange_ladstack.R`, `2.calculate_25m_metrics.R`, `main_viz_norms.R`
- **Dernière MAJ** : 2024-05-27
- **État** : fonctionnel mais brouillon (debug prints laissés dans `myPAD_dsm`, fonctions redondantes)

Fonctions :
- `VCI_local(z)` : Vertical Complexity Index local (wrap de `lidR::VCI`) ; utilisée dans métriques pixel
- `VCI_dsm(z, zmax=40, z_opt=5)` : VCI filtré sur les points proches du sommet DSM
- `VDR(z)` : Vertical Density Ratio = (zmax − median) / zmax
- `CV(z)` : Coefficient of Variation = sd/mean × 100
- `rumple_index_surface(las, res)` : Rumple Index sur surface points, via `pixel_metrics`
- `myPAI(z, zmin, k=0.5)` : Plant Area Index par gap fraction (Beer-Lambert) — version principale, z0 en argument `zmin`
- `wrapper_gap_fraction(z)` : profil de gap fraction DSM-shifted, retourne GF total (usage interne)
- `myPAD(z)` : LAD profiles + PAI cumulatifs via `LAD()`, couches z ∈ [2.5, 39.5], z0 = 2 implicite (convention Bouvier)
- `myPAD_dtm(z)` : idem avec shift DTM (z + 40 − max_z pour aligner au sommet)
- `myPAD_dtm2(z)` : variante avec z0=2 explicite et print() de debug — probablement obsolète
- `myPAD_dsm(z)` : LAD profiles normalisation DSM — **contient print(lad) debug non retiré**
- `myPAD_dsm_new(z)` : version épurée de myPAD_dsm, retourne uniquement LADs_list
- `myPAD_dtm_new(z, class)` : variante DTM utilisant la classification LAS pour le shift
- `end_lad_dtm(z, class)` : callback `catalog_map` pour profils LAD DTM-shifted, API épurée
- `end_lad_dsm(z)` : callback `catalog_map` pour profils LAD DSM
- `end_dtm_norm(las)` : callback `catalog_map` — normalise par DTM TIN, filtre Z>40, reclassifie Z<2 en sol
- `end_dsm_norm(las)` : callback `catalog_map` — normalise par DSM pitfree, shift sol, reclassifie Z<2
- `is.empty(x)` : test si LAS est vide (redéfini localement, conflit possible avec `lidR::is.empty`)
- `fill.na(x, i=5)` : fonction de remplissage NA pour `terra::focal`
- `perform_shadows_analysis(dsm, metrics_dir)` : analyse d'ombres portées (rayshader) — signature ne correspond pas à l'appel dans `3_shadows_analysis.R` (bug potentiel)
- `my_hillshade(dsm)` : génère une hillshade depuis un DSM
- `calculate_solar_angles(lat, long, date_time, slope, azimuth)` : angles solaires pour l'analyse d'ombres

**Notes** : `myPAD_dsm` contient `print(lad)` et `print(sum(...))` laissés en debug — impact performance en production. Plusieurs variantes de `myPAD_*` dont la version active pour le papier est à clarifier (Question 4, section 8). `perform_shadows_analysis` a une signature différente de l'appel dans `3_shadows_analysis.R`.

---

#### `02_CODES/libraries/functions_lidar2.R`

- **Lu par** : `3.calculate_lidar_metrics2.R` uniquement
- **Dernière MAJ** : non datée (pas de header complet)
- **État** : clean, API épurée

Fonctions :
- `myPAI(z, zmin=1, k=0.5)` : Plant Area Index — **même nom que dans `functions_lidar.R`** mais avec z0=1 m par défaut (au lieu de z0=2 m) et valeur par défaut sur zmin
- `myPAD_z1(z)` : LAD profiles via `LAD(z, z0=1)`, couches z ∈ [1.5, 39.5], retourne liste nommée `LAD_Layer_*`

**Notes** : Implémentation plus propre et plus compacte que `functions_lidar.R`. z0=1 m = écart avec la convention Bouvier (z0=2 m). **Question 4** : quel z0 doit être utilisé pour la révision ? Les deux fichiers définissent `myPAI` avec des signatures partiellement incompatibles — risque de collision si les deux sont sourcés dans un même script.

---

#### `02_CODES/libraries/functions_normalize_heights.R`

- **Lu par** : **aucun script identifié**
- **Dernière MAJ** : 2024-05-24 (header contient "functions_chm.R" — titre incorrect)
- **État** : **obsolète** — copie mal renommée de l'ancienne version de `functions_chm.R`

Fonctions (version ancienne de functions_chm.R) :
- `format_time(start, end)` : format du temps écoulé
- `print_raster_size(raster)` : dimensions du raster en chaîne
- `create_chm(data_path, site, results_path)` : ancienne signature **sans** `s2_rast`, `leaf_state`, `resolution`

**Notes** : Ce fichier a le header `title: "functions_chm.R"` — il s'agit d'une copie mal renommée de l'ancienne version de `functions_chm.R`. La signature `create_chm(data_path, site, results_path)` est incompatible avec `functions_chm.R` actuel (6 arguments). **Candidat à supprimer après confirmation.**

---

#### `02_CODES/libraries/functions_chm.R`

- **Lu par** : `0_normalize_heights.R`, `0tmp_normalize_heights.R`, `1_create_chm.R`, `main_viz_norms.R`
- **Dernière MAJ** : 2024-06-03
- **État** : fonctionnel, bien documenté (roxygen2)

Fonctions :
- `format_time(start, end)` : format du temps écoulé (identique à celle de `functions_normalize_heights.R`)
- `print_raster_size(raster)` : dimensions du raster en chaîne
- `dsm_tile(las, resolution)` : calcule DSM 1 m par tuile via `p2r` + focal mean (callback catalog)
- `create_chm(data_path, results_dir, site, s2_rast, leaf_state, resolution)` : pipeline complet DTM (TIN) + DSM (p2r) + CHM, gère les tuiles problématiques par site

**Notes** : `create_chm` a une signature à 6 paramètres — incompatible avec `1_create_chm (copy).R` (4 params). Les tuiles problématiques (Aigoual: LAS_745500_6337000, Blois: LAS_567000_6724500, Mormal: 4 tuiles) sont hardcodées dans une liste — à paramétrer lors du refactoring.

---

#### `02_CODES/libraries/functions_create_masks.R`

- **Lu par** : `2_create_vegetation_forest_masks.R`, `3_calculate_lidar_metrics.R`, `0tmp3_calculate_lidar_metrics.R`, `3_arrange_ladstack.R`, `3b_arrange_ladstack.R`, `4_prepare_heterogeneity_depth_analysis.R`, `3.calculate_lidar_metrics2.R`, `2.calculate_25m_metrics.R`
- **Dernière MAJ** : 2024-05-27
- **État** : fonctionnel, bien documenté (roxygen2 partiel)

Fonctions :
- `mask_forest_composition_area(site_edges, forest_composition, masks_dir, resolution)` : découpe BDForêt par emprise site, crée masques feuillu/résineux strict/flex, rasterise à `resolution` m
- `mask_site_edges(reflectance, site_edges_path, masks_dir, resolution)` : applique masque d'emprise site sur un raster de réflectance
- `mask_clouds(reflectance, cloud_path, site_edges, masks_dir, resolution)` : masque nuages sur raster réflectance
- `save_chm(chm, ...)` : sauvegarde CHM (détails non lus)
- `mask_chm_thresh(chm, ...)` : masque CHM sous un seuil (détails non lus)
- `project_chm_thresh(chm_thresh, ...)` : reprojection du CHM masqué (détails non lus)
- `mask_majority_project_chm_thresh(chm, ...)` : masque CHM par majorité + reprojection
- `create_vegetation_forest_mask(data_dir, site, shapefiles_dir, masks_dir, results_path, resolution)` : fonction principale appelée dans `2_create_vegetation_forest_masks.R`
- `create_heterogeneity_quantiles(error, ...)` : crée les masques de quantiles d'hétérogénéité (appelée dans `4_prepare_heterogeneity_depth_analysis.R`)
- `create_cv_lad_quantiles(error, ...)` : variante de `create_heterogeneity_quantiles` pour CV_LAD
- `apply_and_save_masks(raster, out_filename, masks_dir, metrics_dir)` : applique les masques de composition à un raster et sauvegarde les versions masquées

**Notes** : Fichier volumineux (>15k tokens). Utilise encore `raster::rasterize` dans `mask_forest_composition_area` (API dépréciée — devrait utiliser `terra::rasterize`). Les fonctions `save_chm`, `mask_chm_thresh`, `project_chm_thresh` n'ont pas été lues en détail.

### 3.4 Scripts PROSAIL — `PROSAIL-Optimization/02_CODES/`

> Ce sous-projet a sa propre arborescence `01_DATA/`, `02_CODES/`, `03_RESULTS/`.
> Tous les scripts `Main_*.R` sourcent l'intégralité des 8 fichiers du dossier
> `libraries/` via `lapply(list.files('libraries', full.names=T), source)`.
> Tous utilisent `setwd(dirname(rstudioapi::getSourceEditorContext()$path))` et
> des chemins absolus hardcodés vers `/home/corroyez/Documents/NC_Full/`.

---

#### Fiche — `Main_01_analysis_optim_depth.R`

| Champ | Valeur |
|-------|--------|
| **Rôle** | Analyse de la profondeur optimale d'intégration des profils LAD (d_opt). Pour chaque site et chaque seuil `max_h ∈ {10, 15, 20}` m, filtre les pixels où `max_height > max_h`, effectue un échantillonnage stratifié (n=5 000) par bins LAI_ALS, calcule la corrélation R² entre LAI_ALS_dopt et LAI_S2 selon la profondeur d'intégration. |
| **Sources** | `lapply(list.files('libraries', full.names=T), source)` |
| **Entrées** | Rasters LiDAR : `max`, `mean`, `lcv`, `lskew`, `vci`, `rumple`, `fCover` depuis `../01_DATA/{site}/LiDAR/` ; réflectances S2 B03/B04/B08 (chemins hardcodés) ; profils PAD empilés par couche |
| **Sorties** | Figures de corrélation R²(depth) pour chaque normalization × site ; CSV intermédiaires de métriques |
| **Figures papier** | Fig. 4 (sélection de d_opt) — **à confirmer** |
| **Packages hors CLAUDE.md** | `sgsR` (stratified sampling) |
| **État** | Lu partiellement (120 premières lignes sur ~450). Hardcodé pour `sites = c("Aigoual", "Blois", "Mormal")` |
| **Candidat refactoring** | Oui — chemins absolus, `setwd()`, valeurs `max_h` hardcodées, pas de `here::here()` |
| **Notes** | `max_h_values = c(10, 15, 20)` est une première approximation ; la valeur finale (d_opt) est déterminée par `Main_02`. Dépend en amont des profils PAD produits par `3_arrange_ladstack.R`. |

---

#### Fiche — `Main_02_analysis_best_prosail_config.R`

| Champ | Valeur |
|-------|--------|
| **Rôle** | Optimisation de la configuration PROSAIL LUT. Lit un grand CSV de résultats (`all_results_combined_LIDFa_lai_LMA_BROWN_N_CHL_psoil_q_Final13_05.csv`) contenant les scores R, NRMSE, slope pour chaque combinaison de paramètres PROSAIL testés. Calcule un score composite `Score = (norm_R + norm_NRMSE + norm_slope)/3`, identifie les meilleures configs par site (Best_Indiv) et cross-site (Common), effectue une analyse d'importance via Random Forest. |
| **Sources** | `lapply(list.files('libraries', full.names=T), source)` |
| **Entrées** | CSV `all_results_combined_LIDFa_lai_LMA_BROWN_N_CHL_psoil_q_Final13_05.csv` ; filtre `Norm="DSM_keepTrees"`, `Depth=5` (Aigoual/Blois), `Depth=9` (Mormal) |
| **Sorties** | Figures d'importance de paramètres (RF), distributions des scores, tableau configs optimales (ATBD / Common / Best_Indiv) |
| **Figures papier** | Fig. 3 (comparaison configs PROSAIL) — **à confirmer** |
| **Packages hors CLAUDE.md** | `rPref` (Pareto-like preferences), `randomForest`, `pdp` (partial dependence), `broom`, `forcats`, `Metrics`, `RColorBrewer`, `stringr` |
| **État** | Lu partiellement (120 premières lignes sur ~540). Stratégie nommée `"LIDFa_lai_LMA_BROWN_N_CHL_psoil_q"`. 8 paramètres testés. |
| **Candidat refactoring** | Oui — version précédente, remplacée fonctionnellement par `Main_02b` |
| **Notes** | `Depth=5` pour Aigoual/Blois correspond à h=5 m d'intégration du bas de canopée — **différent de la convention d_opt papier (35.5 m depuis le haut)**. Sens à clarifier. Beaucoup de packages non listés dans CLAUDE.md. |

---

#### Fiche — `Main_02b_final_analysis_best_prosail_config.R`

| Champ | Valeur |
|-------|--------|
| **Rôle** | Version finale de l'optimisation PROSAIL. Même logique que `Main_02` mais avec une stratégie réduite (4 paramètres : ALA, lai, LMA, BROWN au lieu de 8), profondeur d'intégration unifiée (`Depth=4` pour tous sites), renommage LIDFa→ALA. Vraisemblablement la version utilisée pour les figures finales du papier soumis. |
| **Sources** | `lapply(list.files('libraries', full.names=T), source)` |
| **Entrées** | CSV `all_results_combined_LIDFa_lai_LMA_BROWN_Agg_10m_16_09_Try.csv` ; filtre `Norm="DSM_keepTrees"`, `Depth=4` (tous sites) |
| **Sorties** | Mêmes types de figures que `Main_02` (importance RF, distributions, configs optimales) ; vraisemblablement les figures finales |
| **Figures papier** | Probablement Fig. 3 version finale — **à confirmer** |
| **Packages hors CLAUDE.md** | `rPref`, `randomForest`, `pdp`, `broom`, `forcats`, `Metrics`, `RColorBrewer`, `stringr` |
| **État** | Lu partiellement (120 premières lignes sur ~430). `name_strategy = "LIDFa_lai_LMA_BROWN_Agg_10m"` |
| **Candidat refactoring** | Oui — **version active** pour la révision. Priorité haute. |
| **Notes** | `Depth=4` pour tous les sites = convention différente de `Main_02`. LIDFa renommé ALA (Average Leaf Angle), cohérent avec la biophysique PROSPECT-D. Deux versions de CSV impliquent deux rounds d'optimisation LUT distincts. |

---

#### Fiche — `Main_03_student.R`

| Champ | Valeur |
|-------|--------|
| **Rôle** | Scatter plots de comparaison LAI_ALS vs LAI_S2 pour les 3 configurations PROSAIL (ATBD, Common, Best_Indiv) × 2 versions LAI_ALS (full vs dopt). Calcule R, RMSE, slope, bias et tests t appariés par site. |
| **Sources** | `lapply(list.files('libraries', full.names=T), source)` |
| **Entrées** | Rasters LAI_ALS et LAI_S2 par site. d_opt hardcodé : `h="35.5"` pour Aigoual et Blois, `h="31.5"` pour Mormal |
| **Sorties** | Scatter plots LAI_ALS vs LAI_S2 × 3 configs × 3 sites (sauvegardés en PNG/PDF) |
| **Figures papier** | Prédécesseur de `Main_03_student2.R`. Probablement pas les figures finales. |
| **Packages hors CLAUDE.md** | `Metrics`, `broom`, `ggpubr` (dans student2 uniquement) |
| **État** | Lu en entier (369 lignes). `stop()` à la ligne 351 → exécution interrompue avant la fin. |
| **Candidat refactoring** | Non — remplacé par `Main_03_student2.R` |
| **Notes** | d_opt = 35.5 m pour Aigoual et Blois, 31.5 m pour Mormal = intégration des LAD de h jusqu'à 40 m. Ces valeurs sont hardcodées comme chaînes de caractères (`h="35.5"`). |

---

#### Fiche — `Main_03_student2.R`

| Champ | Valeur |
|-------|--------|
| **Rôle** | Version améliorée de `Main_03_student`. Ajoute le calcul du gain en % entre 3 paires de configurations : (ATBD_dopt vs ATBD_ALS), (Best_Indiv_dopt vs ATBD_dopt), (Best_Indiv_dopt vs ATBD_ALS). Test de variance avant la régression lm() pour robustesse. |
| **Sources** | `lapply(list.files('libraries', full.names=T), source)` |
| **Entrées** | Mêmes que `Main_03_student`. d_opt identiques (35.5 / 31.5 m). |
| **Sorties** | Scatter plots + tableau de gains en % entre configurations |
| **Figures papier** | Fig. 5 (comparaison LAI_ALS_full vs LAI_ALS_dopt vs LAI_S2 configs) — **version active** |
| **Packages hors CLAUDE.md** | `Metrics`, `broom`, `ggpubr` |
| **État** | Lu en entier (428 lignes). `stop()` à la ligne 313 — code après le stop() non exécuté. |
| **Candidat refactoring** | Oui — version active pour la révision. |
| **Notes** | `stop()` ligne 313 partitionne le script en deux blocs indépendants. d_opt hardcodés en chaînes de caractères. Utilise `ggpubr` (hors liste CLAUDE.md). |

---

#### Fiche — `Main_04_final_factors_tabs.R`

| Champ | Valeur |
|-------|--------|
| **Rôle** | Figures finales avec axes inversés par rapport à `Main_03_student2` (LAI_S2 en x, LAI_ALS en y). Ajoute une analyse par boxplots en stratifiant les pixels LAI_ALS < 3 vs ≥ 3, avec test de Wilcoxon (`ggpubr::stat_compare_means`). |
| **Sources** | `lapply(list.files('libraries', full.names=T), source)` |
| **Entrées** | Mêmes que `Main_03_student2`. |
| **Sorties** | Scatter plots axes inversés + boxplots stratifiés par classe LAI |
| **Figures papier** | Fig. 6 ou Fig. S* (boxplots par classe LAI) — **à confirmer** |
| **Packages hors CLAUDE.md** | `Metrics`, `broom`, `ggpubr` |
| **État** | Lu en entier (428 lignes). `stop()` à la ligne 313. |
| **Candidat refactoring** | Oui — version active pour la révision. |
| **Notes** | L'inversion des axes (LAI_ALS en y) est intentionnelle : le papier cadre LAI_ALS comme la référence. Le seuil LAI < 3 est hardcodé — à paramétrer. |

---

#### 3.4.1 Fonctions des libraries PROSAIL

##### `define_parm_combinations.R`

| Fonction | Signature | Rôle |
|----------|-----------|------|
| `define_parm_combinations` | `(parms2test, output_dir, name_strategy, overwrite=F)` | Construit la grille de simulation (`expand.grid`) pour tous les paramètres PROSAIL à tester en appelant `get_combination()` pour chacun. Sauvegarde en `.rds` dans `Simulations_Strategy/{name_strategy}/simulation_strategy.rds`. Bloque si le fichier existe déjà (`overwrite=F`). |

##### `get_combination.R`

| Fonction | Signature | Rôle |
|----------|-----------|------|
| `get_combination` | `(parm)` | Retourne la liste des distributions à tester pour un paramètre PROSAIL donné. Pour la stratégie active (4 paramètres) : LIDFa/ALA = 5 configs (ATBD + 4 gaussiennes), lai = 6 configs (ATBD + 3 gaussiennes + LiDAR_LAI + LiDAR_LAI_Best_Site_Depth), LMA = 3 configs, BROWN = 3 configs. Une ancienne version commentée teste aussi N, CHL, psoil, q (8 paramètres). |

##### `get_s2_angles.R`

| Fonction | Signature | Rôle |
|----------|-----------|------|
| `get_s2_angles` | `(path_angles, path_bbox, dateAcq=NULL)` | Lit les rasters d'angles S2 (SAA, SZA, VAA, VZA), recadre sur la bbox du site, retourne les min/max de chaque angle + l'angle azimutal relatif psi. Utilisé pour conditionner la géométrie d'acquisition dans la génération du LUT PROSAIL. |

##### `get_s2_samples.R`

| Fonction | Signature | Rôle |
|----------|-----------|------|
| `get_s2_samples` | `(aoi_path, S2_path, ..., nbSamples, site, method)` | Extrait des échantillons de pixels S2 selon 3 méthodes : `"random"` (spatSample + snap aux centres de pixels), `"stratified"` (stratification sur mean+lskew LAI via sgsR), `"stratified_uniform"` (stratification uniforme sur bins LAI, filtre max_height ∈ [10, 40] m). |
| `align_and_remove_na_for_aoi` | `(aoi_path, S2_path, lidar_lai_path, sampling_path, save=T)` | Aligne les rasters S2 et LiDAR, masque les pixels NA, retourne un AOI nettoyé. |
| `stratified_sampling` | `(mean_path, lskew_path, field_points_path, nbSamples)` | Stratification sur quantiles de mean + lskew via `sgsR::strat_breaks` + `sgsR::sample_strat`. |
| `stratified_sampling_uniform` | `(site, lidar_lai_path, max_path, nbSamples)` | Stratification uniforme sur bins LAI (cible 5 000 points équirépartis). Filtre max_height ∈ [10, 40] m. Inclut un `print(hist(...))` de débogage non nettoyé. |

##### `PROSAIL_train_sensitivity.R`

| Fonction | Signature | Rôle |
|----------|-----------|------|
| `PROSAIL_train_sensitivity` | `(simulations, geom_s2, lidar_lai_rast, ..., nbSamples, filename, nbCPU, S2BandSelect)` | Dispatcher d'entraînement LUT PROSAIL. Si `nbCPU=1` : séquentiel via `mapply`. Si `nbCPU>1` : parallèle via `future.apply::future_mapply` avec cluster `parallel`. Appelle `PROSAIL_train_sensitivity_subset` pour chaque ligne de la grille de simulation. |

##### `PROSAIL_train_sensitivity_subset.R`

| Fonction | Signature | Rôle |
|----------|-----------|------|
| `PROSAIL_train_sensitivity_subset` | `(simuset, filename, combinations, lidar_lai_rast, lidar_lai_site_rast, lidar_lai_common_rast, lidar_lai_3m_rast, geom_s2, nbSamples, p, S2BandSelect)` | Cœur de l'entraînement PROSAIL. Pour chaque simulation : (1) charge la distribution ATBD via `prosail::get_InputPROSAIL(atbd=TRUE)`, (2) remplace les paramètres selon leur distribution (uniforme, gaussienne tronquée, ou échantillonnage depuis les LAI LiDAR), (3) génère le LUT via `Generate_LUT_PROSAIL(SAILversion='4SAIL')`, (4) convole vers les bandes S2 via `applySensorCharacteristics`, (5) ajoute du bruit (`apply_noise_atbd`), (6) entraîne les SVR via `PROSAIL_Hybrid_Train` (`nbEnsemble = nbSamples/100` modèles), (7) sérialise via `liquidSVM::serialize.liquidSVM` et sauvegarde en `.rds`. |

##### `function_analysis.R`

Fichier de fonctions d'analyse et de visualisation (26 fonctions) :

| Fonction | Rôle résumé |
|----------|-------------|
| `create_scatterplot` | Scatter plot LAI_ALS vs LAI_S2 avec annotations R/R²/NRMSE/Bias. Échantillonnage 5 000 pts. |
| `create_density_scatterplot` | Version hexbin (`geom_hex`) du scatter plot. |
| `generate_optimal_plot` | Génère un scatter plot optimal en lisant CSV LAI estimé + CSV LiDAR LAI. |
| `find_optimal_depth` | Trouve la profondeur d_opt maximisant R² moyen, par méthode de normalisation. |
| `ensure_dir` | Crée un répertoire si inexistant. |
| `extract_levels` | Parse les niveaux numériques d'un paramètre depuis des noms de colonnes (regex). |
| `compute_metrics` | Calcule R, RMSE, NRMSE, bias, slope. **Note** : référence `final_df` dans l'environnement global — couplage implicite. |
| `fill_last_seen` | LOCF (Last Observation Carried Forward) ligne par ligne. |
| `plot_config` | Boxplot+jitter par configuration (coordonnées retournées). |
| `plot_depth` | Violin+boxplot par profondeur. |
| `normalize` | Normalisation min-max. |
| `get_param_group` | Extrait le préfixe alphabétique d'un nom de paramètre. |
| `generate_param_colors` | Palette RColorBrewer par groupe de paramètres. |
| `plot_param_freq` | Bar chart de fréquence des paramètres. |
| `plot_score_ranges` | Boxplots de distributions de scores par site. |
| `compute_param_freq` | Fréquence data.table des paramètres dans les noms de colonnes. |
| `compute_param_freq_by_score` | Fréquences filtrées par critère de score. |
| `compute_param_freq_by_score2` | Idem + colonnes Site et Score. |
| `extract_parameters` | Extrait ALA, lai, LMA, BROWN depuis noms de colonnes via regex. |
| `rank_parameters` | Rang composite = frank(-R) + frank(RMSE) + frank(\|Slope-1\|) + frank(\|Bias\|). |
| `plot_metric_density` | Densité normalisée par métrique avec lignes-seuil quantile 20 %. |
| `plot_metric_density2` | Version améliorée de `plot_metric_density` (légende unique). |
| `manual_param_config` | Retourne une tribble ALA×lai×LMA×BROWN (configs numérotées 1-5). |
| `manual_param_config_final` | Idem avec labels ("ATBD", "OPT#1"…). |
| `plot_lai_scatter` | Scatter plot final (`geom_bin2d`) avec statistiques via `broom`/`Metrics`. Gère l'affichage conditionnel des axes selon site et label. |

##### `functions.R`

Fichier **vide de fonctions**. Ne contient qu'un commentaire décrivant les fonctions de sampling et de parallélisation à implémenter (`## Sampling function... ## Parallellization function ????`). Probablement un placeholder non finalisé. **Code mort — candidat à la suppression.**

### 3.5 Scripts à la racine de `02_CODES/`

*Scripts "orphelins" au niveau racine, souvent des scripts d'exploration ou
de finalisation. À inventorier :*
- `3b_apply_masks_to_metrics.R`
- `main_explore_init.R`
- `main_stands.R`
- `plot_initial_lai_rasters_final.R`
- `plot_initial_scatterplots.R`
- `three_factors_analysis_final.R` / `_gam.R` / `_linear*.R` / `_new*.R`
- autres `plot_*.R`

> Beaucoup de variantes avec suffixes (`_final`, `_new`, `_v2`...). Identifier
> la version canonique utilisée pour les figures finales.

---

## 4. Graphe de dépendances

*Dessiner un graphe texte (arbre ou liste adjacente) qui montre l'ordre
d'exécution du pipeline pour reproduire les figures du papier. Format suggéré :*

```
[données brutes]
    |
01_DATA/{site}/LiDAR/*.las           01_DATA/{site}/Sentinel-2/*.SAFE
    |                                         |
LiDAR/0_normalize_heights.R          S2/Main_01_download_S2_sites.R
    |                                         |
LiDAR/1_create_chm.R                 S2/Main_02_produceLAI.R
    |                                         |
LiDAR/3_arrange_ladstack.R                    |
    |                                         |
    +---------------------+-------------------+
                          |
              [jointure sur grille S2]
                          |
         [analyses : d_opt, PROSAIL, hétérogénéité]
                          |
                      [figures]
```

*À compléter avec les vrais noms de scripts et les vraies dépendances.*

---

## 5. Redondances et code mort identifiés

Liste des paires ou groupes de scripts qui semblent faire la même chose :

| Fichier A | Fichier B | Différence apparente | Lequel garder ? |
|-----------|-----------|----------------------|-----------------|
| `3_arrange_ladstack.R` | `3b_arrange_ladstack.R` | A : 3 sites, ladstack.tif ; B : Aigoual only, ladstack_0408.tif, sous-dossier 3107 | ? (Q5) |
| `3_calculate_lidar_metrics.R` | `0tmp3_calculate_lidar_metrics.R` | A : pas de lidRmetrics ; B : lidRmetrics actif | A (probablement) |
| `6b_explain_heterogeneity_mix_sites.R` | `6b2_explain_heterogeneity_mix_sites.R` | A : test_size=50000, sans pls ; B : test_size=3000, avec pls, combinaisons différentes | ? (Q6) |
| `functions_chm.R` | `functions_normalize_heights.R` | A : version actuelle (6 args) ; B : vieille copie (3 args), header incorrect | A — supprimer B |
| `functions_lidar.R` → `myPAI(z, zmin, k)` | `functions_lidar2.R` → `myPAI(z, zmin=1, k)` | même nom, z0 différent (2 vs 1 m), default sur zmin | ? (Q4) |
| `1_create_chm.R` | `1_create_chm (copy).R` | A : nouvelle signature create_chm (6 args) ; B : ancienne (4 args) | A — supprimer B |
| `0_normalize_heights.R` | `0tmp_normalize_heights.R` | A : version active ; B : stop() inline, chemin disque ancien | A — supprimer B |

Liste des fichiers qui semblent être du code mort (non sourcés, non appelés) :

| Fichier | Dernière modification | Raison soupçonnée d'abandon |
|---------|-----------------------|------------------------------|
| `functions_normalize_heights.R` | 2024-05-24 | Copie mal renommée de l'ancienne functions_chm.R, aucun script ne la source |
| `0tmp_normalize_heights.R` | 2024-03-15 | stop() inline, remplacé par 0_normalize_heights.R |
| `1_create_chm (copy).R` | 2024-05-24 | Signature obsolète, remplacé par 1_create_chm.R |
| `grid_metrics_final.Rmd` | 2022-08-22 | Auteur externe (Marianne Laslier), rgdal déprécié, chemin Windows |
| `5_correlation_per_PAI_profiles.R` | 2024-03-18 | Chemins source() incorrects, format ENVI abandonné |
| `3_shadows_analysis.R` | 2023-05-28 | Exploration précoce, rayshader, pas dans le pipeline final |
| `PROSAIL-Optimization/02_CODES/libraries/functions.R` | — | Placeholder vide : uniquement des commentaires de fonctions à implémenter, aucun code R |

---

## 6. Dépendances externes et versions

### Packages R utilisés (union sur tous les scripts actifs)

| Package | Utilisé dans | Version installée | Version critique ? |
|---------|--------------|-------------------|---------------------|
| lidR | | | oui (API change régulièrement) |
| terra | | | oui |
| prosail | | | oui (version Féret à figer) |
| mgcv | | | |
| data.table | | | |
| sf | | | |
| ggplot2 | | | |

### Données externes référencées

| Source | Chemin | Utilisé pour | Accessible ? |
|--------|--------|--------------|--------------|
| ALS Aigoual | `01_DATA/Aigoual/LiDAR/` | LAI_ALS, CHM, DTM | |
| ALS Blois | `01_DATA/Blois/LiDAR/` | | |
| ALS Mormal | `01_DATA/Mormal/LiDAR/` | | |
| S2 L2A | via Copernicus | réflectances | credentials requis |
| BDForêt V2 | `01_DATA/.../Geo_Files/` | masque feuillu | |

---

## 7. Fragilités et dettes techniques repérées

Liste des problèmes techniques à corriger pendant le refactoring :

- [ ] Chemins absolus hardcodés (lister les fichiers concernés)
- [ ] Usage de `setwd()` (lister)
- [ ] Mélange Lambert-93 et UTM 31N sans reprojection explicite
- [ ] Fonctions dupliquées entre `libraries/` et scripts principaux
- [ ] Noms de variables ambigus (ex : `lai` sans préciser ALS ou S2)
- [ ] Paramètres hardcodés qui devraient être des arguments (k = 0.5,
      h_min = 2, fCover = 0.9, dmin ∈ {10, 15, 20})
- [ ] Absence de vérification d'intégrité des fichiers intermédiaires
- [ ] Absence de logs ou de reporting sur les runs longs
- [ ] Autres :

---

## 8. Questions ouvertes pour Nathan

Liste des points que Claude Code ne peut pas trancher seul et sur lesquels
il attend une réponse avant de refactoriser :

1. **Scripts actifs sur les sites du papier** : `1_create_chm.R` et `3_calculate_lidar_metrics.R` sont actuellement configurés pour `sites = c("Hayes", "Reine")`. Les résultats dans `03_RESULTS/Aigoual|Blois|Mormal/` ont été produits par une version antérieure. Quelle version a produit les résultats finaux du papier ?

2. **`3_arrange_ladstack.R` vs `3b_arrange_ladstack.R`** : que signifient `_0408` et le sous-dossier `3107` ? Retraitement partiel ou expérimentation Aigoual ?

3. **Packages `lidRmetrics`, `solaR`, `purrr`, `randomForest`, `doParallel`** : hors liste CLAUDE.md. À conserver ou remplacer lors du refactoring ?

4. ~~**Convention z0**~~ → **RÉPONDU** : z0=2 m est la référence (Bouvier). Les reviewers demandent une sensibilité à h_min ∈ {2, 3, 5 m} et fCover ∈ {80, 90, 95 %}. Tout paramètre h_min / fCover doit être un argument de fonction. `functions_lidar2.R` (z0=1 m) n'est PAS la référence pour le papier.

5. **`3b_arrange_ladstack.R` : suffixe `_0408` et dossier `3107`** — cf. Q2.

6. ~~**`6b` vs `6b2`**~~ → **RÉPONDU** : `6b2` est la version active (plus propre). Mais l'approche RF est **abandonnée pour la révision** — les scripts `6a/b/c` ne sont plus centraux. L'analyse facteurs explicatifs passe par `PROSAIL-Optimization/02_CODES/Main_*.R`.

7. **`myPAD_dtm` vs `myPAD` dans le papier** : laquelle a été utilisée pour les LAD profiles des figures 3 et A.2 ?

8. **`2.calculate_25m_metrics.R`** (MAJ 2025-11-13) : à quelle analyse de révision correspond-il ?

9. **`Depth` dans `Main_02` vs `Main_02b`** : `Main_02` utilise `Depth=5` pour Aigoual/Blois et `Depth=9` pour Mormal ; `Main_02b` utilise `Depth=4` pour tous les sites. Ces valeurs sont-elles des indices dans le vecteur de profondeurs testées, ou des hauteurs en mètres ? La convention "Depth" correspond-elle à l'intégration depuis la base (h_min = Depth m) ou depuis le sommet ?

10. **`Main_03_student.R` vs `Main_03_student2.R`** : lequel a produit les figures finales du papier soumis ? `Main_03_student2` est plus propre mais les deux ont un `stop()` à mi-chemin — le code avant ou après le `stop()` est-il la partie active ?

11. **Packages hors CLAUDE.md dans PROSAIL** : `rPref`, `randomForest`, `pdp`, `broom`, `Metrics`, `ggpubr`, `liquidSVM`, `truncnorm`, `sgsR`, `RColorBrewer`, `forcats` — à conserver ou remplacer lors du refactoring ? `liquidSVM` en particulier est un package inhabituel (SVM vectorisé) — version à figer impérativement.

12. **Score composite PROSAIL** : `Score = (norm_R + norm_NRMSE + norm_slope)/3` — les trois critères ont le même poids. Le reviewer demande une sélection Pareto multi-critères : faut-il remplacer ce score par un front de Pareto, ou l'ajouter en complément ?

---

> **Note reviewer (2026-04-09)** : Les reviewers demandent explicitement de justifier h_min=2 m et fCover>90 %, et de fournir une analyse de sensibilité : h_min ∈ {2, 3, 5 m}, fCover ∈ {80, 90, 95 %}. Ces paramètres doivent être rendus entièrement paramétrables dans le refactoring du module LiDAR (priorité haute).

---

## 9. Proposition d'ordre de refactoring

*À remplir après l'inventaire complet. Format suggéré :*

### Phase 1 — Fondations (module LiDAR)
- `revision/R/lidar_normalize.R` : fonctions de normalisation DTM/CHM
- `revision/R/lidar_chm.R` : génération et smoothing du CHM
- `revision/R/lidar_lad.R` : calcul LAD profiles (Bouvier 2015) paramétré en k
- `revision/R/lidar_lai.R` : intégration en LAI_ALS avec h_min paramétrable
- Tests unitaires : reproduction de LAI_ALS sur un pixel témoin

### Phase 2 — Module Sentinel-2
- `revision/R/s2_preprocess.R`
- `revision/R/s2_masks.R` : masques feuillu et fCover paramétrable
- Tests : reproduction d'une tuile

### Phase 3 — Module PROSAIL
- `revision/R/prosail_lut.R` : génération LUT paramétrée
- `revision/R/prosail_invert.R` : wrapper inversion hybride
- `revision/R/prosail_metrics.R` : r, RMSE, bias, slope, Pareto

### Phase 4 — Analyses de révision
- `revision/R/analysis_dopt.R` : sélection d_opt Pearson + multi-critères
- `revision/R/analysis_heterogeneity.R` : CHM_std + DSM_std
- `revision/R/analysis_sensitivity.R` : k, h_min, fCover, Mormal-phase1

### Phase 5 — Orchestration et reproduction
- `revision/scripts/01_...` à `revision/scripts/NN_...`
- Vignette Quarto sur mini-dataset
- README encadrants

---

## 10. Log de la phase d'inventaire

Journal datée des actions de lecture menées par Claude Code pendant cette
phase, pour traçabilité.

| Date | Action | Fichiers lus | Observations |
|------|--------|--------------|--------------|
| 2026-04-08 | Listage arborescence | `02_CODES/LiDAR/` (hors archive/) | 28 fichiers R/Rmd identifiés + images/rasters de test |
| 2026-04-08 | Lecture scripts preprocessing | `0.change_filenames_L93_UTM31N.R`, `0_convert_l93_into_utm.R`, `0_convert_l93_into_utm_mask.R`, `0_normalize_heights.R` (partial), `0tmp_normalize_heights.R` | Deux variantes normalisation, chemin disque différent (espace vs pas d'espace) |
| 2026-04-08 | Lecture scripts CHM et masques | `1_create_chm.R`, `1_create_chm (copy).R`, `2_create_vegetation_forest_masks.R` | 1_create_chm actif sur Hayes/Reine (hors papier) |
| 2026-04-08 | Lecture scripts métriques | `3_calculate_lidar_metrics.R` (partial), `3.calculate_lidar_metrics2.R`, `0tmp3_calculate_lidar_metrics.R` (partial) | 3 variantes de calcul métriques ; actif sur Hayes/Reine |
| 2026-04-08 | Lecture scripts LADstack | `3_arrange_ladstack.R`, `3b_arrange_ladstack.R` | Headers identiques mais contenu différent ; suffixe _0408 non élucidé |
| 2026-04-08 | Lecture scripts analyse | `4_prepare_heterogeneity_depth_analysis.R`, `5_heterogeneity_depth_analysis.R`, `5_correlation_per_PAI_profiles.R` (partial) | 5_correlation cassé (source() incorrects) ; 5_heterogeneity produit figures papier |
| 2026-04-08 | Lecture scripts RF | `6a_explain_heterogeneity_unique_site.R`, `6b_explain_heterogeneity_mix_sites.R`, `6b2_explain_heterogeneity_mix_sites.R`, `6c_explain_heterogeneity_stability.R` | 6b vs 6b2 : headers identiques mais paramètres différents |
| 2026-04-08 | Lecture scripts viz/exploration | `0_good_parcel_viz.R` (partial), `0_viz_las.R`, `main_viz_norms.R`, `max_hists.R`, `pad_boxplots.R`, `2.calculate_25m_metrics.R`, `3_shadows_analysis.R`, `grid_metrics_final.Rmd` | max_hists et pad_boxplots récents (2025) |
| 2026-04-08 | Lecture libraries | `functions_lidar.R` (partial), `functions_lidar2.R`, `functions_normalize_heights.R`, `functions_chm.R`, `functions_create_masks.R` (partial) | functions_normalize_heights.R = copie obsolète de functions_chm.R ; debug prints dans myPAD_dsm |
| 2026-04-08 | Lecture scripts PROSAIL Main_* | `Main_03_student.R` (complet), `Main_03_student2.R` (complet), `Main_04_final_factors_tabs.R` (complet), `Main_01_analysis_optim_depth.R` (120 lignes), `Main_02_analysis_best_prosail_config.R` (120 lignes), `Main_02b_final_analysis_best_prosail_config.R` (120 lignes) | stop() à mi-chemin dans Main_03/04 ; deux stratégies LUT identifiées ("LIDFa_..." vs "LIDFa_...Agg_10m") |
| 2026-04-08 | Lecture libraries PROSAIL | `define_parm_combinations.R`, `get_combination.R`, `get_s2_angles.R`, `get_s2_samples.R`, `PROSAIL_train_sensitivity.R`, `PROSAIL_train_sensitivity_subset.R`, `function_analysis.R` (complets), `functions.R` (vide) | functions.R = placeholder vide ; liquidSVM hors CLAUDE.md ; compute_metrics référence global final_df |
