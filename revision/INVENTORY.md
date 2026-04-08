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

*À remplir. Scripts prioritaires identifiés à partir de l'arborescence :*
- `0_normalize_heights.R`
- `1_create_chm.R`
- `2_create_vegetation_forest_masks.R`
- `3_arrange_ladstack.R` / `3b_arrange_ladstack.R`
- `3_calculate_lidar_metrics.R` / `3.calculate_lidar_metrics2.R`
- `4_prepare_heterogeneity_depth_analysis.R`
- `5_heterogeneity_depth_analysis.R`
- `5_correlation_per_PAI_profiles.R`
- `6a_explain_heterogeneity_unique_site.R`
- `6b_explain_heterogeneity_mix_sites.R` / `6b2_*`
- `6c_explain_heterogeneity_stability.R`

> Note : plusieurs scripts apparaissent en doublons numérotés (ex :
> `3_calculate_lidar_metrics.R` et `3.calculate_lidar_metrics2.R`). Identifier
> lequel est actif et lequel est une version antérieure à archiver.

### 3.2 Scripts Sentinel-2 — `02_CODES/Sentinel_2/`

*À remplir. Scripts prioritaires :*
- `Main_01_download_S2_sites.R`
- `Main_02_produceLAI.R` / `Main_02_produceLAI_v2.R`
- `0_verif_atbd_biophysical_toolbox.R` / `0.verif_atbd_biophysical_toolbox_v2.R`
- `1_mask_s2_ts.R`
- `3_train_predict_prosail.R` / `3.train_predict_prosail_ts.R`

> Note : idem, identifier les versions actives vs legacy.

### 3.3 Libraries — `02_CODES/libraries/`

*Fonctions utilitaires. À inventorier fonction par fonction pour identifier
celles qui sont effectivement utilisées dans le pipeline actif vs celles qui
sont du code mort.*

Fichiers à inspecter :
- `functions_chm.R`
- `functions_lidar.R` / `functions_lidar2.R`
- `functions_normalize_heights.R`
- `functions_sentinel_2.R`
- `functions_create_masks.R`
- `functions_heterogeneity_depth_analysis.R`
- `functions_explain_heterogeneity.R` / `functions_explain_heterogeneity_improved.R`
- `functions_intercomparisons.R`
- `functions_plots.R`
- `functions_general_tools.R`
- `functions_JBF.R` (fonctions de Jean-Baptiste Féret ?)
- `functions_convert_l93_into_utm.R`
- `functions_shadows_analysis.R`

Pour chaque fichier de `libraries/`, produire une liste :

```
functions_xxx.R :
  - nom_fonction_1(args) : rôle, utilisée par [scripts]
  - nom_fonction_2(args) : rôle, utilisée par [scripts]
  - ...
```

### 3.4 Scripts PROSAIL — `PROSAIL-Optimization/02_CODES/`

*À remplir. Scripts prioritaires :*
- `Main_01_analysis_optim_depth.R`
- `Main_02_analysis_best_prosail_config.R`
- `Main_02b_final_analysis_best_prosail_config.R`
- `Main_03_student.R` / `Main_03_student2.R`
- `Main_04_final_factors_tabs.R`

> Ce sous-projet a sa propre arborescence `01_DATA/`, `02_CODES/`, `03_RESULTS/`.
> Documenter comment il s'articule avec le pipeline principal.

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
| | | | |

Liste des fichiers qui semblent être du code mort (non sourcés, non appelés) :

| Fichier | Dernière modification | Raison soupçonnée d'abandon |
|---------|-----------------------|------------------------------|
| | | |

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

1.
2.
3.

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
| | | | |
