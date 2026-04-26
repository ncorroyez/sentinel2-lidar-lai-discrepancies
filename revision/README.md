# Révision — LiDAR × Sentinel-2 LAI (Corroyez et al., RSE)

Ce dossier contient tout le code refactorisé et les analyses additionnelles
produites pour la révision majeure du manuscrit RSE-D-25-04417.

## Structure

```
revision/
  R/              fonctions (une thématique par fichier)
  scripts/        orchestration numérotée (01… → 20…)
  output/
    figures/      figures produites (sm5/, sm6/, scan_angle/, …)
    tables/       tableaux CSV
    intermediate/ résultats intermédiaires (RDS, CSV lourds — gitignorés)
  reviewers/      commentaires reviewers, notes de réponse
  tests/          tests unitaires
```

## Configuration des chemins (collaborateurs)

Les données brutes et résultats externes sont sur le partage INRAE :
`smb://pnas3.stockage.inrae.fr/mo-mtd-pulse/root/_PROJETS/2023_2026_These_Nathan_Corroyez/`

### 1. Monter le partage SMB

**Linux (GNOME/Nautilus)** : ouvrir Nautilus → Connexion à un serveur →
`smb://pnas3.stockage.inrae.fr/mo-mtd-pulse` → identifiants LDAP INRAE.
Le point de montage sera de la forme :
`/run/user/<uid>/gvfs/smb-share:server=pnas3.stockage.inrae.fr,share=mo-mtd-pulse/`

**macOS** : Finder → Aller → Se connecter au serveur →
`smb://pnas3.stockage.inrae.fr/mo-mtd-pulse` → identifiants LDAP INRAE.
Le partage apparaît dans `/Volumes/mo-mtd-pulse/`.

### 2. Créer `config.yml` à la racine du projet

Copier le modèle et ajuster le chemin :

```bash
cp config.yml.template config.yml
```

Éditer `config.yml` pour pointer vers le dossier projet sur le partage monté,
par exemple :

```yaml
# Linux
data_root: /run/user/1234/gvfs/smb-share:server=pnas3.stockage.inrae.fr,share=mo-mtd-pulse/root/_PROJETS/2023_2026_These_Nathan_Corroyez/NC_Full

# macOS
data_root: /Volumes/mo-mtd-pulse/root/_PROJETS/2023_2026_These_Nathan_Corroyez/NC_Full
```

Le fichier `config.yml` est dans `.gitignore` — ne jamais le committer.

### 3. Vérifier les chemins

En R, depuis la racine du projet :

```r
source("revision/R/paths.R")
str(paths)
# doit afficher les chemins résolus vers 01_DATA/, 03_RESULTS/, etc.
file.exists(paths$raw_data)   # TRUE si le partage est monté
```

**Fallback** : si `config.yml` est absent, `paths.R` utilise `here::here()`
(racine du projet local). Nathan peut donc lancer les scripts sans `config.yml`
si les données sont présentes localement.

## Ordre d'exécution du pipeline

Les scripts sont numérotés dans l'ordre de dépendance :

| Script | Rôle |
|--------|------|
| `01a_download_s2_from_safe.R` | Extraire réflectances S2 depuis archives SAFE |
| `01b_download_s2_from_cdse.R` | Télécharger S2 depuis CDSE (alternative) |
| `02a_train_prosail_atbd.R` | Entraîner SVR ATBD (1 config × 3 sites) |
| `02_dryrun_train_prosail.R` | Dry run — valider l'API et estimer les temps |
| `03a_sample_s2_pixels.R` | Échantillonner les pixels S2 (stratified uniform) |
| `03b_extract_lidar_at_samples.R` | Extraire PAD LiDAR aux positions d'échantillons |
| `04a_apply_prosail_atbd.R` | Appliquer le SVR ATBD aux échantillons |
| `05a_sm5_compute_metrics_atbd.R` | Calculer les métriques ATBD par profondeur |
| `05b_sm5_plot_dopt_metrics.R` | Figures métriques vs profondeur (SM5) |
| `05c_sm5_k_zmin_sensitivity_dopt.R` | Sensibilité de d_opt à k et z_min |
| `06_sm5_select_dopt.R` | Sélectionner d_opt par critère Pareto |
| `07_compute_lai_als_dopt.R` | Calculer LAI_ALS_dopt (rasters) — **reprend si le .tif existe** |
| `08_sm5_scatter_lai_atbd.R` | Figure scatter LAI_ALS vs LAI_S2_ATBD |
| `09_train_prosail_full.R` | Entraîner les 270 SVR (per_site / common) — **reprend au scénario/site manquant** |
| `10_apply_prosail_full.R` | Appliquer les 270 SVR aux pixels S2 — **reprend si les CSV existent** |
| `11_sm5_compute_metrics_full.R` | Métriques pour toutes les configs — **reprend les chunks manquants** |
| `12_sm5_select_prosail_opt.R` | Sélectionner la config PROSAIL optimale |
| `12b_sm5_plot_prosail_pareto.R` | Figure Pareto R² vs RMSE (config PROSAIL) |
| `13_sm5_predict_lai_raster.R` | Prédire LAI_S2 en mode raster — 3 scénarios : `per_site`, `common`, `fixed_4` |
| `14_sm6_compute_heterogeneity.R` | Calculer les rasters d'hétérogénéité (SM6) |
| `15_sm6_analysis.R` | Analyser LAI vs hétérogénéité (SM6) |
| `16_sm6_analysis_sweep.R` | Sweep multi-échelle SM6 |
| `17_k_sensitivity.R` | Sensibilité LAI_ALS au coefficient k |
| `18_fcover_sensitivity.R` | Sensibilité au seuil fCover |
| `18_scan_angle_correction.R` | Correction angle de visée (complet) |
| `18a_scan_angle_tiles.R` | Correction par tuile LAS |
| `18b_scan_angle_finalize.R` | Finaliser la correction (mosaïque) |
| `19_h_min_sensitivity.R` | Sensibilité au seuil h_min |
| `20_sm6_plot_heterogeneity.R` | Figures hétérogénéité (SM6) |

Le script `run_pipeline.sh` à la racine du projet enchaîne les étapes
principales dans l'ordre correct.

### Reprise sur interruption (idempotence partielle)

Les scripts 07, 09, 10 et 11 détectent automatiquement les sorties déjà
produites et les sautent :

- **07** — skip par `(site × scénario)` si le `.tif` LAI_ALS_dopt existe.
- **09** — skip scénario entier si tous les 270 RDS existent pour les 3 sites ;
  sinon skip par site si tous les RDS du site sont présents (le pool cross-site
  est toujours recalculé pour les sites restants).
- **10** — skip par `(site × scénario × h_min)` si les deux CSV (mean + SD) existent.
- **11** — charge le CSV existant au démarrage, identifie les chunks
  `(scénario × h_min)` déjà calculés, et ne recalcule que les manquants avant
  de réécrire le fichier fusionné.

En cas d'interruption, il suffit de relancer le script ou `run_pipeline.sh` —
les étapes terminées sont ignorées.

## Données nécessaires

Les entrées (lecture seule) sont résolues via `paths.R` à partir de `config.yml` :

| Clé `paths` | Contenu |
|-------------|---------|
| `paths$raw_data` | `01_DATA/` — données brutes (S2 SAFE, géoréférences, masques) |
| `paths$ext_results` | `03_RESULTS/` — réflectances S2 prétraitées, PAD profiles, masques |
| `paths$prosail_lidar` | `PROSAIL-Optimization/01_DATA/` — LAD stacks, rasters PAD |
| `paths$prosail_codes` | `PROSAIL-Optimization/02_CODES/` — fonctions PROSAIL (`get_s2_angles`, etc.) |

Les sorties vont dans `revision/output/` (figures, tables, intermédiaires).

## Packages R requis

```r
install.packages(c(
  "here", "terra", "sf", "lidR",
  "prosail", "prospect",
  "data.table", "readr",
  "ggplot2", "patchwork", "scales",
  "cli", "yaml", "sgsR", "dplyr",
  "e1071"   # SVR via prosail
))
```
