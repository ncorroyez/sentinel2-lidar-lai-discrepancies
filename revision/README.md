# Révision — LiDAR × Sentinel-2 LAI (Corroyez et al., RSE)

Code refactorisé et analyses additionnelles pour la révision majeure du
manuscrit RSE-D-25-04417.

---

## Lancer le pipeline

```r
source("revision/scripts/00_run_all.R")
```

---

## Reproduire le pipeline sans accès au serveur INRAE

1. Cloner le code depuis GitHub
2. Obtenir le bundle de données (`revision/data.zip`) et l'extraire dans `revision/data/`
3. Copier `config.yml.template` en `config.yml` et choisir `profile: local`
4. Ajuster `data_root` et `output_root` dans le bloc `local` de `config.yml`
5. `source("revision/scripts/00_run_all.R")` (~3–6 h)

Le bundle est généré depuis le serveur INRAE avec :

```r
source("revision/bundle_inputs.R")  # nécessite SMB monté
```

---

## Structure du dossier

```
revision/
├── R/                      fonctions R (une thématique par fichier)
├── scripts/
│   ├── 00_run_all.R        orchestration maître — source les 9 phases
│   ├── 01_prepare.R        Phase 1 : SVR ATBD + échantillonnage S2 + PAD
│   ├── 02_dopt_compute.R   Phase 2 : métriques ATBD + sélection d_opt + rasters
│   ├── 02_dopt_figures.R   Phase 2 : figures r vs profondeur + scatter
│   ├── 03_prosail_compute.R Phase 3 : 270 SVR + LAI_S2_opt + rasters ATBD
│   ├── 03_prosail_figures.R Phase 3 : figures Pareto + param + validation ATBD
│   ├── 04_het_compute.R    Phase 4 : DSM_sd/CHM_sd + régression SM6
│   ├── 04_het_figures.R    Phase 4 : figures hétérogénéité + profils PAD
│   ├── 05_sensitivity.R    Phase 5 : sensibilités k, fCover, h_min, joint
│   ├── 05_sensitivity_figures.R  Phase 5 : figures de sensibilité
│   ├── steps/              scripts constitutifs (détail d'implémentation)
│   └── archive/            scripts non-pipeline (téléchargements S2, scan angle)
├── output/
│   ├── figures/            figures produites (sm5/, sm6/, reviewers/, …)
│   ├── tables/             tableaux CSV finaux
│   └── intermediate/       résultats intermédiaires (gitignorés)
├── reviewers/              commentaires reviewers + notes de réponse
└── tests/                  tests unitaires
```

---

## Configuration initiale

### 1. Créer `config.yml`

```bash
cp config.yml.template config.yml
```

Choisir le profil dans `config.yml` :

```yaml
profile: local   # données dans revision/data/  (défaut)
# profile: smb   # données sur le serveur INRAE partagé
```

Le fichier `config.yml` est dans `.gitignore` — ne jamais le committer.

### 2. Profil `local` — données en local

Les données sont dans `revision/data/` (non versionnées, lourdes).
Peupler ce dossier via `bash sync_to_smb.sh --outputs` depuis le serveur,
ou copier manuellement depuis `03_RESULTS/`, `01_DATA/`, `PROSAIL-Optimization/01_DATA/`.

### 3. Profil `smb` — serveur INRAE partagé

**Linux (GNOME/Nautilus)** : Connexion à un serveur →
`smb://pnas3.stockage.inrae.fr/mo-mtd-pulse` → identifiants LDAP INRAE.
Point de montage : `/run/user/<uid>/gvfs/smb-share:server=pnas3.stockage.inrae.fr,share=mo-mtd-pulse/`

**macOS (Finder)** : Aller → Se connecter au serveur →
`smb://pnas3.stockage.inrae.fr/mo-mtd-pulse`.
Point de montage : `/Volumes/mo-mtd-pulse/`

Adapter l'UID dans `config.yml` si nécessaire (Linux : `id -u`).

### 4. Vérifier les chemins

```r
source("revision/R/paths.R")
str(paths)
file.exists(paths$raw_data)    # TRUE si données accessibles
```

---

## Synchronisation SMB ↔ local

```bash
bash sync_to_smb.sh            # tout (entrées + sorties)
bash sync_to_smb.sh --inputs   # données d'entrée seulement (~3 Go)
bash sync_to_smb.sh --outputs  # sorties revision/output/ seulement
```

Le script ne synchronise que les sous-dossiers nécessaires au pipeline
(pas les 113 Go complets).

---

## Phases du pipeline

| Script | Durée approx. | Rôle |
|--------|--------------|------|
| `01_prepare.R` | ~30 min | SVR ATBD, échantillonnage S2, extraction PAD LiDAR |
| `02_dopt_compute.R` | ~5 min | Métriques ATBD × profondeur, sélection d_opt Pareto, rasters LAI_ALS_dopt |
| `02_dopt_figures.R` | ~2 min | Figures r vs profondeur (SM5), scatter LAI_ALS vs LAI_S2 |
| `03_prosail_compute.R` | ~3–5 h | 270 SVR PROSAIL × 2 scénarios, métriques, sélection Pareto, rasters LAI_S2_opt + ATBD_T/F |
| `03_prosail_figures.R` | ~3 min | Pareto front, distributions paramètres, validation ATBD |
| `04_het_compute.R` | ~10 min | DSM_sd, CHM_sd, régression SM6b, sweep configurations |
| `04_het_figures.R` | ~5 min | Figures hétérogénéité, scatter densité, histogrammes S2, profils PAD |
| `05_sensitivity.R` | ~15 min | Sensibilités k, fCover ∈ {80,90,95 %}, h_min ∈ {2,3,4,5} m, joint |
| `05_sensitivity_figures.R` | ~3 min | Figures de sensibilité (réponse reviewers R3, R4) |

Exclusions de `run_all()` (lancer manuellement si nécessaire) :
- `archive/01a_download_s2_from_safe.R` — extraction réflectances S2 depuis SAFE
- `archive/01b_download_s2_from_cdse.R` — téléchargement S2 via CDSE
- `archive/18_scan_angle_*` — correction angle de visée (nécessite les LAS sur SMB)

---

## Chemins exposés par `paths.R`

| Variable | Profil `local` | Profil `smb` |
|----------|---------------|-------------|
| `paths$raw_data` | `revision/data/raw/` | `01_DATA/` |
| `paths$ext_results` | `revision/data/results/` | `03_RESULTS/` |
| `paths$prosail_lidar` | `revision/data/prosail_lidar/` | `PROSAIL-Optimization/01_DATA/` |
| `paths$prosail_codes` | `PROSAIL-Optimization/02_CODES/` | `PROSAIL-Optimization/02_CODES/` |
| `paths$output` | `revision/output/` | `revision/output/` |

---

## Packages R requis

```r
install.packages(c(
  "here", "terra", "sf", "sgsR", "dplyr",
  "prosail", "prospect",
  "data.table", "readr",
  "ggplot2", "patchwork", "scales",
  "cli", "yaml",
  "e1071"
))
```

---

## Qualité des figures

Toutes les figures sont sauvegardées en **PDF** (`device = cairo_pdf`) et
**PNG 600 dpi** (`bg = "white"`).
