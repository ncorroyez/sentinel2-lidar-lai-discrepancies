# HARMONISATION des deux dépôts de thèse — document unique (MAJ 2026-09-14)

**Source unique = cette copie (`NC_Full/HARMONISATION.md`).** La copie dans rmusica est un miroir.
Remplace `HARMONISATION_carto_2026-09-14.md` et `ARBORESCENCE_CIBLE_2026-09-14.md` (les deux, dans les
deux dépôts), archivés dans `_archive_harmonisation_2026-09-14/` de chaque dépôt.

**Dans ce dépôt (NC_Full)** : Ch2 (RSE), Ch3, Ch4/GEDI. Code dans `scripts/{ch2,ch3,ch4}/` + `R/`,
manuscrits dans `manuscripts/{ch2_rse,ch3,ch4}/`, plans FR dans `plans/`, bib canonique
`Bib_nathan_corroyez_10Sep.bib`. L'autre dépôt (Ch1, intro, discussion, moteur MuSICA) est
`~/Documents/z_Example_rmusica_31012025` (« rmusica »).

<!-- corps commun : identique dans les deux copies -->

Conventions : `NC_Full/` = `~/Documents/NC_Full`, `rmusica/` = `~/Documents/z_Example_rmusica_31012025`.
Tous les chemins ci-dessous sont préfixés par le dépôt, les deux copies sont donc justes.
Zones read-only NC_Full jamais touchées : `01_DATA 02_CODES 03_RESULTS 04_FIGURES paper1_Figures PROSAIL-Optimization`.

## 1. ÉTAT FINAL RÉEL (vérifié sur disque le 2026-09-14 soir)

| Chapitre | Manuscrit CANONIQUE | Code vivant | Figures/tables lues par le manuscrit, compile | Dépôt |
|---|---|---|---|---|
| **Ch1** LiDAR × MuSICA | `rmusica/chapter1/manuscript/manuscript_chap1_FINAL_coherence_2026-09-09.md` (+ `_2026-09-11.pdf/.docx`) | `rmusica/chapter1/producers/` (175 `c1_*.R`) ; orchestration `chapter1/orchestration/run_chapter1_chs41.R` (vivant : Fig 4/5/6/7/B1/S1) + `run_chapter1.R` (SUPERSEDED en partie : Fig 1/2/3/8/A1/A2/C1/D1) ; ledger `chapter1/ledger/make_article_figure_set.R` (14 figures, 0 MISSING) ; libs **à la racine** : `rmusica/R/` (78), `rmusica/pipeline/` (34), binaire `rmusica/in_files/model-3.2.3/musica` | 14 figures embarquées : `chapter1/manuscript/figures/` (7, chaîne chs41) + `figures/article_v323/` (7). Compile depuis `chapter1/manuscript/`, YAML `bibliography: ../../Bib_these.bib`, `apa.csl` local | rmusica |
| **Intro générale** | `rmusica/thesis_frame/intro_generale_these_EN.md` (+pdf), YAML `bibliography: ../Bib_these.bib` | `rmusica/thesis_frame/fig_intro_*.R` (12) + `_intro_fig_style.R` → écrivent `thesis_frame/figures/` | 16 figures dans `thesis_frame/figures/`. Plan EN `thesis_frame/plan_intro_generale_these.md` ; plan FR `NC_Full/plans/PLAN_intro_generale_FR.md` | rmusica (+ plan FR NC_Full) |
| **Discussion générale** + Annexe G | `rmusica/thesis_frame/discussion_generale_these_EN.md` (+pdf ; pas de YAML → `--bibliography ../Bib_these.bib`) ; Annexe G = `rmusica/thesis_frame/annexe_GEDI_EN.md` (+pdf) | — (Annexe G s'appuie sur le code Ch4) | Annexe G lit `thesis_frame/figs_gedi/` (4 figures). Plan FR `NC_Full/plans/PLAN_discussion_generale_FR.md` | rmusica (+ plan FR NC_Full) |
| **Ch2** RSE S2 × LiDAR | `NC_Full/manuscripts/ch2_rse/submission_Sep26/Manuscript_Corroyez_2025_RSE_r3.docx` (Word, r3, resoumission 15/09) ; mode d'emploi `manuscripts/ch2_rse/README_Sep26.md` | `NC_Full/scripts/ch2/00_run_all.R` → `01..05_*.R` → `scripts/ch2/steps/` (44) + `NC_Full/R/` ; legacy `02_CODES/` (read-only) | Sorties `NC_Full/output/` ; figures/tables soumises `manuscripts/ch2_rse/Figures_Sep26/`, `Tables_Sep26/`. `NC_Full/paper1_sub_August26/` (r2) GARDÉ = ancre du skill `nathan-scientific-voice` | NC_Full |
| **Ch3** forçage LAI → microclimat | `NC_Full/manuscripts/ch3/Chapter3_article_standalone_EN.md` (+pdf du 14/09). **PIÈGE** : `Chapter3_article_plan_EN.md` a un mtime PLUS RÉCENT mais c'est la lignée 3 sites abandonnée (voir `Chapter3_VERSION_TRIAGE.md`) | `NC_Full/scripts/ch3/` (32 : `c3_*`, `make_ch3_*`, `make_fig*`, `fig_h2_table3f`), sourcent `NC_Full/scripts/_article_style.R` et `_figure_label_check.R` ; **la partie MuSICA du Ch3 vit à la racine de rmusica** (65 `c3_*.R`, écrivent en absolu vers `NC_Full/manuscripts/ch3/`) | `manuscripts/ch3/figures/` (116) et `tables/` (153), écrits par 29 scripts de `scripts/ch3/` ; consomme les intermédiaires Ch2 (sm5/sm6) dans `NC_Full/output/`. Compile : pas de YAML → bib + csl racine NC_Full | NC_Full (+ rmusica pour MuSICA) |
| **Ch4 / GEDI** (perspective) | `NC_Full/manuscripts/ch4/Chapter4_prototype_plan_EN.md` (plan seul) ; le texte canonique est l'Annexe G (rmusica) | `NC_Full/scripts/ch4/` (8 `c4_*.R`) + racine rmusica (20 `c4_*`/`*gedi*.R`) | `manuscripts/ch4/figures/`, `tables/` ; `rmusica/thesis_frame/figs_gedi/` | les deux |
| **Thèse assemblée** | `NC_Full/thesis.md` (+pdf) — YAML `bibliography: Bib_nathan_corroyez_10Sep.bib`, `csl: apa.csl`. Copie DIVERGÉE du Ch1 (pas une inclusion) | — | lit `NC_Full/figures/article_v323` = symlink → `rmusica/chapter1/manuscript/figures/article_v323` (repointé le 14/09) | NC_Full |
| **Plans FR** | `NC_Full/plans/` : `PLAN_{ch3,intro_generale,discussion_generale}_FR.md` (+docx), `PLAN_chapitre_M_R_D.md`, `ANALYSES_Ch{1,3}_*`, `Retours_encadrants_2026-09-02_FR.md`, `PLANS_FR_voice_review.md` (14 fichiers, 0 référence script) | — | — | NC_Full |
| **Bibliographies** | NC_Full : **une seule**, `NC_Full/Bib_nathan_corroyez_10Sep.bib` (965 entrées ; les 136 clés citées Ch3 + thèse y sont). rmusica : `rmusica/Bib_these.bib` (compile Ch1/intro/disc, régénéré par `rmusica/fix_bibliography.R`) — les exports `Bib_nathan_corroyez{June13,12Aug}.bib` et `bib_fixes.bib` restent à la racine rmusica | | | une par dépôt |

Arborescences résultantes :

```
rmusica/                                   NC_Full/
  chapter1/                                  scripts/
    manuscript/   (FINAL_coherence + figures/)   ch2/ (00..05 + steps/ 44)   ch3/ (32)   ch4/ (8)
    producers/    (175 c1_*.R)                   _article_style.R, _figure_label_check.R, 3 utils
    orchestration/(run_chapter1*.R, RUN_*.md)   R/  output/  tests/  reviewers/ (interne)
    ledger/       (make_article_figure_set.R)  manuscripts/  ch2_rse/  ch3/  ch4/
    _reference/   (README, RECAP, config_musica, plans/        (14 plans FR)
                   comparaison_versions, reu_20260626) Bib_nathan_corroyez_10Sep.bib  apa.csl  thesis.md
  thesis_frame/   (intro + discussion + annexe_GEDI  figures/article_v323 -> symlink rmusica
                   + fig_intro_*.R + figures/ + figs_gedi/  _archive_harmonisation_2026-09-14/ (177 Mo)
                   + plans EN, crib, meetings, crosschecks) _archive_redaction_2026-09-14/ (12 Mo)
  R/ pipeline/ in_files/ out_files/ transfer/  # RO : 01_DATA 02_CODES 03_RESULTS 04_FIGURES
  prep_site_forcing/ scripts/ (116 make_*)     #      paper1_Figures PROSAIL-Optimization
  ~85 c3_*/c4_* MuSICA à la racine
  _archive_harmonisation_2026-09-14/ (322 Mo)
  _archive_redaction_2026-09-14/ (1 Mo)
  chapter1/manuscript/_ARCHIVE_perime_2026-09-12/ (137 Mo)
```

**Règle rmusica non négociable** : `R/`, `pipeline/`, `in_files/` (dont le binaire MuSICA) RESTENT à la
racine — 321 `source("R/…")` relatifs, `list.files("R")` dynamiques, binaire appelé par chemin.
**Règle NC_Full** : le code utilise `here::here()` depuis la racine projet → un `mv` de script ne casse
que les `source()` intra-groupe et les appelants par chemin (méthode validée : audit → `mv` → sed →
`Rscript -e parse+resolve`).

## 2. RESTE / DÉCISIONS OUVERTES (par priorité)

1. **rmusica n'a AUCUN remote git** (`git remote -v` vide) : 5 fichiers suivis, 361 entrées non
   commitées, branche `fix/lad-hmax-coherence-rescaling`, dernier commit = ère Shapley (périmée). Le Ch1,
   l'intro, la discussion et le moteur MuSICA n'existent que sur ce disque. **Sauvegarde = priorité 1** :
   remote forge INRAE + push du code/prose ; `in_files/` (23 Go), `out_files/` (305 Go), `SAFRAN/` (12 Go)
   hors git → SMB.
2. **NC_Full git.** `.gitignore` corrigé le 14/09 (lignes 194-212 : `manuscripts/ plans/ _archive_* thesis*
   *.bib *.csl outputs/ figures/ …` ignorés, `origin` = dépôt GitHub PUBLIC du pipeline Ch2). Restent :
   (a) `scripts/ch3/` et `scripts/ch4/` sont `??` (non ignorés) → un `git add -A` publierait le code Ch3/Ch4
   sur `origin` — décider public ou à ignorer ; (b) le déplacement Ch2 n'est pas commité (`D scripts/0X.R`,
   `?? scripts/ch2/`) ; (c) `HARMONISATION.md` est `??` (contient la carte interne) — l'ignorer ou l'assumer ;
   (d) un second remote `inrae` existe, non poussé depuis la réorg.
3. **Consolidation « en un endroit » — verdict : OMBRELLE, PAS FUSION.** Créer `~/Documents/These_Corroyez/`
   avec des symlinks vers les deux dépôts + `shared/` (bib + csl + style) + `build_thesis.sh`
   (`thesis.md` est aujourd'hui une copie divergée du Ch1). Pas de fusion physique : mémoires Claude keyées
   par chemin, ~342 Go + ~170 Go de données, disciplines de chemins différentes (`here::here` vs
   `source("R/")` relatif). À faire après le 15/09.
4. **`rm` des `_archive_*`** : NON fait, attend le go explicite de Nathan. 5 dossiers, ~650 Mo
   (`rmusica/_archive_harmonisation_2026-09-14/` 322 Mo, `rmusica/chapter1/manuscript/_ARCHIVE_perime_2026-09-12/`
   137 Mo, `NC_Full/_archive_harmonisation_2026-09-14/` 177 Mo, les deux `_archive_redaction_2026-09-14/`).
   Tout est réversible par `mv` inverse tant que rien n'est supprimé.
5. **Producteurs c1_* non référencés.** `chapter1/producers/_unreferenced/` N'EXISTE PAS (idée notée, non
   exécutée). État réel : 175 producers gardés, dont **18 orphelins HELD** délibérément (6 `c1_diag_*`
   demandes encadrants, `fig3_attribution_native`, `cmp_operating_slope`, tests de provenance,
   `*_CHS41` à risque) ; 67 `c1_*` archivés en 7 sous-dossiers de `rmusica/_archive_harmonisation_2026-09-14/c1_*`.
   Décision : créer `_unreferenced/` pour les 18 ou les laisser ; par lots + `run_chapter1_chs41.R --check`.
6. **`rmusica/SAFRAN/` 12 Go stale**, non suivi, hors périmètre des chaînes vivantes → sortir du dépôt (SMB).
7. **Bib rmusica non consolidée** : `Bib_these.bib` est le compilé, mais 3 exports + `bib_fixes.bib`
   (+ 2 `.bak`) restent à la racine et une copie périmée dans `chapter1/manuscript/`. Piège connu : citekey
   `hawinkelOutofsample$R^2$…` malformée (nettoyer sur COPIE, jamais l'export).
8. **Références périmées restantes dans la prose/doc (non bloquantes)** : `rmusica/thesis_frame/annexe_GEDI_EN.md`
   et 5 autres `.md` de `thesis_frame/` citent encore `NC_Full/chapter3_S2_LAI|chapter4_GEDI` ; 7 `.md` de
   `thesis_frame/` citent `review/` ou `Chapitre1/` ; le skill `these-manuscript` (§Compile) pointe encore
   `$RMUSICA/review/` ; `NC_Full/CLAUDE.md` renvoie à `ARBORESCENCE_CIBLE_2026-09-14.md` (désormais archivé →
   pointer ce fichier, modification de CLAUDE.md soumise à accord).
9. **Ledger Ch1** : le self-check `grepl("NOT LOCATED")` est vide de sens (colonne label) → ajouter un
   `file.exists()` réel (chantier code, pas doc).
10. **Nommage** `chapter1/` (rmusica) vs `ch2/ch3/ch4/` (NC_Full) à unifier — différé, 🔴 (chemins réécrits).
11. **Racine NC_Full** : ~24 `thesis_*` backups (.md/.pdf), `pipeline.RData`, 12 logs, `Retours_encadrants`
    docx/pdf en double avec `plans/` — non triés.
12. **Ch1 figures** : deux `Fig6` dans `rmusica/out_files/Chapter1/figures/` (gaussian = embarqué,
    h2_controlled = non) ; PNG B1/S1 régénérés le 14/09 par les producteurs réécrits, en attente de
    validation visuelle avant d'écraser les originaux (`_bak_handmade_2026-09-14/`).

## 3. JOURNAL condensé du 2026-09-14

- **Cartographie** read-only des deux dépôts ; identification des frictions : (A) l'orchestrateur Ch1
  documenté ne produisait pas le Ch1 vivant (chaîne CHS41 lancée à la main, ledger périmé) ; (B) intro/
  discussion éclatées sur deux dépôts ; (C) code Ch3/Ch4 MuSICA dans rmusica.
- **Stage 1 (archivage réversible)** : rmusica 321 Mo (`run_article_v323.R`, `align_ncdf.R`, `debug_musica_9Feb/`,
  `review/{_archive_stale,_bak_2026-08-18,Reu_17Jun}`) ; NC_Full 165 Mo (`new_code/`, `lidar-s2-lai-discrepancies/`,
  `paper1_sub_June26/`, 11 `.bak_PLAN_*`, 9 `CSI3_*_fable.md`). Vérification `source()` avant chaque `mv` :
  `R/h1_*`, `R/lovb_*`, `pipeline/` entier, `functions.R`, `transfer/`, `main.R`, `paper1_sub_August26/`
  sont BRANCHÉS → gardés.
- **Ch1 réparé** : carte producteur→figure à rebours des 14 figures embarquées (`PRODUCER_FIGURE_MAP_2026-09-14.md`) ;
  ledger réécrit (14 figures, 0 MISSING) ; producteurs B1/S1 perdus ré-écrits ; `run_chapter1_chs41.R`
  câblé (`--check` / `--run`, MuSICA gaté par `CH1_RUN_MUSICA=TRUE`) ; bannière SUPERSEDED sur `RUN_CHAPTER1.md`.
- **Tri c1_*** : 67 archivés en 5 lots + backups + hobo-validation (→ Ch3, tranché par Nathan) ; 18 HELD.
- **Tri rédaction** : 76 fichiers archivés (`_archive_redaction_2026-09-14/`, 51 Ch3, 20 intro/disc, 5 plans) ;
  canoniques fixés par chapitre (tableau §1).
- **NC_Full chapitré** : `scripts/ch4/` (8) → `scripts/ch3/` (32, dont 8 cachés sous `make_*`/`fig_*`) →
  `scripts/ch2/` (10 + steps 44, `phase()/step()` réécrits) ; `plans/` (14) ; `manuscripts/{ch2_rse,ch3,ch4}/`
  (30 + 8 réfs réécrites, 57 fichiers parsent) ; **bib unique** `Bib_nathan_corroyez_10Sep.bib` (5 autres
  archivées) ; `CLAUDE.md` NC_Full mis à jour (section Structure) avec accord Nathan.
- **rmusica chapitré** : `chapter1/{manuscript,producers,orchestration,ledger,_reference}` (ex-`Chapitre1/redaction`
  + 175 c1_* réunis), `Chapitre1/` supprimé ; `review/` renommé `thesis_frame/` (27 fichiers réécrits) ;
  `R/ pipeline/ in_files/` laissés à la racine (délibéré). Tests : ledger 0 MISSING, chs41 `--check` OK,
  175 producers parsent, compile bib intact.
- **Passe Fable (4 agents) + correctifs** : symlink `NC_Full/figures/article_v323` repointé ; 106 scripts
  rmusica écrivant en absolu vers `chapter3_S2_LAI|chapter4_GEDI` → `manuscripts/ch3|ch4` ; `run_chapter1.R`
  (`setwd` 2 niveaux trop bas) → `setwd(here::here())` ; guards `c1_check_*` repointés. `.gitignore` NC_Full
  étendu (thèse hors remote public).
- **Doc** : les 4 docs datés (2 par dépôt, divergés) fusionnés dans ce `HARMONISATION.md` unique + miroir ;
  README dans `NC_Full/manuscripts/ch3/`, `rmusica/chapter1/manuscript/`, `rmusica/thesis_frame/` ;
  `rmusica/CLAUDE.md` doté d'une section Structure (persona comité conservée en dessous).
- Aucun `rm`, aucun `git push`, aucune écriture dans les zones read-only.
