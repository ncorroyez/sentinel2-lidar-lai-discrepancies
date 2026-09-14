# Synthèse — Méthodes, résultats, et rôle de d_opt (Ch2 ↔ Ch3)

*Consolidation de la session de juin 2026. Sites : Blois (chênaie), Aigoual
(hêtraie montagnarde), Mormal (chênaie-hêtraie). Validation micro-climat : HOBO
sub-canopée vs MuSICA. LAI : ALS (LiDAR, Bouvier 2015) vs S2 (PROSAIL hybride SVR).*

---

## 0. La question

`d_opt` = profondeur de canopée effectivement « vue » par S2 (Ch2, métrique de
**consistance inter-capteur**, ~6/7/8/10 m selon site). **Améliore-t-il
l'estimation du LAI et/ou le forçage micro-climatique de MuSICA (Ch3) ?**

---

## 1. Cause racine établie : **saturation S2**

Le LAI S2 (toute paramétrisation : ATBD, opt) **sature** en forêt dense :

| | sd | cor(LAI_plein) | cor(Hmax) | plateau > LAI 4 |
|---|---|---|---|---|
| S2 (ATBD/opt) | 0.50–0.69 | 0.19–0.29 | −0.33 à −0.45 | oui (moy ≤4 ≈ >4) |
| LiDAR (réf) | **1.57** | — | — | non |

→ S2 seul prédit le LAI plein à R²≈0.02 ; anti-corrélé au LiDAR aux HOBOs (−0.27).
**Conséquence structurante : la magnitude DOIT venir du LiDAR ; S2 n'apporte que
le _timing_ (phénologie).** Aucun LUT ne corrige (limite optique).

---

## 2. Méthodes construites — Pipeline Part 1 (correction LAI temporelle)

Série S2 LAI annuelle Blois 2021 (le facteur limitant — 10 dates été — a été levé) :

```
Réflectance S2 annuelle (19 dates ; fév/avr/nov/déc téléchargées + processées)
  ↓ inversion SVR existants (ATBD + opt)            — pas de ré-entraînement
Série LAI_S2 par date × plot
  ↓ Whittaker cyclique pondéré (λ≈800)              — lisseur satellite retenu
Forme journalière lisse + gap-fill                   (> GAM, SavGol, double-log)
  ↓ fraction de verdissement                         — fix plancher hiver (→0)
  ↓ × magnitude { LAI_ALS · LAI_ALS_dopt · RF-ML }   — design matrix
LAI journalier corrigé par-plot, vraie phéno annuelle
```

- **Whittaker** = LE lisseur (tunable, gap-fill cyclique). SavGol colle au bruit,
  double-logistique trop rigide (mais donne les dates : leaf-out DOY 120,
  sénescence DOY 290).
- **Fix plancher hivernal** : recalage multiplicatif gonflait l'hiver (1.0–1.4) →
  résolu par _fraction de verdissement_ (hiver≈0.2, été≈LAI_ALS).
- **ML** reformulé : RF pour la **magnitude** (structure+S2, OOB R²≈0.69) ×
  **forme** Whittaker. Aux plots, `ML ≈ rescale_ALS` (on a le LiDAR) → **valeur
  du ML = hors-LiDAR uniquement**.
- Artefacts : `output/intermediate/blois_perplot_daily_v2_2021.rds`,
  figures `output/figures/blois_*_2021.png`, scripts download/inversion.

---

## 3. Résultats micro-climat clés (indépendants de d_opt)

- **Le buffering est piloté par le LAI TOTAL.** Slope horaire observée
  `Tmicro ~ Tmacro` : **cor(slope, LAI_plein) = −0.92** (vs LAI_dopt +0.24).
  MuSICA + LiDAR plein la reproduit à **r = 0.92**.
- **La fenêtre de validation change le verdict.** Été seul (canopée plate) →
  STATIC gagne (le dynamique = bruit). **Année complète → `DYN_RF` gagne**
  (RMSE 2.19 < tous statiques) : S2 capte le timing réel leaf-on/off.
- **Meilleur forçage dynamique = `DYN_RF`** (RF cible=LAI_ALS plein, garde l'obs
  S2 journalière). = **magnitude LiDAR + timing S2**.
- **Biais MuSICA structure-dépendant** (cor(biais, LAI)=+0.82) = levier
  non-d_opt (post-correction structure).

---

## 4. d_opt — ce qu'il **NE FAIT PAS** (forçage micro-climat)

Testé sous **13 angles** → aucun apport indépendant au forçage MuSICA :

| # | angle | résultat |
|---|---|---|
| 1 | cible de forçage (STATIC_ALS_DOPT) | ❌ magnitude 1.6 trop basse, biais +2-3 °C |
| 2 | ancre de rescale (d7 / d20) | ❌ monotone, optimum = pleine magnitude |
| 3 | RF cible = LAI_ALS_dopt | ❌ hérite magnitude basse |
| 4 | feature in-sample (LOO HOBO) | ❌ redondant (contrôle bruit = idem) |
| 5 | feature out-of-sample (cLHS→HOBO) | ⚠️ −8 % récup LAI MAIS **nul en MuSICA** |
| 6 | paire cohérente (S2_DOPT × LAI_dopt) | ❌ pire tiers, toutes fenêtres |
| 7 | physique S2 × Hmax/d_opt | ❌ d_opt constant **absorbé** par la calib |
| 8 | LAD tronqué à d_opt (LADOPT) | ❌ mid-pack |
| 9 | résidu Tmax ~ d_opt | ❌ secondaire, dominé par structure brute |
| 10 | résidu slope ~ d_opt | ❌ redondant (compacité = 1/LAI_plein) |
| 11 | split 2-couches (slope) | ❌ haut/bas buffent pareil (F p=0.58) |
| 12 | jours les + chauds | ❌ s'enfonce (écart se creuse) |
| 13 | cross-site **harmonisé** (d_opt constant) | ❌ absorbé |

**Raison de fond :** d_opt optimise la **consistance inter-capteur** (objectif
Ch2), **orthogonale** à l'**exactitude de forçage micro-climat** (besoin Ch3).
Le micro-climat (offset ET slope) dépend du **LAI total** ; la faible magnitude
d_opt (canopée sous-d_opt jetée) coûte. Et avec **LiDAR partout**, on lit la
magnitude LiDAR direct → d_opt redondant pour le forçage.

---

## 5. d_opt — ce qu'il **FAIT**

### (a) Diagnostic Ch2 — il quantifie la limite d'observabilité de S2
d_opt mesure **jusqu'où S2 voit** (~7 m du haut) et où il sature. Découpe la
canopée en **visible-S2** (top d_opt, LAI≈1.6) vs **invisible-S2** (sous d_opt,
LAI≈2.5, > la moitié).

### (b) Il PROUVE que la partie S2-invisible compte pour le micro-climat
Régression `slope ~ LAI_top + LAI_below` : **β_below = −0.066 (p<1e-4) ≈
β_top = −0.057** → la couche **sous-d_opt buffe autant** que le haut et en
contient plus. Donc S2 omet structurellement une couche buffer-active →
**justifie d'utiliser le LiDAR pour la magnitude.** C'est le **pont Ch2→Ch3
conceptuel**.

### (c) Descripteur structurel **transférable cross-site** (résultat positif)
Leave-one-site-out, prédire le LAI plein du site tenu à l'écart :

| site test | structure seule | **+d_opt** |
|---|---|---|
| **Aigoual** (hêtraie) | RMSE 2.23 / R² 0.36 | **1.36 / 0.65** |
| Blois | 1.50 / 0.08 | **1.15 / 0.34** |
| Mormal | 1.31 / 0.30 | **1.17 / 0.49** |

d_opt améliore le transfert **dans les 3 plis**, fort sur la structure la plus
différente (Aigoual). Mécanisme : il porte le **niveau absolu LiDAR**, et le
**ratio plein/d_opt est cross-site-stable** (vs relations S2→LAI site-spécifiques).
⚠️ **Caveat** : d_opt reste LiDAR-dérivé (circulaire pour un vrai « sans LiDAR »).
**Piste hors-circularité = GEDI** (top-canopée ≈ d_opt sans ALS dense) — perspective,
non développée.

---

## 6. Verdict — l'histoire d_opt, complète et honnête

| niveau | rôle de d_opt |
|---|---|
| **Forçage micro-climat (Ch3)** | **aucun** rôle indépendant (13 angles, null) |
| **Récup LAI in-situ (avec LiDAR partout)** | redondant (on a le LAI plein) |
| **Carte spatiale / transfert cross-site** | **rôle quantitatif démontré** (R² ↑, surtout cross-structure) |
| **Cadre conceptuel** | **central** : explique la saturation S2 → justifie le LiDAR-magnitude (pont Ch2→Ch3) |

**En une phrase :** d_opt est un **résultat de consistance inter-capteur (Ch2)**
qui **diagnostique pourquoi S2 ne suffit pas** (sature, ne voit que le haut) et
**justifie le pipeline LiDAR-magnitude + S2-timing** ; il **n'améliore pas le
forçage** (le LAI total, lu sur LiDAR, domine), mais c'est un **descripteur
structurel transférable** utile pour la cartographie cross-site.

---

## 7. Suite

- **Part 2 = validation MuSICA** de la série annuelle (vraie phéno observée) →
  le dynamique annuel bat-il enfin le statique aux HOBOs ?
- Mineur : 3 dates S2 manquantes (Blois mars), passage raster-domaine.
- Perspective : lien d_opt ↔ GEDI (non prioritaire).
