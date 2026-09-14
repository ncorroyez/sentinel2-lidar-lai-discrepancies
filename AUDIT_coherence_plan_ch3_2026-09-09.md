# Cohérence PLAN_ch3_FR.md ↔ Chapter3_article_standalone_EN.md

*2026-09-09. Aucune modification appliquée : la divergence principale demande un arbitrage.
Source des chiffres : `$RMUSICA/out_files/Chapter3/tables/junsep_full_by_archetype.csv`,
qui contient les deux fenêtres côte à côte.*

## 1. Le plan et le chapitre ne sont pas sur la même fenêtre d'analyse

Le plan fixe **juin-septembre** et dit pourquoi : c'est la seule fenêtre couverte à la fois
par les loggers, le vol LiDAR et les scènes S2. Le chapitre §2.1 énonce la même couverture
logger (1er juin - 30 septembre) mais **tous ses r de validation sont ceux d'avril-octobre**.
Vérifié ligne à ligne : P4 Combinaison 0.49 / LiDAR fixe 0.48 et P2 S2 0.51 sont les valeurs
`apriloct` de la table ; les valeurs `junsep` sont autres.

## 2. Le changement de fenêtre inverse deux verdicts sur quatre

Le plan affirme « le changement de fenêtre n'inverse aucun verdict ; il les clarifie ».
**C'est faux, et sa propre table le montre.** Pente de couplage, r inter-placettes :

| | scénario | avril-oct (chapitre) | juin-sept (plan) |
|---|---|---:|---:|
| **P4 dense** | LiDAR fixe | 0.48 | **0.71** |
| | Combinaison | **0.49** | 0.70 |
| | S2 seul | 0.32 | 0.34 |
| | S2 opt | −0.15 | 0.12 |
| **P2 mi-densité** | S2 seul | 0.51 | **0.71** |
| | LiDAR fixe | **0.52** | 0.63 |
| **P3** | meilleur | Combi 0.33 | Combi 0.29 |
| **P1** | S2 seul | 0.73 | 0.66 |

- **P4** : en avril-octobre la combinaison passe devant le LiDAR d'un cheveu, ce qui autorise
  la phrase du chapitre « *Only the combination matched LiDAR* ». En juin-septembre le LiDAR
  fixe reprend la tête (0.71 contre 0.70) et l'apport de S2 en dense disparaît. Le plan
  l'écrit d'ailleurs noir sur blanc : « Ne pas écrire "seule la combinaison tient" ». Le
  chapitre écrit exactement cela. **C'est le cœur de H2.**
- **P2** : en avril-octobre S2 est très légèrement *derrière* le LiDAR (0.51 contre 0.52), et
  le chapitre écrit « as well as ». En juin-septembre S2 passe nettement *devant* (0.71 contre
  0.63). Le sens de l'écart change.
- **P3** : « aucun produit n'atteint 0.35 » tient dans les deux fenêtres.
- Sur ΔTmax en juin-septembre le classement dense est le même (LiDAR 0.71, Combinaison 0.70,
  S2 seul 0.29), donc le basculement de P4 n'est pas un artefact d'une seule métrique.

**Décision à prendre.** L'argument du plan pour juin-septembre est solide et le chapitre
énonce lui-même cette couverture logger : rester en avril-octobre revient à ajuster la pente
sur six mois dont deux sans observation. Mais le passage à juin-septembre change la
conclusion dense de « la combinaison est nécessaire » à « le LiDAR fixe suffit », ce qui
touche H2, l'abstract, §3.3, §4.1 et la conclusion. C'est ton arbitrage, pas le mien.

## 3. Deux chiffres du chapitre ne correspondent à aucune des deux fenêtres

- P2 LiDAR fixe : le chapitre écrit **0.53**, la table donne 0.52 (avril-oct) et 0.63 (juin-sep).
- P4 S2 opt : le chapitre écrit **−0.21**, la table donne −0.15 (avril-oct) et +0.12 (juin-sep).

Ils viennent probablement d'un run antérieur. À resourcer quelle que soit la fenêtre retenue.

## 4. Le mécanisme d_opt du plan n'est pas dans le chapitre

Le plan pose comme mécanisme central « la différence entre le LAI LiDAR total et le LAI
ramené à la couche haute (au-dessus de d_opt = 7 m, prorata sous 10 m) ». Le chapitre n'a
aucun scénario de ce type : ses quatre scénarios varient la *magnitude*, jamais la
*profondeur*, et « Combinaison » désigne la saisonnalité S2 sur magnitude LiDAR. d_opt
n'apparaît qu'en §2.3 (une phrase qui fixe 7 m) et en Fig. 4 comme covariable descriptive.
La décomposition existe pourtant dans le code (`scripts/ch3/c3_hybrid_lai.R` : S2 haut + LiDAR
bas) et le produit tronqué `STATIC_S2_DOPT` existe aussi. Soit le chapitre l'intègre, soit
§2.3 doit dire explicitement que la question haut/bas est traitée par la Fig. 4.

Deux précisions que le plan demande de porter au texte et qui n'y sont pas :
- d_opt = 7 m est la valeur inter-sites ; **Blois seul donne 6 m** (le chapitre dit seulement
  « within the 6-10 m companion range »).
- **Le prorata sous 10 m n'est pas implémenté** (`build_real53_df.R:26`) : pour Hmax ≤ 7 m la
  fraction vaut 1. Rien dans le chapitre ne le signale.

## 5. Points mineurs

- **Strate de la RF** : le plan précise que le n = 25 est l'ancienne strate LAI ≥ 3.86, **pas
  P4**. Le chapitre donne n = 25 sans le dire, dans un texte par ailleurs entièrement structuré
  par archétype. Une clause suffit. (Le contrôle du jeu de bandes du 09/09 hérite de la même
  strate.)
- **Groupé 0.81 / 0.70** : le plan les marque « tables avril-octobre, à recalculer sur
  juin-septembre » ; le chapitre les présente comme un résultat établi (§3.3, §4.1).
- **Journal du plan** : la section « Corrections déjà appliquées au manuscrit » ne liste que
  les deux phrases du matin ; les modifications de l'après-midi (contrôle jeu de bandes,
  paragraphe saisonnalité §4.3) sont décrites ailleurs dans le plan mais pas là.

## Ce qui est cohérent

Objectif et question · H1/H2 · les quatre scénarios et leur définition (structure LiDAR dans
tous les scénarios S2) · P1 = trouées, LAI S2 faux de ~4 m²/m², deux placettes sans série S2
sur huit · P3 stratum faible · métrique primaire = pente, ΔTmax secondaire · robustesse
densité de points (ρ ≥ 0.99) · FORMS-H et GEDI en discussion · modèles 3D/bords en
perspective · jeu de bandes B03/B04/B08 (corrigé des deux côtés le 09/09).
