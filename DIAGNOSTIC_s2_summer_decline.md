# Diagnostic — le déclin estival du LAI Sentinel-2 (Blois 2021)

*Établi le 2026-08-14. Trouvé en cherchant une figure de phénologie pour l'introduction
générale, pas dans le cadre d'une analyse demandée. Non résolu : ce document décrit ce qui
est vérifié et ce qui ne l'est pas.*

## Le point le plus urgent, avant même le creux

**Aucun script du dépôt ne produit `blois_s2_lai_ts_2021.rds` ni
`blois_s2_lai_corrected_daily_2021.rds`.** `grep -rn "saveRDS.*s2_lai_ts_2021"` sur
`scripts/` et `R/` ne renvoie rien ; les huit scripts qui les mentionnent ne font que les
lire. Ils viennent d'un producteur perdu ou d'une session interactive. Ces deux fichiers
alimentent des figures du Ch3 et des tables du Ch4 d'un manuscrit à déposer dans trois mois,
et ni leur pré-traitement ni leurs réglages d'inversion ne sont vérifiables en l'état. Un
creux se diagnostique ; une entrée irreproductible, non.

## Le fait

Dans `output/intermediate/blois_s2_lai_ts_2021.rds` (60 placettes × 19 dates), le LAI
restitué par Sentinel-2 culmine à la mi-juin puis **décline d'environ 40 % jusqu'à fin
août**, avant la sénescence réelle.

| DOY | date | LAI_S2_opt | LAI_S2_ATBD | ATBD/opt |
|----:|------|-----------:|------------:|---------:|
| 165 | 14 juin | 2.37 | 4.69 | 1.98 |
| 200 | 19 juil. | 1.77 | 3.79 | 2.14 |
| 238 | 26 août | 1.44 | 3.37 | 2.35 |
| 265 | 22 sept. | 1.40 | 3.54 | 2.53 |

## Ce qui est vérifié

**Ce n'est pas du bruit.** Entre le pic (DOY 175) et le creux (DOY 225), **59 placettes sur
60** baissent, de 0.91 en médiane. Aucun pixel n'est masqué à aucune date (60/60 valides
partout).

**Ce n'est pas un artefact de mélange de provenances.** La série combine les dates
mensuelles de `scripts/download_s2_fullyear.R` et des dates d'été supplémentaires. Le
déclin subsiste **dans la série mensuelle seule** : −39 % du 14 juin au 26 août.

**Les deux inversions divergent entre elles.** C'est le diagnostic le plus solide. Le
rapport ATBD/opt passe de 1.98 (mi-juin) à 2.53 (fin septembre). Une perte réelle de surface
foliaire est commune aux deux inversions et s'annule dans un rapport ; un rapport qui dérive
signale une sensibilité de la **restitution**, pas de la canopée.

**La première phase du déclin échappe à l'angle solaire.** Entre DOY 175 et 198, le SZA ne
bouge que de 30.1° à 32.0° tandis que le LAI chute de 34 %. La géométrie d'éclairement ne
peut donc pas être la cause principale.

*(Sur l'ensemble de la fenêtre feuillée, LAI_opt ~ SZA donne r = −0.71 et le rapport
ATBD/opt r = +0.83. Ces corrélations sont confondues : après août, SZA et sénescence réelle
augmentent ensemble. Elles ne sont pas une preuve — la marche du DOY 175→198 en est une.)*

**L'explication physique est peu plausible.** Blois est une chênaie sessile ; l'été 2021 en
France a été frais et humide, sans épisode de sécheresse susceptible de provoquer une perte
de 40 % de surface foliaire en août.

## Hypothèse principale, non testée

Dégénérescence LAI–chlorophylle dans l'inversion PROSAIL. La teneur en chlorophylle décline
et les pigments bruns augmentent au fil de l'été ; si le LUT ne couvre pas cette dérive,
l'inversion peut la compenser en abaissant le LAI. Cela produirait exactement un déclin
post-pic s'accélérant vers l'automne, et une divergence croissante entre deux LUT ayant des
a priori différents sur ces pigments.

**Pour trancher** : rejouer l'inversion sur les mêmes dates avec un LUT dont la gamme de
chlorophylle et de pigments bruns est élargie, et vérifier si le déclin s'aplatit. Sinon,
examiner les réflectances par bande : si le NIR reste stable pendant que le rouge dérive,
c'est bien la pigmentation et non la structure.

## Ce qui en dépend, et ce qui n'en dépend pas

**En dépend — la fin de saison du Chapitre 4.** Vérifié en rejouant la fonction `trs50()` de
`scripts/ch4/c4_phenofit_metrics.R` sur la série observée puis sur trois contrefactuels à
plateau plat (DOY 165–265) :

| série ATBD | SOS | EOS | écart S2–GEDI |
|---|---:|---:|---:|
| observée (avec le creux) | 125 | **279** | 25 j |
| plateau tenu au pic (max inchangé, 4.69) | 132 | 279 | 25 j |
| plateau tenu à la valeur du DOY 200 (3.79) | 125 | **286** | 18 j |
| plateau tenu à la moyenne estivale | 125 | **293** | 11 j |

Le premier contrefactuel ne déplace rien, mais c'est un test faible : TRS50 fixe son seuil à
`min + (max − min)/2`, or remplir le plateau avec le pic laisse `max` inchangé. Dès que le
contrefactuel abaisse aussi le maximum — ce que ferait un plateau réaliste — **l'EOS recule
de 7 à 14 jours**. Autrement dit, entre 30 % et 55 % de l'écart S2–GEDI de 25 jours
(`Table19_c4_phenofit_trs50.csv`) pourrait provenir du creux plutôt que d'une différence
réelle de sénescence perçue. La conclusion qualitative du Ch4 (S2 décroche avant GEDI)
survit ; son ampleur chiffrée, non.

**En dépend — le forçage dynamique S2 de MuSICA.** `blois_s2_lai_corrected_daily_2021.rds`
transmet le déclin au pas journalier (2.07 au DOY 181 → 1.37 au DOY 261). Les scénarios
`DYN_S2_*` font donc voir à MuSICA une canopée qui perd un tiers de sa surface foliaire
pendant le plateau. Leur moindre performance (`Table18_c4_plateau.csv`, automne :
dT_R² = 0.10 contre 0.30 pour `DYN_ALS_*`) est cohérente avec ce déclin.

Nuance importante : si le produit S2 opérationnel se comporte ainsi, ce n'est pas un bug du
pipeline, c'est une propriété du produit, et le verdict reste valide. Mais la cause invoquée
change — ce n'est alors pas la saturation optique, c'est une dérive saisonnière de la
restitution. Cela mérite d'être dit explicitement plutôt que laissé implicite.

**Est déjà exposé.** `scripts/ch3/c3_figs_narrative.R` trace cette courbe en
`FigS1_annual_series.png` (Ch3, §3.5). Le creux est donc déjà visible dans le manuscrit,
sans explication. Un jury peut le voir.
