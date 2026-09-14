# Deux tests, deux verdicts opposés

*2026-09-09. Scripts : `$RMUSICA/c3_rf_bandset_test.R`, `$RMUSICA/c3_s2_decline_bands.R`.
Tables : `manuscripts/ch3/tables/TableB_rf_bandset*.csv`,
`TableB_s2_decline_bands.csv`, `TableB_s2_decline_canopy_dependence.csv`,
`TableB_s2_decline_ndvi.csv`, `TableB_s2_decline_scenecorrected.csv`. Lecture seule sur `03_RESULTS/`.*

---

## Test A — le jeu de bandes n'explique rien. La phrase d'origine tient.

La forêt aléatoire du chapitre voit dix bandes, l'inversion trois. J'avais borné la
conclusion hier pour ce motif. **Le contrôle l'invalide** : refaite sur B03/B04/B08
seules, recette identique (LOO 53 placettes, ntree = 500, seed 42), la RF fait
*au moins aussi bien* que sur dix bandes.

| Jeu de features | mtry | R² groupé | R² dense (n=25) | R² P4 (n=20) |
|---|---:|---:|---:|---:|
| 10 bandes | 2 | 0.585 | 0.355 [0.07, 0.63] | 0.311 [0.02, 0.65] |
| **3 bandes (B03/B04/B08)** | 2 | 0.528 | **0.469** [0.24, 0.70] | **0.484** [0.11, 0.77] |
| 3 bandes | 1 | 0.549 | 0.421 [0.16, 0.66] | 0.387 [0.03, 0.73] |
| 3 bandes | 3 | 0.475 | 0.471 [0.24, 0.71] | 0.484 [0.11, 0.78] |
| 10 bandes | 3 | 0.580 | 0.396 [0.10, 0.67] | 0.360 [0.05, 0.69] |

ΔR² apparié (mêmes placettes rééchantillonnées), 10 bandes − 3 bandes, strate dense :
**−0.114 [−0.37, +0.09]**. Négatif, IC incluant zéro, stable sur tous les mtry
(−0.07 à −0.12). Les sept bandes supplémentaires n'apportent rien.

**Conséquence.** Avec exactement les trois bandes que reçoit l'inversion, un
apprenant supervisé restitue le classement en canopée dense à R² ≈ 0.47, quand le
LAI restitué par PROSAIL donne 0.01. L'information est donc **présente dans les
trois bandes** et c'est bien l'inversion qui la jette. La formulation d'origine du
chapitre était juste ; ma réserve d'hier portait sur un confondant qui n'opère pas.

**La limite qui subsiste, et qui est la vraie.** La RF est entraînée directement sur
ΔTmax observé ; l'inversion passe par un LAI physique. Le contraste RF/restitution
mélange donc toujours « supervisé sur la cible » et « inversion ». C'est cette
réserve-là qu'il faut écrire, pas celle des bandes. Le contraste 10 − 3 est le seul
qui isole proprement le jeu de bandes, et il est nul.

---

## Test B — le déclin estival n'est pas une perte de feuilles.

Réflectances des dix bandes aux 53 placettes, quatre scènes 2021 (14/06, 19/07,
26/08, 22/09), aucune donnée manquante. Variation relative médiane 14 juin → 26 août :

| | B02 | B03 | B04 | B05 | B06 | B07 | B08 | B8A | B11 | B12 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| médiane % | +70 | +0.5 | **+33** | −10 | −21 | −21 | **−22** | −18 | −18 | −20 |
| accord de signe | 100 % | 51 % | 96 % | 91 % | 100 % | 100 % | 100 % | 100 % | 94 % | 94 % |

Rouge en hausse, PIR en baisse, sur la quasi-totalité des placettes : à première vue
la signature d'une perte de surface foliaire. **Deux contrôles la démentent.**

### Le facteur de scène, mesuré sur des surfaces qui ne peuvent pas changer

Entre deux passages, tout ne change pas dans la canopée : la géométrie
d'éclairement, la géométrie de visée et la correction atmosphérique changent aussi.
Des surfaces non végétalisées aux deux dates (NDVI < 0.30 le 14/06 **et** le 26/08,
B04 > 500 pour écarter l'eau) mesurent ce facteur commun directement — toits,
routes, sol nu ; les cultures sont exclues par construction, elles sont vertes en
juin. Sur 990 pixels de ce type, la baisse est **quasi identique dans les dix
bandes : −15.8 % en médiane, étalement de 6.7 points seulement**. C'est un facteur
multiplicatif de scène, pas un signal de surface.

Une fois ce facteur divisé :

| | B02 | B03 | B04 | B05 | B06 | B07 | B08 | B8A | B11 | B12 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| forêt brut % | +70 | +0.5 | +33 | −10 | −21 | −21 | −22 | −18 | −18 | −20 |
| scène % | −14 | −16 | −16 | −16 | −17 | −16 | −16 | −15 | −10 | −11 |
| **canopée seule %** | +97 | **+20** | **+58** | +7 | −6 | −6 | **−7** | −4 | −9 | −10 |

**Le visible monte fortement, le PIR ne bouge quasiment pas.** Rouge +58 %, vert
+20 %, PIR −7 % et −4 %. C'est la signature d'une perte de pigmentation, pas d'une
perte de surface foliaire : une canopée qui perd 39 % de ses feuilles fait chuter le
PIR, qui est précisément la bande où la structure se lit. Le bleu à +97 % dépasse ce
qu'une dérive de pigments explique et garde une part de brume résiduelle
(THIN_CIRRUS passe de 6·10⁻⁴ % en juin à 1.17 % en août), mais le bleu ne participe
pas à l'inversion.

### Le changement ne dépend pas de la quantité de feuilles

Régression de la variation relative sur LAI_ALS, par bande : toutes p > 0.18,
|r| ≤ 0.18. Pour le PIR, pente = 0.01 % par unité de LAI, r = 0.00. La placette
41_39 (LAI_ALS = 0.00, Hmax 4.0 m, fCover 0.00 — une trouée réelle, pas un pixel
masqué) perd 16.3 % de PIR ; la plus dense (LAI_ALS = 6.74) en perd 15.2 %. Sur une
gamme de LAI de 0 à 6.7, une perte réelle de surface foliaire ne peut pas être
indépendante de la surface foliaire présente.

*Ce second argument est le plus faible des deux et ne se suffit pas à lui-même* : au
LAI médian de 3.80 le PIR est déjà saturé, donc une perte confinée au haut du profil
donnerait aussi une pente plate. C'est le contrôle sur cibles invariantes qui porte
la conclusion ; celui-ci la corrobore.

### Ce que l'inversion en fait

Le NDVI reste dans son régime saturé (0.931 → 0.889 en médiane, baisse sur 53/53
placettes). Un déplacement radiométrique modeste y suffit à produire un grand écart
de LAI. C'est là que le déclin de 39 % se fabrique.

**Lecture.** Un facteur de scène de −16 %, plus une hausse réelle du visible à PIR
constant (pigments), sont amplifiés par une restitution qui travaille dans son
régime saturé et mal conditionné. Cela réconcilie les deux pièces du diagnostic
d'août : la dérive du rapport ATBD/opt (1.98 → 2.53 : deux LUT amplifient
différemment le même petit décalage) et la marche DOY 175→198 (SZA quasi immobile,
LAI −34 % : une amplification, pas une réponse linéaire au SZA).

**Ce que ça change.** Le creux estival n'est ni du bruit ni de la sénescence. C'est
une propriété de la restitution, dans le régime où le Ch2 la dit déjà saturée. Le
verdict du Ch3 sur les scénarios `DYN_S2_*` survit, mais sa cause doit être
renommée : dérive saisonnière de la restitution en canopée fermée, et non perte de
feuilles vue par le capteur.

**Test définitif, non lancé** (>10 configurations PROSAIL, demande ton accord) :
rejouer l'inversion sur ces mêmes dates avec un LUT dont la gamme de chlorophylle et
de pigments bruns est élargie, et voir si le déclin s'aplatit. Le test B prédit qu'il
s'aplatit.

---

## Les deux tests ne se contredisent pas

Le test A exploite la variation spectrale **entre placettes** sur une seule scène ;
le test B la variation **entre dates**. Le premier dit que les trois bandes portent
une information de structure que l'inversion jette ; le second dit que la dérive
temporelle de ces mêmes bandes n'est pas de la structure. Les deux pointent la même
chose : l'inversion, dans son régime saturé, lit mal ces trois bandes — dans un cas
elle perd un classement qui y est, dans l'autre elle fabrique une variation qui n'y
est pas.
