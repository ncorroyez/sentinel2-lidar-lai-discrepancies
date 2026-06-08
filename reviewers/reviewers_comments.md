# RSE-D-25-04417 — Reviewers' Comments

**Decision date**: 2026-04-03  
**Decision**: Major revision  
**Due date for revised manuscript**: 60 days (i.e. ~2026-06-03)  
**Editor**: Jing M. Chen (Editor-in-Chief)  

Manuscript title: Assessing and Reducing Discrepancies Between Sentinel-2 
and LiDAR-Derived LAI in Temperate Deciduous Forests

---

## Reviewer #2

Overall assessment: thorough investigation, interesting findings, but 
substantial issues to address before acceptance.

### R2.1 — Research gap and contributions
The research gap and the contributions of the study need to be 
strengthened. Authors are encouraged to articulate more explicitly which 
aspects of Sentinel-2 vs LiDAR LAI inconsistencies have already been 
addressed and which underlying mechanisms have not yet been systematically 
quantified. Contributions should be clearly stated.

### R2.2 — Rationale for temperate deciduous forests
The rationale for focusing on temperate deciduous forests requires further 
justification. Are these forests more representative, more challenging, or 
particularly suitable for the research questions?

### R2.3 — Study sites selection
The reasons for selecting the three specific study sites are not 
sufficiently explained. More detailed information should be provided in 
the study area section.

### R2.4 — CHM-based normalization rationale
DTM-based normalization is common practice and serves to remove terrain 
effects. However, the rationale for CHM-based normalization requires 
further clarification. LAD is typically derived either using DTM-based 
normalization or without any normalization; the physical basis for 
CHM-based normalization should be more clearly justified.

### R2.5 — Mormal temporal mismatch
In the Mormal study area, there is a relatively large temporal mismatch 
between the ALS and Sentinel-2 acquisitions. The potential impact of this 
mismatch on the results should be discussed.

### R2.6 — Section numbering inconsistency
Manuscript structure appears inconsistent: section 2 is Materials and 
Methods, followed directly by Results as section 4. Clarify whether 
section 3 is missing or incorrectly numbered.

### R2.7 — 270 LUT configurations
The use of 270 LUT configurations may raise concerns regarding trial-
and-error parameter optimization. Authors should explain the prior 
knowledge or theoretical basis used to define the parameter ranges.

### R2.8 — CHM std vs DSM std as heterogeneity metric
CHM standard deviation was used as the metric for horizontal heterogeneity. 
Why not DSM? In principle, the standard deviation of DSM may better capture 
structural heterogeneity. The criteria used to define low, medium, high 
heterogeneity classes should be clearly explained.

### R2.9 — Blois weaker relationships
The Blois forest consistently exhibits weaker relationships than the other 
study sites. This deserves further discussion.

### R2.10 — LMA cross-site transferability
LMA shows cross-site transferability. It would be valuable to identify 
which PROSAIL parameters require site-specific calibration and highlight 
this as a study conclusion. A dedicated discussion on the stability of 
individual PROSAIL parameters across sites is recommended.

### R2.11 — Local parameter adjustment interpretation
For parameters that require local adjustment, a brief interpretation is 
encouraged (species composition, leaf traits, stand age structure).

### R2.12 — Heterogeneity mechanisms
Current discussion interprets heterogeneity effects as violations of the 
homogeneous canopy assumption. Elaborate further on how horizontal 
heterogeneity affects reflectance through shadowing and gap fraction 
variations. Areas with strong height heterogeneity often coincide with 
complex vertical structure; these two factors may jointly amplify 
Sentinel-2 retrieval errors.

### R2.13 — Three-line table style
All tables should be reformatted using a three-line table style.

### R2.14 — Future research directions
Future research directions could be made more specific by outlining clear 
next steps for extending or validating the proposed framework.

### R2.15 — Data and code availability
Authors are encouraged to disclose data availability and consider sharing 
relevant code.

---

## Reviewer #3

Overall assessment: important and timely topic, well-designed workflow, 
clearly written. One major concern on LAI_ALS robustness, several minor 
comments.

### R3.MAJOR — LAI_ALS robustness
A key interpretation is that S2 underestimates LAI primarily because 
passive optical reflectance is mainly sensitive to an upper canopy 
thickness in dense stands. This relies on LAI_ALS acting as a robust 
reference. However, the manuscript reports LAI_ALS values reaching ~15 
m²/m² in Aigoual and Mormal, while LAI_S2_ATBD saturates at ~6–7 m²/m². 
Such high LAI_ALS values are unusual for temperate deciduous stands and 
raise the possibility that part of the apparent S2 underestimation may 
reflect methodological inflation in LAI_ALS.

Specific concerns:
- Fine-support, hit-based nature of ALS LAD/LAI retrieval: at 100 m² 
  support, finite footprint and geolocation uncertainty allow small 
  canopy elements to generate returns. With hit-based counting, sparse 
  intercepts can disproportionately reduce estimated gap fraction.
- The log transform in Eq. 1 propagates small errors in estimated 
  transmittance into systematically large LAD values.
- Wide scan angle and high point density increase detection of small 
  components.
- Extinction coefficient fixed at k = 0.5 without sensitivity analysis. 
  k is sensitive to canopy architecture and clumping.
- The manuscript should demonstrate that main conclusions are robust to 
  reasonable uncertainty in k and to support/aggregation choices in LAD 
  computation.

Without these robustness checks, it remains unclear whether the reported 
mismatch and derived d_opt primarily reflect S2 penetration/saturation 
or ALS methodological sensitivity.

### R3.minor.1 — ALS not ground truth
The manuscript states ALS "leads to accurate LAI estimates" in dense 
vegetation. Balance this by explicitly stating limitations of ALS-based 
LAI/PAI (scan angle, occlusion, clumping, return-processing choices), 
not only temporal coverage. Avoid implying ALS LAI is truth.

### R3.minor.2 — Phenological state
Given deciduous forests, phenological state matters for both S2 
reflectance and LiDAR foliage conditions (even if acquisitions are 
close) in Figure 1.

### R3.minor.3 — LUT parameterization mechanism
Clarify why LUT parameterization is considered a key source of ALS-S2 
discrepancy. Add 1–2 sentences explaining the mechanism (forests violate 
PROSAIL assumptions; LUT priors constrain inversion).

### R3.minor.4 — Thresholds sensitivity
Justify the "2 m" vegetation threshold and the ">90%" cover filter. 
Provide sensitivity analyses (e.g., 3–5 m for h_min; 80–95% for fCover).

### R3.minor.5 — Name PROSAIL explicitly
In methods, explicitly name "PROSAIL" when describing PROSPECT-D + 4SAIL.

### R3.minor.6 — LAD_DTM vs LAD_CHM explanation
Add 3–4 lines explaining "ground-up vs top-down depth coordinate," and 
why complex topography or rough canopy-top makes them diverge.

### R3.minor.7 — Report R² alongside r
Reporting R² alongside r would help many readers.

### R3.minor.8 — Sites combined d_opt = 4 m
Explain why "sites combined" yields d_opt = 4 m, smaller than individual 
sites. Explain whether the combined result is dominated by one site's 
distribution, sampling scheme, or height filtering.

### R3.minor.9 — d_opt sampling clarification
In 2.5.4 you describe sampling 5000 pixels with uniform LAI_ALS 
distribution; clarify if this is done per site or pooled, and how it 
is handled for "sites combined."

### R3.minor.10 — Citation Wan et al. 2024 RSE
Cite and discuss Wan et al. (2024, RSE) on canopy structure, biochemistry 
and soil background confounding in Sentinel-2 LAI.

### R3.minor.11 — Shadowing citations
If you argue residual discrepancies come from crown geometry and shadowing 
effects, provide appropriate references.

### R3.minor.12 — Heterogeneity metric window
CHM std may not isolate canopy-top heterogeneity from terrain/surface 
effects. Report the effective window explicitly (10 m unit) and clarify 
the CHM smoothing effect. Consider testing DSM variability or other 
surface roughness metrics to separate topography from canopy.

### R3.minor.13 — d_opt multi-criteria optimization
d_opt is chosen by maximizing r, and later you note trade-offs among 
metrics. Additionally test d_opt defined by minimizing RMSE or |bias|, 
or multi-objective selection, especially because your discussion 
acknowledges correlation-driven trade-offs.

---

## Reviewer #4

Overall assessment: interesting study, tries to explain S2 variations 
using upper-canopy structure, but weaknesses need to be addressed.

### R4.1 — No field validation
The biggest concern: the study focuses on improving agreement between two 
remote-sensing-derived products, without validating either against 
independent field measurements. The study demonstrates improved 
consistency between two proxy estimates, but we don't know how much of a 
clear improvement in LAI retrieval accuracy.

### R4.2 — Conclusions too broad
Conclusions are stated too broadly relative to the evidence. The analysis 
is based on only three temperate deciduous forests in France under leaf-on 
conditions. Not yet sufficient to support broader claims regarding forest 
LAI retrieval in general.

### R4.3 — Abstract and introduction too long
Abstract and introduction parts are currently too long. Prioritize 
research question, current gaps, and scientific significance. The 
introduction contains a large amount of general background material; the 
specific research gap could be presented more directly.

### R4.4 — Vertical canopy structure importance
The importance of vertical canopy structure for interpreting LAI_S2 should 
be explained more clearly. The introduction states this is important for 
interpreting LAI_S2, but the explanation and discussion are not yet 
sufficiently clear.

### R4.5 — Pixel/site/forest framing
Main analyses are at 10 m pixel level and then compared across three 
study sites. But parts of the manuscript are framed in terms of forest-
level interpretation. Clarify the link between pixel-level analysis, 
site-level comparisons, and forest-level conclusions. Explain what is 
specifically gained by framing the study at the forest scale.

### R4.spec.1 — L376 reference-driven optimization
The LUT design includes LAI priors derived from LAI_ALS and LAI_ALS_dopt, 
and model performance is then evaluated against LAI_ALS_dopt. This makes 
the optimization partly reference-driven; the reported improvement is 
better interpreted as improved consistency with ALS, rather than clear 
evidence of improved LAI accuracy.

### R4.spec.2 — L356 d_opt criteria
The optimal depth is selected solely from the maximum Pearson correlation. 
Could combining other criteria lead to different depths?

---

## Editorial requirements (from decision letter)

- E.1: Cover letter summarizing major changes
- E.2: Itemized reply to all comments with original comments quoted verbatim
- E.3: Track-changes version of revision
- E.4: Clean revision
- E.5: Research Highlights max 85 characters per bullet
- E.6: High-resolution figure files (eps, tiff, jpeg, high-quality pdf)
- E.7: List of figure captions after references
- E.8: Optional KML files for geographic visualization

---

## Cross-reviewer overlaps

Some reviewer comments converge on the same issues. These are priority 
items because addressing them satisfies multiple reviewers at once:

- **Multi-criteria d_opt selection**: R3.minor.13 + R4.spec.2
- **DSM std vs CHM std for heterogeneity**: R2.8 + R3.minor.12
- **Reference-driven / consistency vs accuracy framing**: R4.1 + R4.spec.1
- **Heterogeneity mechanisms (shadowing, gap fraction)**: R2.12 + R3.minor.11
- **Sensitivity thresholds (h_min, fCover)**: R3.minor.4 (unique but mentioned here for tracking)
- **k sensitivity**: R3.MAJOR (unique but the most technically demanding)