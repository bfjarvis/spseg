# T001 — Method and provenance inventory

Status: in review
Owner: science officer. Reviewers: validation/testing and management.
Starting revision / worktree: `de48d12a3e9d29ca09eb289873f1877ba5582117`,
canonical `spseg` checkout, inspected 2026-10-26. The framework and user's
unfinished getting-started vignette were committed during assessment as
`e2b8acac868cdca8d0119f79c5439557e6a41dff`. That commit did not change
the package source, methods vignette, network vignette, or tests inspected here.
Management's current edits to central records were left untouched.

## Brief

- Read `AGENTS.md`, `project/GOALS.md`, and `project/TEAM.md` in the current
  `spseg` checkout; inherited discussion is background, not an accepted specification.
- Map the implemented measure families to the four original papers in `project/references/`
  and definitions in the multiple-variable/network vignettes. Reuse prior research.
- Separate established definitions, implementation conventions, package supplements,
  and proposed rank-pair/network extensions. Flag unsupported claims and questions
  about asymmetric distances, disconnected populations, and normalization.
- Deliver a concise source/behavior/evidence-gap table and prioritized follow-up
  questions. Identify what needs a deeper audit rather than attempting every proof.
- Edit only this report. Do not change package code, vignettes, or manuscripts.

## Results

**Outcome.** The existing methods vignette mostly labels provenance correctly.
The essential distinction for introductory examples is between Reardon's
nonoverlapping variation ratios, their spatial local-environment analogues,
the package's finite-data conventions, and the package's additional measures.
The package computes spatial summaries from focal local compositions; it does
not currently compute a rank-pair surface, rank-distance correlation, or dyadic
KNN statistic. This is a source-and-code inventory, not a full correctness
certification.

Sources below are the four original PDFs in [`project/references`](../references/),
read with page boundaries. Page numbers refer to PDF pages, including covers:
[Reardon and Firebaugh 2002](../references/Reardon_et_al_2002_Measures_of_Multigroup_Segregation.pdf),
[Reardon and O'Sullivan 2004](<../references/Reardon and O’Sullivan - 2004 - Measures of Spatial Segregation.pdf>),
[Reardon 2009](../references/Reardon_2009_Measures_of_ordinal_segregation.pdf),
and [Reardon 2011](../references/Reardon_2011_Measures_of_Income_Segregation.pdf).
Also inspected [`spseg-multiple-variables.Rmd`](../../vignettes/spseg-multiple-variables.Rmd),
[`spseg-networks.Rmd`](../../vignettes/spseg-networks.Rmd),
`src/spseg_engine.cpp`, `src/spseg_neighbors.h`, `R/spseg_variables.R`,
`R/spseg_network.R`, and selected tests. The current `spseg` HEAD is
`e2b8acac868cdca8d0119f79c5439557e6a41dff`; all inspected package
source and principal methods/network vignettes match the starting revision.

| Family or behavior | Original source and established definition | Current behavior, classification, and evidence gap |
| --- | --- | --- |
| Nominal multigroup information, relative diversity, and dissimilarity | Reardon and Firebaugh 2002, PDF pp. 9–14, derive these as three distinct multigroup formulations. Relative diversity uses Simpson interaction; it is **one** index that reduces to the binary variance ratio. | `Model::update` and `finish` calculate entropy and Simpson ratios, plus absolute deviations. The methods vignette labels the diversity alias and binary reductions. The existing two- and three-group tests are useful reference checks, but do not certify all multigroup axioms or spatial cases. |
| Nominal square root | Hutchens' binary index is discussed in Reardon 2009, PDF p. 18, and Reardon 2011, PDF p. 22; it is **not** one of the nominal multigroup indices derived in Reardon and Firebaugh 2002. | The code uses `1 − Σᵢ αᵢ GM(pᵢ)/GM(P)` over nonempty regional groups. The vignette attributes this multigroup Atkinson form to Frankel and Volij 2011 and states its binary Hutchens reduction. The binary and three-group examples in `test-multiple-variables.R` support implementation behavior; primary Frankel/Hutchens papers were not among the four supplied PDFs and need direct source verification before a provenance claim is finalized. |
| Ordinal information, variance ratio, square root | Reardon 2009, PDF pp. 13–19, integrates binary variation over ordered cumulative category cuts and explains the three standard binary reductions. The rejected absolute-difference variant is distinct from ordinary dissimilarity (PDF pp. 18, 24–26). | C++ sums variation at cumulative cuts. Documentation correctly distinguishes these three established indices. Ordinal dissimilarity is the package's separately defined supplement: average absolute local/global CDF deviation divided by twice the sum of regional binary variances. It should never be called Reardon's ordinal absolute-difference index or assigned its axioms. |
| Rank-order information, variance ratio, square root | Reardon 2011, PDF pp. 24–30, gives threshold averages, variation-ratio representations (eqs. 15–16), the squared CDF deviation for `R^R` (eq. 18), and spatial threshold analogues (eqs. 19–21). | Grouped continuous input joins the observed local CDF knot values linearly in **global percentile space**; C++ integrates each segment rather than fitting Reardon's polynomial segregation curve. This is an explicit estimator convention. The denominators `1/2`, `1/6`, `π/4` correspond to code's entropy, `p(1−p)`, and `2√[p(1−p)]` scalings. A single populated interval still has latent rank variation under this model. |
| Exact observations and income ties | Reardon 2011 specifies a continuous population rank formulation; it does not uniquely specify finite-sample rank treatment. | Exact values get population-weighted midranks, including ties. Empirical regional CDF variation supplies the denominator; no random tie breaking occurs. Same-value count and equivalent weighted-record cases are tested. This is the package's finite-data extension, and a small data example should distinguish genuine ties from censored income intervals. Rank-order dissimilarity is likewise a package supplement, defined by integrated absolute local/global empirical CDF deviation. |
| Spatial environments, exposure and isolation | Reardon and O'Sullivan 2004, PDF pp. 9–10 and 17–20, define local composition through a nonnegative **symmetric** proximity function, then average local exposure across source-group residents and local variation across residents. The paper explicitly allows negative spatial evenness values for hyperintegration (PDF pp. 12, 28). | Current `spseg()` uses Euclidean radius or population KNN, plus local kernels; focal mass enters, KNN ties are fractionally selected, and eligible population may differ by variable. Exposure is the raw source-to-target local share matrix, including diagonal isolation; unlike Reardon and Firebaugh 2002's *normalized* multigroup exposure index, it is not a normalized evenness scalar. The package correctly avoids clamping spatial indices to `[0,1]`. Literature-to-code agreement for every kernel and denominator needs targeted reference cases. |
| Network neighborhoods | Reardon and O'Sullivan 2004 permits theoretically chosen social proximity, discusses barriers and transportation, and specifies a symmetric proximity function for the analyzed indices; it does not establish results for arbitrary directed shortest paths. | `spseg_network()` accepts prepared edges and nonnegative additive costs; traversal is outward on directed graphs; disconnected nodes are unreachable. Radius searches leave other components out of local environments but retain them in the regional reference. Population KNN errors when a target cannot be reached. These are package conventions and a proposed extension of the spatial framework. Directed paths make the proximity asymmetric, so published symmetry-dependent exchange results cannot be transferred without proof. Tests exercise direction, zero-cost ties, disconnected KNN errors, transit nodes, and a simple independent reference; a dedicated methods audit remains needed. |
| Pairwise rank-slice surface and geographic-distance dyads | None of the four papers establishes the new rank-pair slice surface, beta weighting, dyadic rank-gap/geographic-distance correlation, or its KNN profile as a package index. Reardon 2011's cumulative *one-threshold* construction and 2004's local-environment construction are the nearest foundations. | These remain research proposals in `project/papers`, outside `spseg()` and its documented outputs. An exact additive rank-cell decomposition of a proposed slice average does not by itself establish spatial normalization, bounds, interpretation, or portability to directed graphs. Treat manuscript identities and toy cases as hypotheses/evidence specific to their assumptions, not as validated package functionality. |

### Prioritized follow-up questions

1. **Tiny income cases for science and validation:** What is the intended
   population estimand for a single reported interval, a tied exact score,
   and a top-coded interval? Demonstrate each with two or three supports,
   exact expected local CDFs, scalars, and exposure output. A single interval
   should yield zero rank-order segregation under linear interpolation;
   an observed constant should yield an undefined scalar with `constant`
   status. Verify package output in an isolated build rather than inferring
   correctness from the vignette.
2. **Tiny network cases:** Use a three-node chain with population at the ends
   and a transit-only middle node to show hop radius one versus two,
   direction reversal, reachability, and global reference population. State
   the focal-inclusion and kernel-boundary rule. With the default biweight
   kernel, a one-hop boundary gets zero kernel weight; an introductory
   adjacency example should request `weighting="unweighted"`.
3. **Spatial interpretation and normalization:** Specify which 2004 spatial
   properties are claimed for radius and KNN kernels, and which require
   symmetric proximity or geographic symmetry. Keep a negative
   hyperintegration example available. For directed networks, establish a
   separate interpretation of outward exposure and decide whether an
   evenness scalar is meaningful under asymmetric access; do not silently
   symmetrize the paths. For disconnected graphs, decide whether the
   reference is the whole supplied population or a connected component;
   the current code uses the whole population while local environments may
   omit unreachable observations.
4. **Primary supplement sources:** Verify the multigroup root attribution
   against Frankel and Volij, the binary reduction against Hutchens, and
   the ordinal/continuous dissimilarity formulas against independent cases.
   Check any claimed bounds and exchange properties under *both*
   nonoverlapping and overlapping environments. The current tests are not
   a substitute for those scientific proofs.
5. **Later rank-pair methods:** Before adding a fitted rank-slice surface or
   dyadic-distance family, define the target population, pair weights,
   behavior under ties and interpolation, reference normalization, and
   interpretation at each bandwidth or KNN threshold. In overlapping or
   directed environments, the mixture of local CDFs need not equal the
   regional CDF; the paper's partition-based decomposition and `[0,1]`
   claims need a new spatial argument. This is later research, not a
   prerequisite for the bounded getting-started examples above.

**Checks and limits.** Read the four original PDFs, methods and network
vignettes, implementation branches, and named tests. No package execution,
clean-session vignette render, exhaustive property proof, or independent
network-normalization experiment was run for this report. The local
`project/references/` sources are analysis inputs, not evidence that
all attribution and code paths have been independently validated. Next:
validation should execute the tiny reference cases and developer/science
should resolve any scientific interpretation that disagrees with them;
management can then assign focused source, documentation, or code fixes.
