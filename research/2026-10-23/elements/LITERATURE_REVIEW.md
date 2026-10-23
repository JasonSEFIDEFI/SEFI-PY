# SEFI:ELEMENTS — focused literature review

Prepared October 23, 2026. This review identifies relevant prior work and the standards a geometric nuclear prediction must meet. It is a focused review, not a systematic review with exhaustive database coverage.

## Evaluated data and what a table entry means

Huang et al., *The AME2020 atomic mass evaluation (I)*, describe the evaluation and adjustment process; Wang et al., Part II, provide recommended masses and derived quantities. These are evaluated estimates built from experimental information and constraints, rather than a set of statistically independent direct measurements. The electronic table distinguishes extrapolations with `#`. Our analysis excludes those extrapolated mass entries from the predictive scoring set. Evaluation uncertainties and their correlations must not be confused with model-prediction errors. [Part I](https://doi.org/10.1088/1674-1137/abddb0), [Part II](https://doi.org/10.1088/1674-1137/abddaf)

Kondev et al., *The NUBASE2020 evaluation of nuclear physics properties*, supply decay and state information, including ground states, isomers, and estimates for some unobserved nuclides. A database row is therefore not automatically evidence of observation. We retain state and estimation flags and use the ground-state identifier in the join. A `stbl` label is handled as no observed decay in this evaluation, not an infinite numerical training target. Bounded, estimated, and finite lifetimes need distinct handling. [Evaluation](https://doi.org/10.1088/1674-1137/abddae)

The data sources explicitly request citation of their original evaluation papers. The present package records the downloaded versions and file hashes. Using AME2020 does not imply a claim that its contents exhaust knowledge available in 2026. [AMDC data guide](https://amdc.impcas.ac.cn/web/masseval.html)

## Conventional structural models already contain geometry and shells

Duflo and Zuker constructed microscopic mass formulae using nuclear scaling and shell-model structure. Their work establishes that a simple liquid-drop formula is an introductory baseline, not the strongest conventional competitor. We fit a transparent liquid-drop model here to establish a working protocol; an improvement against it does not demonstrate competitive nuclear-mass prediction or a novel physical interaction. Supplied shell closures in a feature vector cannot be presented as rediscovery of those same closures. [Original paper](https://arxiv.org/abs/nucl-th/9505011)

Nuclear binding and chemical periodicity are different targets. A nucleus is identified by proton number, neutron number, and state, while a periodic-table element fixes proton number and admits many isotopes. Predicting electron configurations, chemical groups, or periodic lengths requires an electronic structure argument. A successful nuclear mass fit alone would not derive those chemical properties. We therefore frame the first computational work as an isotope-chart study, rather than a derivation of the entire periodic table.

## Statistical residual models and honest extrapolation

Neufcourt et al. used Bayesian residual models for nuclear observables and evaluated extrapolation to later-measured nuclei, including uncertainty coverage. This is directly relevant to the proposed geometry: a residual surface is already a conventional statistical construction. The required comparison is against a model with the same information and tuning resources, not only against an uncorrected formula. Calibration coverage on unfamiliar regions matters alongside average error. [Paper](https://arxiv.org/abs/1806.00552)

Yüksel, Soydaner, and Bahtiyar investigated support-vector and Gaussian-process mass prediction with physics-based inputs and extrapolation tests. Their study reinforces the need to distinguish interpolation among measured neighbors from prediction beyond the measured region. Our random, proton-blocked, and neutron-frontier splits answer different questions; favorable random-split results cannot be transferred directly to superheavy nuclei. [Paper](https://arxiv.org/abs/2401.02824)

Huang et al. studied physics-informed networks and evaluated AME2016-to-AME2020 predictions. We therefore add an actual earlier-evaluation test: AME2016 fitting targets and 73 newly non-extrapolated AME2020 entries, with no nuclide overlap between the fitting pool and later test cohort. The separate discovery-year split remains a proxy, because nuclide discovery time is not mass-measurement time. Neither retrospective protocol is an independently blinded prospective discovery. [Preprint](https://arxiv.org/abs/2501.01352)

## Directly relevant recent latent-geometry work

Richardson, Trifinopoulos, and Williams, published October 21, 2026, report interpretable learned representations for nuclear masses, including structured proton and neutron embeddings. This makes latent geometric organization a relevant existing research topic. A SEFI manuscript must explain what its equations add beyond learned representation, regularization, or re-expression of conventional models. The present prototype does not reproduce or benchmark their full model and makes no performance comparison against it. [Published paper](https://www.nature.com/articles/s44488-026-00017-9)

## Superheavy realization requires formation and decay physics

Ackermann reviews experimental decay spectroscopy of heavy and superheavy nuclei. This work anchors the distinction between a theoretical chart location and an experimentally identified nuclear species. Mapping an empty coordinate does not supply a measured decay chain or identify proton number. We do not claim a new element from the present data analysis. [Review](https://arxiv.org/abs/2501.04053)

Simenel, Guo, and Sekizawa review microscopic and hybrid descriptions of superheavy formation and decay. The necessary physics includes dynamical and structural effects absent from our kernel smoother. A new-element candidate would need more than a mass estimate: relevant reaction accessibility, decay competition, deformation and barrier uncertainties also matter. The current package contains no reaction synthesis prescription or prediction of a production rate. [Review](https://arxiv.org/abs/2609.25621)

## Consequences for this paper

The defensible first contribution is a transparent reconstruction, a mathematical audit of its inference, and an openly reproducible held-out assessment. A later advantage claim must beat modern, matched conventional models under predefined extrapolation tests, with validated uncertainty. A later existence claim needs a physical bound-state or resonance argument and ultimately observation. Neither follows just from a smooth embedding or a density gap.

The original personal SEFI periodic assignment was not saved. We therefore cannot establish that the current reconstruction reproduces it, nor that it predates related published representations. No priority claim is made.
