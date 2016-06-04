# Gastric-Emptying-BE
Exploring the Effects of Gastric Emptying and MMC-dependent Transit on Bioequivalence for BCS Class I and III Drug Products

© 2016 A. Talattof

Fasted-state cyclical contractions associated with gastrointestinal motility transport content from the stomach into the small bowel in a pulsatile manner, influencing the concentration presented at the absorption site. In the case of high solubility/permeability (BCS Class I) drugs, this physiological factor presents a source of inter- and intra-subject variability affecting plasma levels, thus having bioequivalence (BE) implications.

A physiologically-based model was constructed to account for cyclical gastric emptying rates. The time-dependence was represented by a periodic, piecewise-continuous function that increased from phase I through III. Assuming fasted state and linear metabolism, simulated emptying rates corresponded to in vivo studies [Oberle et al., 1990]. The range of volumetric emptying was also evaluated against results from in vivo studies [Mudie et al., 2013]. The variations in Cmax and Tmax were calculated relative to the randomly chosen dosing times. The population reference was calculated using 10000 simulations. Samples of 6, 12, or 24 virtual subjects were randomly chosen and evaluated as pilot or BE trials.

Approximately a quarter of the subjects in the volumetric emptying study displayed non-first order emptying and 20% of the simulations showed similar kinetics. For BE studies, approximately 47% of the FDA-required mean Cmax 90% confidence intervals (CI) of the samples exceeded the reference simulation 80-125 range for a 7-min. half-life drug; 33% for a 14-min. half-life; and 7% for a 30-min. half-life. For 1-hr. and 4-hr. elimination half-lives the sample CI constituted 28% and 10% of the 80-125 population ranges, respectively.

Files:
* massBalanceBE.py
  - Class that defines the mass balance setup: serial array of 8 continuously-stirred reactor tanks. Transit depends on the MMC cycle relative to dosing time (t0).
  - Create class with the following parameters:
    - Kpel: elimination rate
    - PermRates: permeation rates for each of the compartments
    - KgeParams: parameters that define MMC-dependent transit rates
    - LagParams: parameters that define the MMC-associated lag-times
    - vol: initial dose volume
    - M0: initial value of matrix (initial state for the fluid volume, solid dosage form, disintegrated dosage, and dissolved drug)
    - title: string for naming of the output files
* BCS[#][#].py: Analysis of BCS Class I or III compound administered with either 50 or 200 mL fluid volumes.
  - This will run 5 simulations of 24 subjects, each iterating over the 120-minute dosing time range.
  - Outputs concentration time-profiles, as well as population Cmax, AUC, and Tmax values.

Note: These scripts rely on the following libraries
* numpy, scipy, simpy, time, sys, joblib, multiprocessing, csv, pickle, bz2, contextlib

Further reading:
* Gastrointestinal Motility Variation and Implications for Plasma Level Variation: Oral Drug Products.
Talattof A, Price JC, Amidon GL.
* http://www.ncbi.nlm.nih.gov/pubmed/26692042
