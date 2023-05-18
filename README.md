## About
This repo contains the data analylsis pipeline used to process the electrophysiological data found in the following publication:

**Huang JS**, Kunkhyen T, Rangel AN, Brechbill TR, Gregory JD, 
    Winson-Bushby ED, Liu B, Avon JT, Muggleton RJ,
    Cheetham CEJ. (2022). Immature olfactory sensory neurons provide behaviorally useful sensory input to the
    olfactory bulb. *Nature communications, 13, 6194 (2022).*

The project asks the question:
> **Do immature olfactory sensory neurons provide functional input to the olfactory bulb?**

Below is a brief overview of the analysis steps performed. For experimental results and source data, please refer to the paper cited above.

## Data collection
Electrophysiological recordings were obtained using the [Igor](https://www.wavemetrics.com/)-based [MIES](https://github.com/AllenInstitute/MIES) software from superficial tufted cells (TC) in the olfactory bulb during optogenetic stimulation of olfactory sensory neurons (OSN).

## Data analysis

### Strength of OSN input onto individual TCs
The amplitude of excitatory currents recorded in TCs following optogenetic activation of OSNs were quantified.

![https://github.com/janeswh/nwb_current_analysis/blob/main/figs/single_cell_1.png]

Then, the excitatory response properties were quantified.

![https://github.com/janeswh/nwb_current_analysis/blob/main/figs/single_cell_2.png]

### Average strength of OSN input in a population of TCs
Excitatory response properties were quantified and averaged across a population of TCs recorded during the same developmental timepoint

![https://github.com/janeswh/nwb_current_analysis/blob/main/figs/timepoint_summary.png]

### Comparing the strength of OSN input between mature (OMP) and immature (Gg8) OSNs
The excitatory response properties recorded in TCs were compared between two genotypes of mice, in which either mature or immature OSNs were activated.

![https://github.com/janeswh/nwb_current_analysis/blob/main/figs/genotype_comparison.png]