# Effective Fluorine Conversion Factor

The **Effective Fluorine Conversion Factor (EFCF)** provides a statistically robust relationship between organofluorine measurements and total PFAS mass concentration. Using more than **800,000 environmental samples** from Europe and North America, this work demonstrates that a single global conversion factor reconstructs total PFAS concentrations within approximately **±7%** across diverse environmental compartments, providing a practical bridge between bulk organofluorine analytical techniques (EOF, AOF, CIC) and regulatory PFAS metrics.

Code and reproducibility materials for:

> D.-A. Mendels, “Estimating Total PFAS Concentration from
> Organofluorine Measurements Using an Effective Fluorine
> Conversion Factor.”

**Status**

This repository accompanies the manuscript

> *Estimating Total PFAS Concentration from Organofluorine Measurements Using an Effective Fluorine Conversion Factor*

currently submitted for peer review.

## Overview

This repository implements the Effective Fluorine Conversion Factor
(EFCF), a statistical relationship between organofluorine concentration
and total PFAS mass concentration.

For a sample containing PFAS compounds \(i\),

`C_F = Σ C_i f_i`

`EFCF = C_PFAS / C_F`

The analysis combines molecular fluorine fractions with harmonized
environmental monitoring datasets.

## Main result

Across more than 800,000 environmental samples, the global median EFCF is:

\[
\mathrm{EFCF} = 1.548
\]

Using this value, 90% of internally reconstructed total-PFAS
concentrations fall within approximately -6.7% to +6.5% of the
reported values.

## Repository contents

src/              Python package
scripts/          Reproducible analysis workflow
data/             Compound reference library
results/          Numerical outputs
figures/          Publication figures
supplementary/    Supporting material

## Installation

```bash
git clone <repository-url>
cd PFB_Simulator

python3 -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Tested with Python 3.12.

Obtaining the environmental datasets

The environmental datasets are obtained from the PFAS Data Hub and
additional public US sources. Raw third-party datasets are not
necessarily redistributed in this repository.

Follow the instructions in:

```bash
datasets/README.md
```

Place the downloaded files in:
```bash
datasets/raw/
```

using the filenames listed in the dataset inventory.

## Reproducing the analysis

1. Run all dataset-level analyses

```bash
PYTHONPATH=src python scripts/run_all_individual_experiments.py
```

2. Run the global analysis

```bash
PYTHONPATH=src python scripts/global_analysis.py
```

3. Generate Figure 1

```bash
PYTHONPATH=src python scripts/plot_efcf_boxplots.py
```

4. Generate Figure 3

```bash
PYTHONPATH=src python scripts/plot_efcf_prediction.py
```

The final outputs are written to results/ and figures/.

### Expected headline results

* Global median EFCF: 1.547922
* Global mean EFCF: 1.548811
* Global 5th percentile: 1.452999
* Global 95th percentile: 1.658404
* Prediction-ratio interval: 0.933–1.065

Small differences may occur across software versions because of
floating-point calculations and percentile implementations.

### Data provenance

A complete dataset inventory, original sources, download locations,
sample counts, and preprocessing notes are provided in:

```bash
supplementary/dataset_inventory.csv
datasets/README.md
```

## Reproducibility note

Figure 3 is an internal reconstruction analysis. Reported total PFAS
and calculated fluorine burdens are derived from the same targeted
analyte concentrations. It evaluates the uncertainty introduced by
replacing sample-specific PFAS composition with a single global EFCF;
it is not an independent experimental validation against measured EOF
or AOF.

## Citation

See CITATION.cff.

## Licence

Source code is distributed under the BSD 3-Clause License.

Original documentation and derived research outputs are available under
CC BY 4.0 unless otherwise stated. External datasets remain subject to
their original terms and licences.

## Contributing

Bug reports, reproducibility improvements and independent validation
using additional PFAS datasets are welcome.

Scientific discussions and collaborations are encouraged.

![Python](https://img.shields.io/badge/Python-3.12-blue.svg)
![License](https://img.shields.io/badge/License-BSD--3--Clause-green.svg)
