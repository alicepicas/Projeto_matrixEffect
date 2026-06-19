# Molecular Dynamic Simulations of Molecules at Teflon Surface

> *** Towards understanding antibody adsorption in the microcapillary film (MCF) rapid diagnostic platform***

---

## Project Overview
This project establishes a robust computational framework to investigate molecular-level interactions at fluoropolymer bio-interfaces, specifically focusing on Teflon Fluorinated Ethylene Propylene (FEP). FEP microcapillary films (MCFs) are key substrates in miniaturized point-of-care (POC) diagnostics, relying on the passive adsorption of diagnostic IgG antibodies. 

The primary objective of this work is the systematic thermodynamic validation of coarse-grained molecular dynamics (MD) models against experimental wetting behaviors to enable reliable future scale-up simulations of macromolecular protein adsorption.

---

## Objectives
* **Force Field Validation:** Validate the coarse-grained representation of the water–FEP solid-liquid interface.
* **Wetting and Thermodynamics:** Quantify the static water contact angle across a temperature range (0°C to 40°C).
* **Interfacial Characterization:** Analyze interfacial water structural properties through mass density profiles.
* **Diagnostic Foundation:** Provide the essential parameterization baseline required to simulate full-length antibody orientation dynamics on polymer surfaces.

---

## Computational Methodology
The molecular simulation workflow is implemented in GROMACS utilizing advanced coarse-grained potentials to model the multi-component bio-interface:

1. **System Composition:** Solid slab of amorphous fluorinated ethylene propylene (FEP) in contact with an explicit water/solvent reservoir.
2. **Energy Minimization and Equilibration:** Rigorous NVT and NPT ensembles to stabilize density profiles and structural backbones.
3. **Contact Angle Calculations:** Evaluation of droplet geometries at specified temperature constraints to track the characteristic line tension and wetting parameters.
4. **Analysis Metrics:**
   * **Mass Density Profiles:** Tracking water layer positioning and packing near the solid boundary.
   * **Wetting Envelope Validation:** Direct comparison against experimentally established FEP contact angles (108ºC to 120ºC).

---

## Repository Structure
Based on the current project architectural layout:

```text
├── CoarseGrained_ContactAngle/
│   └── PureWater/                                          # Raw simulation data and metrics
│       ├── Sim_0C/                                         # Contact angle tracking at 0°C
│       ├── Sim_10C/                                        # Contact angle tracking at 10°C
│       ├── Sim_20C/                                        # Contact angle tracking at 20°C
│       ├── Sim_30C/                                        # Contact angle tracking at 30°C
│       └── Sim_40C/                                        # Contact angle tracking at 40°C
├── Docs/                                                   # Project Documentation
│   ├── Apresentação_Alice_Picas.pdf                        # Presentation PowerPoint
│   ├── Projeto_AdsorptionTeflon_final_Alice_ Picas.pdf     # Final Article Report
│   └── Projeto_matrixEffect_Alice_Picas.pdf                # State-of-Art Article Report
└── README.md                                               # Project documentation and landing page
```

---

## Author
Alice Picas
University of Minho - Center of Biological Engineering
