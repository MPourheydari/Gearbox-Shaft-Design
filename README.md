# Gearbox Shaft Design: Fatigue-Based Sizing in MATLAB

Design of three parallel shafts for a can-making machine gearbox, with a MATLAB solver that sizes the fatigue-critical shaft sections (modified Goodman criterion) and updates shoulder-fillet radii automatically.

> **Course project** – Machine Design I, K. N. Toosi University of Technology (KNTU), Winter 2025–2026. Individual project, July 2026.
> This is a **preliminary design** for educational purposes, not a production-ready design.

---

## Project at a glance

| Item | Value |
|---|---|
| System | 3-shaft spur-gear gearbox (S1, S2, S3) with belt-pulley input and coupling |
| Motor | 21 kW at 1882 rpm |
| Power split | 40 % through coupling C1 (steady); 30 % each to gears G5 and G6, varying periodically (period 21 s) |
| Material | S<sub>ut</sub> = 630 MPa, S<sub>y</sub> = 370 MPa |
| Design factor of safety | n = 2 |
| Input shaft S1 | Static design, von Mises criterion (see report, Section 6.1) |
| Shafts S2 and S3 | Fatigue design, modified Goodman criterion (solved in MATLAB) |
| Tools | MATLAB (base, no toolboxes) |

## Repository contents

```
gearbox-shaft-fatigue-design/
├── README.md                          This file
├── Gearbox_Shaft_Design_Report.pdf    Full report (18 pages): loads, criteria, method, discussion
├── Gearbox_Shaft_Design.m             MATLAB solver
├── Shaft_Design_Input_Data.csv        Section loads and initial guesses (solver input)
└── Shaft_Design_Results.pdf           Command Window output of the final run
```

## Results

Final run of the solver (shafts S2 and S3). `D_std` is the next standard diameter at or above the calculated value; `r_std` is the preferred fillet radius, never larger than the geometric limit |D<sub>left</sub> − D<sub>right</sub>| / 2.

| Physical shaft | Section | D_calc (mm) | D_std (mm) | r_std (mm) |
|:--:|:--:|--:|--:|--:|
| S3 | B5_G4 | 28.74 | 30 | 1.2 |
| S3 | G4_B6 | 31.15 | 32 | 0 |
| S3 | B6_G5 | 31.15 | 32 | 0 |
| S2 | C1_B3 | 29.10 | 30 | 2.5 |
| S2 | B3_G2 | 34.74 | 35 | 0.8 |
| S2 | G2_G3 | 36.50 | 40 | 1.5 |
| S2 | G3_B4 | 33.28 | 35 | 3.0 |
| S2 | B4_G6 | 27.12 | 28 | 0 |

Numerical summary: outer loop converged in 7 passes (final max. diameter change 0.073 mm); every section needed 10 bisection iterations (80 in total). The complete tables (loads, K<sub>t</sub>, K<sub>f</sub>, q, k<sub>b</sub>, S<sub>e</sub>, iteration history) are in `Shaft_Design_Results.pdf`.

> The CSV labels `Shaft1` and `Shaft2` are computational group names. `Shaft1` corresponds to physical shaft **S3** and `Shaft2` to physical shaft **S2**.

## Method in brief

**1. Loads.** Gear and belt forces are resolved into shear-force, bending-moment and torque diagrams for each shaft. Mean and alternating components (M<sub>a</sub>, M<sub>m</sub>, T<sub>a</sub>, T<sub>m</sub>) at each critical section are the solver input (Report, Section 5).

**2. Failure criterion.** For the fatigue-loaded shafts the modified Goodman relation with von Mises equivalent stresses is used:

$$
\frac{1}{n}=\frac{16}{\pi d^{3}}\left[\frac{\sqrt{4(K_f M_a)^2+3(K_{fs} T_a)^2}}{S_e}+\frac{\sqrt{4(K_f M_m)^2+3(K_{fs} T_m)^2}}{S_{ut}}\right]
$$

**3. Stress concentration.** Theoretical factors K<sub>t</sub> and K<sub>ts</sub> come from the polynomial correlations for shoulder fillets in Peterson's *Stress Concentration Factors*, specialised to t/r = 1. Notch sensitivity *q* uses the Neuber relation; K<sub>f</sub> = 1 + q (K<sub>t</sub> − 1).

**4. Numerical solution.**
- *Inner loop:* bisection on the Goodman residual. The size factor k<sub>b</sub>, endurance limit S<sub>e</sub> and the concentration factors all depend on the diameter being solved for, so the equation is not closed-form.
- *Outer loop:* each fillet radius is tied to the geometric step between neighbouring diameters; diameters and radii are updated (under-relaxed) until the diameter change falls below 0.1 mm.
- *Standardisation:* radii are rounded **down** to a preferred series, diameters **up** to a standard series, and the design is re-run with the standard radii.

Full derivations and assumptions are in the report (Sections 6–8).

## How to run

1. Open MATLAB and set the current folder to the project folder.
2. Run:

```matlab
Results = Gearbox_Shaft_Design('Shaft_Design_Input_Data.csv');
```

The command prints the formatted tables to the Command Window and returns the results as a MATLAB table. No files are written.

**Input format** (`Shaft_Design_Input_Data.csv`, one row per section, ordered left to right within each shaft group):

| Column | Meaning |
|---|---|
| `Section` | Section label (component pair) |
| `Shaft` | Computational shaft-group label |
| `r_initial_mm` | Initial fillet radius guess (mm) |
| `d_guess_mm` | Initial diameter guess (mm) |
| `Ma_Nm`, `Mm_Nm` | Alternating and mean bending moment (N·m) |
| `Ta_Nm`, `Tm_Nm` | Alternating and mean torque (N·m) |

Tested with MATLAB `[R2025b]`.

## Assumptions and limitations

- **Fillet correlations:** K<sub>t</sub> and K<sub>ts</sub> use the t/r = 1 specialisation of Peterson's polynomials, and are evaluated with the trial diameter of each row rather than with both shoulder diameters. The actual t/r can differ slightly after fillet radii are standardised. The factors are therefore an approximate implementation of the chart method.
- **Notch sensitivity:** the same *q* is used for bending and torsion.
- **Keys:** key lengths in the report (Section 9) were retained from the earlier calculation and were **not** recomputed for the final diameters.
- **Scope:** the solver covers the fatigue-controlled sections of S2 and S3 only; S1 is sized separately (static von Mises) in the report. Shaft deflection, critical speeds, bearing selection and finite-element validation are not included.
- **Status:** results have not been verified experimentally or with FEA.


## References

1. R. G. Budynas and J. K. Nisbett, *Shigley's Mechanical Engineering Design*, 10th ed., McGraw-Hill, 2015.
2. W. D. Pilkey, *Peterson's Stress Concentration Factors*, 2nd ed., Wiley, 1997 (Chapter 3, shoulder fillets).

## Author

**Mohammadmahdi Pourheydari** – B.Sc. Mechanical Engineering, K. N. Toosi University of Technology
Contact: `[mm.pourheydari@gmail.com]` 
