# Structure of Neutron Stars — TOV + Python

## Aim
Numerically integrate the Tolman–Oppenheimer–Volkoff (TOV) equations for a spherical, non-rotating star and plot its mass–radius relationship under a simplified polytropic equation of state (EOS).

## What is included
- `neutron_star_tov.py` — educational single-polytrope TOV solver
- `neutron_star_sly.py` — SLy literature-calibrated piecewise-polytropic TOV solver
- `sly_convergence_check.py` — sensitivity check for the SLy model
- `requirements.txt` — Python dependencies
- `report_template.md` — report structure and questions to answer
- `outputs/` — generated CSV, graph, and run summary after running the script

## Run it
Use Python 3.10+ in a terminal:

```bash
python -m venv .venv
# Windows:
.venv\Scripts\activate
# macOS/Linux:
source .venv/bin/activate
pip install -r requirements.txt
python neutron_star_tov.py
```

The scripts write results into the `outputs` folder. To run the literature-calibrated extension, use `python neutron_star_sly.py`; to test its numerical sensitivity, use `python sly_convergence_check.py`.

## Model and equations
The TOV equations are

\[
\frac{dm}{dr} = 4\pi r^2\frac{\epsilon}{c^2}
\]

\[
\frac{dP}{dr} =
-\frac{G(\epsilon+P)\left(m+4\pi r^3P/c^2\right)}
{c^2r^2\left(1-2Gm/(rc^2)\right)}
\]

where `m(r)` is enclosed gravitational mass, `P(r)` is pressure, and `epsilon` is total energy density. The chosen EOS is

\[
P=K\rho_0^\Gamma,\qquad
\epsilon=\rho_0c^2+\frac{P}{\Gamma-1}.
\]

The stellar surface is defined by `P = 0`; the mass and radius at that point are recorded.

## Models and scientific limitations
The original `neutron_star_tov.py` uses a deliberately simple single polytrope for teaching. The extension `neutron_star_sly.py` uses published SLy piecewise-polytrope parameters from Read et al. (2009) and a matched low-density crust approximation. `physical_validation.py` checks causality over the star interior, solver sensitivity, and scalar comparison with the official CompOSE RG(SLY4) properties. It is still an analytic approximation, not the original tabulated SLy4 EOS; describe the scope explicitly.

## Suggested project workflow
1. Run the code unchanged and inspect the CSV and mass–radius plot.
2. Explain every term in the TOV equations in your own words.
3. Change the number/range of central pressures and confirm the curve converges.
4. Change `Gamma` and `K` one at a time; explain the scaling and limitations.
5. Add a convergence study (for example, compare solver tolerances).
6. Compare the SLy piecewise-polytropic output against Read et al. (2009), noting the approximate crust model and finite grid.
7. Write up assumptions, validation, limitations, and reproducibility.

Do not list this as completed research until you have run it, checked it, and can explain the code and results yourself.


## Physical validation results (read before using results)

Run `python physical_validation.py` to reproduce the completed approximation-level checks. The causal branch of the Read et al. SLy piecewise-polytropic approximation gives **M = 2.0462668 solar masses**, **R = 9.982195 km** at the causal boundary, and **R(1.4 solar masses) = 11.754450 km**. Relative to the rounded properties published on the official CompOSE RG(SLY4) page (2.06 solar masses, 10.02 km, 11.7 km), differences are **-0.667%, -0.377%, and +0.465%**. The analytic EOS is checked for sound-speed causality across the stellar interior, and the tested numerical changes are below 0.1%.

**Scope caveat:** this is a validated *analytic approximation-level study*, not a full execution of the official tabulated SLy4 EOS. The table-driven code `neutron_star_compose.py` still requires the actual CompOSE `eos.thermo` and `eos.nb` files and has not passed an end-to-end run in this environment. Do not claim an exact tabulated-SLy4 reproduction until that stage is run and its full causality, convergence, and mass-radius curve checks pass.

## Tabulated-EOS physical-validation stage (new)

`neutron_star_compose.py` and `compose_convergence_check.py` provide the next validation workflow using the official CompOSE RG(SLY4) cold beta-equilibrium EOS table. Source: https://compose.obspm.fr/eos/134 . The CompOSE page lists approximate reference properties `Mmax = 2.06 Msun`, `R_Mmax = 10.02 km`, and `R_1.4 = 11.7 km`.

To run this stage, download the RG(SLY4) table from the CompOSE page and extract `eos.thermo` and `eos.nb` into a folder named `data/` beside the scripts. Then run:

```bash
python neutron_star_compose.py --data-dir data
python compose_convergence_check.py
```

The solver writes a mass-radius CSV/plot, a table-wide sound-speed causality CSV, a benchmark comparison summary, and a convergence summary under `outputs/`. The script refuses to claim success merely because integrations complete; inspect its causality and benchmark verdicts. The low-pressure extension below the finite-density EOS table is explicitly approximate and must be sensitivity-tested.

**Status caveat:** the CompOSE table bytes were not available to the current execution environment during preparation, so this workflow has been added and syntax-checked but has not been run against the actual tabulated EOS here. Therefore the project must not yet be described as fully physically validated. Once the EOS files are supplied and the scripts run, review the generated diagnostics before using any numerical result on a CV.


## Reproduce the physical validation

```bash
python physical_validation.py
```

Review `outputs/physical_validation_summary.txt` and `outputs/physical_validation_convergence.csv`. The full table-driven route is separate and requires the official cold EOS files: `python neutron_star_compose.py --data-dir data`, followed by `python compose_convergence_check.py`.
