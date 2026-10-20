# Physical-validation status — neutron-star structure project

## Completed and locally executed: analytic SLy approximation

The `physical_validation.py` script has been run successfully against the Read et al. (2009) SLy piecewise-polytropic approximation in `neutron_star_sly.py`.

- Full stellar-interior causality check for the analytic EOS: **PASS** on the reported causal branch. The EOS reaches `c_s = c` at approximately `2.008568e15 g/cm^3`; the reported model is kept just below that threshold, with maximum interior `c_s^2/c^2 = 0.9999999702`.
- Causal-boundary model: `M = 2.0462668 Msun`, `R = 9.982195 km`.
- Interpolated radius at `1.4 Msun`: `11.754450 km`.
- Compared with rounded official CompOSE RG(SLY4) scalar properties (`2.06 Msun`, `10.02 km`, `11.7 km`), differences are `-0.667%`, `-0.377%`, and `+0.465%`.
- Numerical sensitivity study: maximum observed differences are `8.01e-7%` in mass and `2.46e-5%` in radius, below the predeclared `0.1%` tolerance.

Outputs: `outputs/physical_validation_summary.txt`, `outputs/physical_validation_convergence.csv`, and `outputs/physical_validation_causal_mass_radius.csv`.

## Important boundary of the claim

These checks support an **approximation-level physical and numerical validation** of the causal branch of the analytic SLy piecewise-polytropic model. They do **not** prove that this approximation is identical to the tabulated SLy4 EOS, and the scalar comparison is a cross-model consistency check rather than an exact same-EOS reproduction.

## Still pending: table-driven CompOSE run

`neutron_star_compose.py` and `compose_convergence_check.py` implement a workflow for the official cold, beta-equilibrated CompOSE RG(SLY4) table. The actual `eos.thermo` and `eos.nb` bytes were not accessible to the local runtime, so this parser/solver has **not** been run end-to-end against the real table. To complete that separate validation:

1. Download the cold RG(SLY4) dataset from https://compose.obspm.fr/eos/134.
2. Extract `eos.thermo` and `eos.nb` to a `data/` directory beside the scripts.
3. Run `python neutron_star_compose.py --data-dir data`.
4. Run `python compose_convergence_check.py`.
5. Inspect the generated EOS causality data, low-pressure surface-extrapolation sensitivity, convergence summary, and mass-radius curve comparison against CompOSE's `eos.mr` reference.

Until these steps pass, do not describe the project as an exact, fully validated tabulated-SLy4 reproduction. CV wording should identify the completed model as a causal SLy piecewise-polytropic approximation with convergence testing and published benchmark comparison.
