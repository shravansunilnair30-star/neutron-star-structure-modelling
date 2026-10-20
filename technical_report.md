# Structure of Neutron Stars Using the TOV Equations
**Computational modelling project | Python, NumPy, SciPy, Matplotlib**

## Abstract
This project implements a numerical solution of the Tolman–Oppenheimer–Volkoff (TOV) equations for a static, spherically symmetric star. A simplified polytropic equation of state (EOS), \(P=K\rho_0^\Gamma\), is used to relate pressure to rest-mass density, with total energy density approximated by \(\epsilon=\rho_0c^2+P/(\Gamma-1)\). The code integrates from a small radius near the centre until pressure reaches zero, recording the stellar mass and radius for a sequence of central pressures. It generates a mass–radius curve and tabulated results. Numerical sensitivity is checked by changing the ODE solver tolerance and maximum radial step. This is an educational baseline model: its outputs depend strongly on the chosen EOS and should not be interpreted as realistic observational predictions without further EOS validation.

## 1. Introduction
Neutron stars are compact stellar remnants whose internal densities make relativistic hydrostatic equilibrium important. The mass–radius relationship links a model's equation of state to observable global properties. The objective here is to implement the TOV structure equations, produce a mass–radius curve, and establish a reproducible workflow for numerical verification.

## 2. Governing equations
The TOV equations in cgs units are
\[
\frac{dm}{dr}=4\pi r^2\frac{\epsilon}{c^2},
\]
\[
\frac{dP}{dr}=-\frac{G(\epsilon+P)\left(m+4\pi r^3P/c^2\right)}
{c^2r^2\left(1-2Gm/(rc^2)\right)}.
\]
Here \(r\) is radial distance, \(m(r)\) is enclosed gravitational mass, \(P(r)\) is pressure, \(\epsilon(r)\) is energy density, \(G\) is the gravitational constant, and \(c\) is the speed of light. The model assumes a static, spherical, non-rotating star and isotropic pressure.

## 3. Equation of state and assumptions
The starter model uses
\[
P=K\rho_0^\Gamma,\qquad \epsilon=\rho_0c^2+\frac{P}{\Gamma-1},
\]
with \(\Gamma=2\) and \(K=10^5\) in the code's cgs convention. This is a simplified polytropic prescription, not a tabulated nuclear-matter EOS. Therefore its mass and radius values are useful for demonstrating the method, not for claiming agreement with observed neutron stars.

## 4. Numerical method
1. Select a central pressure from a logarithmically spaced grid.
2. Infer the central density from the EOS.
3. Set initial mass and pressure near \(r=0\), avoiding direct evaluation at the coordinate singularity.
4. Integrate the coupled ODEs with `scipy.integrate.solve_ivp`.
5. Stop when pressure crosses zero; record that radius and the enclosed mass.
6. Repeat across central pressures and plot mass against radius.

The script outputs `mass_radius_results.csv`, `mass_radius_curve.png`, and `run_summary.txt`.

## 5. Baseline results
The initial run generated a mass–radius curve for the selected toy EOS. In the sampled central-pressure grid, the reported maximum was approximately **1.357 solar masses** at a radius of approximately **9.314 km**. These are model-specific results, not a validated estimate of the true maximum neutron-star mass. Refer to `outputs/mass_radius_results.csv` and `outputs/run_summary.txt` for the run's values.

## 6. Numerical convergence / sensitivity testing
The companion script `convergence_check.py` repeats the central-pressure sequence with different ODE relative tolerances and maximum radial steps. It writes `outputs/convergence_test.csv` and `outputs/convergence_test_summary.txt`.

Report the comparison using the generated CSV. For each configuration, compare the maximum sampled mass and its corresponding radius with the baseline. Small changes indicate that these particular solver settings are numerically stable for this implementation; they do not prove the EOS is physically realistic. Also note that the maximum of a finite pressure grid can shift when the grid density changes. A stronger extension is to increase the number of central-pressure samples (for example 40, 80, 160) and compare the sampled peak.

## 7. Literature-calibrated EOS extension and quantitative comparison

A second solver (`neutron_star_sly.py`) uses the high-density piecewise-polytrope parameters associated with the SLy fit in Read et al. (2009): \(\log_{10} p_1=34.384\) in cgs units and \(\Gamma_1=3.005\), \(\Gamma_2=2.988\), and \(\Gamma_3=2.851\). A low-density polytropic segment is matched at \(10^{14.3}\,\mathrm{g\,cm^{-3}}\). Energy-density offsets are used to maintain continuity at the segment boundaries. This is an approximation to the fitted EOS and is not the original tabulated SLy4 microphysical EOS.

A refined sweep was run at 220 logarithmically spaced central pressures using DOP853, relative tolerance \(10^{-9}\), and maximum radial step \(10^4\) cm. All 220 requested integrations succeeded. The unrestricted sampled peak was \(2.046263\,M_\odot\) at \(R=9.97695\) km. Linear interpolation between the two low-central-pressure models bracketing \(1.4\,M_\odot\) gives \(R_{1.4}=11.7544\) km. Compared with the often-quoted Read et al. fit values of approximately \(M_{\max}=2.06\,M_\odot\) and \(R_{1.4}=11.4\) km, the peak mass is 0.67% lower, while the radius is 3.11% higher. The comparison is quantitative, but the discrepancy means the model has not passed a strict validation gate.

### Physical-validity check

For a barotropic EOS, the local sound speed obeys \(c_s^2/c^2=dP/d\epsilon\). For each polytropic segment in this implementation, this derivative is evaluated as \(\Gamma P/(\epsilon+P)\). The check finds that the central sound speed exceeds \(c\) in the high-central-pressure region containing the unrestricted sampled mass peak. Consequently, that peak cannot be treated as a physically valid maximum-mass prediction. The approximate EOS is therefore not suitable for claiming a validated SLy4 result. A central-only check is also not a complete proof of causality throughout each stellar interior; the full EOS and each star's density range must be checked.

The detailed sweep is saved in `outputs/sly_refined_mass_radius_results.csv`, the calculated curve in `outputs/sly_mass_radius_curve.png`, and numerical results plus the validation verdict in `outputs/sly_refined_results_and_validation.txt`.

## 8. Physical validation of the causal SLy approximation

A new validation script, `physical_validation.py`, checks the sound-speed bound over the complete density interval sampled by each star in the causal branch of the analytic piecewise-polytropic EOS. For a segment with
\[
\epsilon(\rho)=(1+A_i)\rho c^2+\frac{P}{\Gamma_i-1},\qquad P=K_i\rho^{\Gamma_i},
\]
the derivative is evaluated from the analytic derivatives, \(c_s^2/c^2=(dP/d\rho)/(d\epsilon/d\rho)\). The EOS is monotone within each segment; segment boundaries are checked from both sides. The highest density at which the approximation reaches \(c_s=c\) is approximately \(2.008568\times10^{15}\,\mathrm{g\,cm^{-3}}\), corresponding to a central pressure of \(1.392931\times10^{36}\,\mathrm{dyn\,cm^{-2}}\). The reported causal-boundary model uses a fractional safety margin of \(10^{-7}\), and its maximum interior \(c_s^2/c^2\) is 0.9999999702.

At this causal limit, the model gives \(M=2.0462668\,M_\odot\) and \(R=9.982195\) km. Interpolation along the causal sequence gives \(R_{1.4}=11.754450\) km. The official CompOSE cold RG(SLY4) page lists rounded values \(M_{\max}\simeq2.06\,M_\odot\), \(R_{M_{\max}}\simeq10.02\) km, and \(R_{1.4}\simeq11.7\) km. Relative to these rounded targets, the approximation differs by -0.667% in maximum mass, -0.377% in radius at maximum mass, and +0.465% in \(R_{1.4}\). All three are within the project's declared 2% scalar-comparison tolerance. These numbers are cross-model consistency checks, not an exact same-EOS reproduction: the analytic fit and crust prescription are not identical to the CompOSE tabulated EOS.

The convergence study repeats the causal-boundary integration with four combinations of solver tolerance and maximum radial step. The largest changes relative to the baseline are approximately \(8.01\times10^{-7}\)% in mass and \(2.46\times10^{-5}\)% in radius, below the predeclared 0.1% tolerance. This supports numerical stability under the tested integration settings; it does not eliminate EOS model-form error.

The outputs are `outputs/physical_validation_summary.txt`, `outputs/physical_validation_convergence.csv`, and `outputs/physical_validation_causal_mass_radius.csv`. Reproduce them with `python physical_validation.py`.

### What this does and does not establish

**Passed for the analytic SLy approximation:** the tested causal branch satisfies the sound-speed bound across each stellar interior; the maximum-mass/radius values are stable under the tested integration settings; and its three global properties agree with the rounded CompOSE RG(SLY4) scalar values to within 0.7%.

**Not yet completed:** the solver `neutron_star_compose.py` has not been run against the actual tabulated `eos.thermo` and `eos.nb` data files. Therefore this is not a full validation of the original tabulated SLy4 EOS implementation, and the close scalar agreement must not be presented as proof that the two EOS models are identical. A full table-driven reproduction still requires running the official table, inspecting the complete causality curve and interpolation/surface treatment, testing convergence, and comparing the computed mass-radius sequence with CompOSE's provided `eos.mr` reference curve.

## 9. Conclusion

The project now has a physically screened and numerically converged causal branch for a literature-calibrated SLy piecewise-polytropic approximation. The causal-boundary solution is \(M=2.0462668\,M_\odot\), \(R=9.982195\) km, with \(R_{1.4}=11.754450\) km. Compared with the rounded values published on the official CompOSE RG(SLY4) page, the relative differences are -0.667%, -0.377%, and +0.465%, respectively. The EOS causality test and solver sensitivity checks pass for the range reported here.

The defensible conclusion is that **the analytic SLy approximation passes the stated approximation-level causality, convergence, and scalar-reference checks**. This does not establish an exact reproduction of the original tabulated SLy4 EOS. The official tabulated-EOS solver is implemented but remains unexecuted against the actual data files in this environment; its full table-driven validation is therefore still pending. On a CV, describe the completed work as *TOV modelling with a causal SLy piecewise-polytropic approximation, numerical convergence testing, and comparison against published neutron-star benchmarks*—not as a fully validated tabulated SLy4 solver.

## Reproducibility

Install dependencies with `pip install -r requirements.txt`. Run `python neutron_star_tov.py` for the educational toy EOS, `python neutron_star_sly.py` for the SLy piecewise-polytrope approximation, and `python physical_validation.py` for the physical-causality, benchmark, and convergence checks. Run `python neutron_star_compose.py --data-dir data` and `python compose_convergence_check.py` only after obtaining the official cold CompOSE `eos.thermo` and `eos.nb` files. All generated results and plots are saved in `outputs/`.

## 10. Tabulated CompOSE EOS workflow (implemented; execution pending)

`neutron_star_compose.py` reads the standard CompOSE `eos.thermo` and `eos.nb` files for the cold, beta-equilibrated RG(SLY4) table. It reconstructs pressure and energy density from the documented per-baryon quantities, converts MeV fm\(^{-3}\) to cgs units, interpolates \(\epsilon(P)\), and calculates \(c_s^2/c^2=dP/d\epsilon\). `compose_convergence_check.py` supplies the planned solver sensitivity test. The official dataset and its documentation are available at [CompOSE RG(SLY4)](https://compose.obspm.fr/eos/134) and in the [CompOSE reference manual](https://doi.org/10.1140/epja/s10050-022-00847-y).

The local execution environment could not retrieve the raw CompOSE table bytes, so this separate table-driven implementation has not been run against the actual files. Do not imply that its parser or outputs have passed end-to-end testing until those files have been used and the generated causality, convergence, and mass-radius comparison results have been reviewed.
