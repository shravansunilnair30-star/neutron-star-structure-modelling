# Structure of Neutron Stars using the TOV Equations

**Author:** [Your name]  
**Date:** [Date]  
**Software:** Python, NumPy, SciPy, Matplotlib

## 1. Abstract
Write 100–150 words summarising the objective, method, main numerical result, and limitations. Insert your actual result from `outputs/run_summary.txt`. Do not claim that a toy EOS represents observed neutron stars.

## 2. Introduction
- What is a neutron star?
- Why does general relativity matter for hydrostatic equilibrium?
- What is the mass–radius relationship, and why is it useful?

## 3. Governing equations
Present the two TOV equations and define every variable and unit. Explain the assumptions: spherical symmetry, static equilibrium, and no rotation or magnetic-field effects in this model.

## 4. Equation of state
State the polytropic EOS:
- \(P=K\rho_0^\Gamma\)
- \(\epsilon=\rho_0c^2+P/(\Gamma-1)\)

Explain what the parameters \(K\) and \(\Gamma\) control. State clearly that this is a simplified educational model.

## 5. Numerical method
Describe:
- The central pressure grid.
- The initial conditions near the centre.
- The ODE integrator (`scipy.integrate.solve_ivp`).
- The surface event at \(P=0\).
- The method used to extract mass and radius.

## 6. Results
Include:
- The mass–radius graph.
- A table of selected central pressures, central densities, masses, and radii.
- The maximum mass in the sampled grid and its corresponding radius.

Use the phrase **“maximum in the sampled grid”**, not “the true maximum mass,” unless you have demonstrated convergence and explored the complete stable branch.

## 7. Verification and convergence
Run at least two checks:
1. Increase the number of central-pressure samples and see whether the curve's shape and sampled maximum change.
2. Tighten solver tolerances and compare selected mass/radius values.

Record numerical differences and discuss whether they are acceptably small.

## 8. Discussion
Explain the broad trends in the computed curve. Discuss why a maximum may appear and why the physical interpretation depends on the EOS. Avoid comparing directly to observations without using a realistic EOS and appropriate reference data.

## 9. Limitations
Discuss the simplified EOS, non-rotating assumption, idealised model, numerical resolution, and the fact that this is not a full nuclear-matter calculation.

## 10. Conclusion
Summarise what you implemented, what you verified, what you learned, and the next improvement you would make.

## References
Add the textbooks or papers you actually consulted. Verify bibliographic details yourself; do not cite a source you have not read.
