"""
Neutron-star structure with the Tolman–Oppenheimer–Volkoff (TOV) equations.

Educational starter project. Uses a simple relativistic polytropic equation of state:
    P = K rho_0**Gamma
    epsilon = rho_0 c**2 + P/(Gamma - 1)

This is a deliberately simplified model, not a precision neutron-star EOS.
Units are cgs throughout the integration.
"""

from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Physical constants (cgs)
G = 6.67430e-8                 # cm^3 g^-1 s^-2
C = 2.99792458e10               # cm s^-1
M_SUN = 1.98847e33              # g
KM = 1.0e5                      # cm per km

# Model parameters (explore these in your project)
K = 1.0e5                       # polytropic constant; cgs for Gamma=2
GAMMA = 2.0
CENTRAL_PRESSURES = np.logspace(33, 36, 40)  # dyn cm^-2


def integrate_star(p_c, k=K, gamma=GAMMA):
    """Integrate one stellar model; return (mass in solar masses, radius in km)."""
    if p_c <= 0:
        raise ValueError("Central pressure must be positive.")

    rho_c = (p_c / k) ** (1.0 / gamma)
    r0 = 1.0  # cm; small nonzero starting radius avoids the r=0 singularity
    eps_c = rho_c * C**2 + p_c / (gamma - 1.0)
    m0 = (4.0 / 3.0) * np.pi * r0**3 * eps_c / C**2

    # Leading-order central expansion for pressure at r0
    denom0 = 1.0 - 2.0 * G * m0 / (r0 * C**2)
    dpdr0 = (
        -G * (eps_c + p_c)
        * (m0 + 4.0 * np.pi * r0**3 * p_c / C**2)
        / (C**2 * r0**2 * denom0)
    )
    p0 = p_c + 0.5 * dpdr0 * r0**2

    def rhs(r, y):
        m, p = y
        # The event stops the integration at the surface (P=0); guard solver
        # trial steps from evaluating a fractional power of a negative pressure.
        p_eval = max(p, 1.0)
        rho0 = (p_eval / k) ** (1.0 / gamma)
        epsilon = rho0 * C**2 + p_eval / (gamma - 1.0)
        compactness_factor = 1.0 - 2.0 * G * m / (r * C**2)
        compactness_factor = max(compactness_factor, 1e-12)

        dm_dr = 4.0 * np.pi * r**2 * epsilon / C**2
        dp_dr = (
            -G * (epsilon + p_eval)
            * (m + 4.0 * np.pi * r**3 * p_eval / C**2)
            / (C**2 * r**2 * compactness_factor)
        )
        return [dm_dr, dp_dr]

    def surface_event(r, y):
        return y[1]

    surface_event.terminal = True
    surface_event.direction = -1

    sol = solve_ivp(
        rhs,
        (r0, 5.0e7),                 # maximum radius: 500 km
        [m0, p0],
        events=surface_event,
        rtol=1e-7,
        atol=[1e15, 1e15],
        max_step=1.0e4,
    )

    if not sol.success or len(sol.t_events[0]) == 0:
        return None

    radius_cm = sol.t_events[0][0]
    mass_g = sol.y_events[0][0][0]
    if not np.isfinite(mass_g) or not np.isfinite(radius_cm):
        return None

    return {
        "central_pressure_dyn_cm2": float(p_c),
        "central_density_g_cm3": float(rho_c),
        "mass_solar": float(mass_g / M_SUN),
        "radius_km": float(radius_cm / KM),
    }


def main():
    out = Path(__file__).resolve().parent / "outputs"
    out.mkdir(exist_ok=True)

    rows = []
    for pc in CENTRAL_PRESSURES:
        result = integrate_star(pc)
        if result is not None:
            rows.append(result)

    if len(rows) < 5:
        raise RuntimeError("Too few successful integrations. Check parameters and solver.")

    csv_path = out / "mass_radius_results.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    masses = np.array([row["mass_solar"] for row in rows])
    radii = np.array([row["radius_km"] for row in rows])
    pressures = np.array([row["central_pressure_dyn_cm2"] for row in rows])

    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.plot(radii, masses, marker="o", markersize=3, linewidth=1.6)
    imax = int(np.argmax(masses))
    ax.scatter([radii[imax]], [masses[imax]], marker="*", s=130, zorder=3,
               label=f"Maximum in sampled grid: {masses[imax]:.3f} M☉")
    ax.set_xlabel("Radius (km)")
    ax.set_ylabel("Gravitational mass (solar masses)")
    ax.set_title("Neutron-star mass–radius curve (toy polytropic EOS)")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    plot_path = out / "mass_radius_curve.png"
    fig.savefig(plot_path, dpi=180)
    plt.close(fig)

    summary_path = out / "run_summary.txt"
    with summary_path.open("w", encoding="utf-8") as f:
        f.write("Educational TOV solver run summary\n")
        f.write(f"EOS: P = K rho_0^Gamma, K={K:g}, Gamma={GAMMA:g}\n")
        f.write(f"Successful models: {len(rows)} of {len(CENTRAL_PRESSURES)}\n")
        f.write(f"Maximum mass in sampled grid: {masses[imax]:.4f} solar masses\n")
        f.write(f"Radius at sampled maximum: {radii[imax]:.3f} km\n")
        f.write(f"Central pressure at sampled maximum: {pressures[imax]:.4e} dyn/cm^2\n")
        f.write("\nIMPORTANT: these are model-dependent toy-EOS results, not observational predictions.\n")

    print(f"Saved {csv_path}")
    print(f"Saved {plot_path}")
    print(f"Maximum mass in sampled grid: {masses[imax]:.4f} M_sun at R={radii[imax]:.3f} km")
    print("Interpret this only as a result for the chosen simplified polytropic EOS.")


if __name__ == "__main__":
    main()
