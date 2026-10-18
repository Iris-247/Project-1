import numpy as np
import matplotlib.pyplot as plt

# DM
mass = np.loadtxt("mass.dat")
r1 = mass[:, 0]
numerical_dm = mass[:, 1]
theoretical_dm = mass[:, 2]

mask_dm = (
    np.isfinite(mass).all(axis=1)
    & (r1 > 0)
    & (numerical_dm > 0)
    & (theoretical_dm > 0)
)

fig1, ax1 = plt.subplots(figsize=(7, 7))
ax1.loglog(r1[mask_dm], numerical_dm[mask_dm], "k-", lw=1.5,
           label="Numerical DM")
ax1.loglog(r1[mask_dm], theoretical_dm[mask_dm], "r-", lw=1.5,
           label="Analytical DM")

ax1.set_xlabel(r"Radius $r$ [kpc]")
ax1.set_ylabel(r"Enclosed mass $M_{\mathrm{DM}}(<r)$ [$M_\odot$]")
ax1.set_title("Dark Matter Enclosed Mass Profile")
ax1.legend()
ax1.grid(True, which="both", linestyle=":")

fig1.tight_layout()
fig1.savefig("mass.png", dpi=300, bbox_inches="tight")

# Gas
gas = np.loadtxt("gas.dat")
r2 = gas[:, 0]
numerical_gas = gas[:, 1]
theoretical_gas = gas[:, 2]

mask_gas = (
    np.isfinite(gas).all(axis=1)
    & (r2 > 0)
    & (numerical_gas > 0)
    & (theoretical_gas > 0)
)

fig2, ax2 = plt.subplots(figsize=(7, 7))
ax2.loglog(r2[mask_gas], numerical_gas[mask_gas], "k-", lw=1.5,
           label="Numerical HSE")
ax2.loglog(r2[mask_gas], theoretical_gas[mask_gas], "r-", lw=1.5,
           label="Analytical HSE")

ax2.set_xlabel(r"Radius $r$ [kpc]")
ax2.set_ylabel(r"Gas density $\rho_g$ [g cm$^{-3}$]")
ax2.set_title("Isothermal Gas Density Profile")
ax2.legend()
ax2.grid(True, which="both", linestyle=":")

fig2.tight_layout()
fig2.savefig("gas.png", dpi=300, bbox_inches="tight")

plt.show()