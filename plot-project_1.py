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
rho_dm = gas[:, 3]

mask_gas = (
    np.isfinite(gas).all(axis=1)
    & (r2 > 0)
    & (numerical_gas > 0)
    & (theoretical_gas > 0)
)

mask_dm = (
    np.isfinite(r2) 
    & np.isfinite(rho_dm) 
    & (r2 > 0) 
    & (rho_dm > 0)
)

fig2, ax2 = plt.subplots(figsize=(7, 7))
ax2.loglog(r2[mask_gas], numerical_gas[mask_gas], "k-", lw=1.5, label="Numerical HSE")
ax2.loglog(r2[mask_gas], theoretical_gas[mask_gas], "r-", lw=1.5, label="Analytical HSE")
ax2.loglog(r2[mask_dm], rho_dm[mask_dm], "b--", lw=1.5, label="DM density")

ax2.set_xlabel(r"Radius $r$ [kpc]")
ax2.set_ylabel(r"Gas density $\rho_g$ [g cm$^{-3}$]")
ax2.set_title("Isothermal Gas Density Profile")
ax2.legend()
ax2.grid(True, which="both", linestyle=":")

fig2.tight_layout()
fig2.savefig("gas.png", dpi=300, bbox_inches="tight")

# Gas_BCG_comparation
gas_no_bcg = np.loadtxt("gas.dat")
gas_bcg = np.loadtxt("gas_BCG.dat")

r_no_bcg = gas_no_bcg[:, 0]
rho_no_bcg = gas_no_bcg[:, 1]

r_bcg = gas_bcg[:, 0]
rho_bcg = gas_bcg[:, 1]

rho_dm = gas_bcg[:, 2]

mask_no_bcg = (
    np.isfinite(r_no_bcg)
    & np.isfinite(rho_no_bcg)
    & (r_no_bcg > 0)
    & (rho_no_bcg > 0)
)

mask_bcg = (
    np.isfinite(r_bcg)
    & np.isfinite(rho_bcg)
    & (r_bcg > 0)
    & (rho_bcg > 0)
)

mask_dm = (
    np.isfinite(r_bcg) 
    & np.isfinite(rho_dm) 
    & (r2 > 0) 
    & (rho_dm > 0)
)


fig3, ax3 = plt.subplots(figsize=(7, 7))
ax3.loglog(r_no_bcg[mask_no_bcg], rho_no_bcg[mask_no_bcg], "k-", lw=1.5, label="Without BCG")
ax3.loglog(r_bcg[mask_bcg], rho_bcg[mask_bcg], "r-", lw=1.5, label="BCG")
ax3.loglog(r_bcg[mask_dm], rho_dm[mask_dm], "b--", lw=1.5, label="DM density")

ax3.set_xlabel(r"Radius $r$ [kpc]")
ax3.set_ylabel(r"Gas density $\rho_g$ [g cm$^{-3}$]")
ax3.set_title("Isothermal Gas Density Profile with BCG and non-BCG")
ax3.legend()
ax3.grid(True, which="both", linestyle=":")

fig3.tight_layout()
fig3.savefig("gas_BCG_comparison.png", dpi=300, bbox_inches="tight")
