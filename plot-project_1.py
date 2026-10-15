import numpy as np
import matplotlib.pyplot as plt

# read data
data = np.loadtxt('mass.dat')

r = data[:,0]
numerical_dm = data[:,1]
theoretical_dm = data[:,2]

# plot density profile
fig, map1 = plt.subplots(figsize=(7,7))

map1.loglog(r, numerical_dm, 'k-', lw=1.5, label='Numerical DM')
map1.loglog(r, theoretical_dm, 'r-', lw=1.5, label='Theoretical DM')

map1.set_xlabel(r"Radius $r$ [kpc]")
map1.set_ylabel(r"Enclosed Mass $M_{\rm DM}(<r)$ [$M_\odot$]")
map1.set_title("Dark Matter Mass Profile")
map1.legend()

fig.tight_layout()
fig.savefig("mass.png", dpi=300, bbox_inches="tight")
