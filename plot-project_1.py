import numpy as np
import matplotlib.pyplot as plt

# read data
data = np.loadtxt('mass.dat')

r = data[:,0]
numerical_dm = data[:,1]
theoretical_dm = data[:,2]

# take log10 of radius and masses
log_r = np.log10(r)
log_numerical_dm = np.log10(numerical_dm)
log_theoretical_dm = np.log10(theoretical_dm)

# plot density profile
plt.fig, ax = plt.subplots(figsize=(7,7))
ax.plot(log_r, log_numerical_dm, 'k-', lw=1.5, label='Numerical DM')
ax.plot(log_r, log_theoretical_dm, 'r-', lw=1.5, label='Theoretical DM')
ax.set_xlabel("log(R) [kpc]")
ax.set_ylabel("log(M) [Msun]")
ax.set_title("Dark Matter Mass Profile (log-log)")
ax.legend()
plt.tight_layout()
plt.savefig("mass.png")