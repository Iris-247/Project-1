import numpy as np
import matplotlib as plt

# read data
data = np.loadtxt("Null")

r = data[:,1]
rho = data[:,2]
mass = data[:,3]


# plot density profile
plt.figure()
plt.loglog(r,rho,label="Density")
plt.xlabel("Radius")
plt.ylabel("Density")
plt.title("Dark Matter Density Profile")
plt.legend()
plt.tight_layout()
plt.savefig("density_profile.png")

# plot mass profile
plt.figure()
plt.loglog(r,mass,label="Mass")
plt.xlabel("Radius")
plt.ylabel("Mass")
plt.title("Dark Matter Distribution")
plt.legend()
plt.tight_layout()
plt.savefig("mass_distribution.png")