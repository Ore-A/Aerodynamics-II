#-----Prandtl-Meyer Calculator-----
import numpy as np
from scipy.optimize import fsolve

# Constants
M1 = float(input("Enter Mach number: "))
theta = float(input("Enter θ: "))   
gamma = 1.4  # Specific heat ratio for air

nu_M1 = np.degrees((np.sqrt((gamma+1)/(gamma-1))*np.arctan(np.sqrt(((gamma-1)/(gamma+1))*((M1**2)-1))))-np.arctan(np.sqrt((M1**2)-1)))

max_theta = np.degrees(np.arcsin(1 / M1) + np.pi / 2)

def prandtl_meyer(M1, theta, gamma):
    # Converting angle from degrees to radians
    theta_rad = np.radians(theta)

    # Solve for nu(M2) using Newton's method
    nu_M2 = theta + nu_M1

    # Mach number M2
    def f(M2, gamma, nu_M2):
        return -nu_M2 + np.degrees((np.sqrt((gamma+1)/(gamma-1))*np.arctan(np.sqrt(((gamma-1)/(gamma+1))*((M2**2)-1))))-np.arctan(np.sqrt((M2**2)-1)))
    
    M2 = fsolve(f, 1.2, args=(gamma, nu_M2))[0]  # Accessing the first element of the array returned by fsolve

    # Mach angles mu1 and mu2
    mu1 = np.degrees(np.arcsin(1 / M1))
    mu2 = np.degrees(np.arcsin(1 / M2))

    # Isentropic ratios
    T0_T1 = (1 + ((gamma - 1) / 2) * (M1 ** 2))
    P0_P1 = (1 + (gamma - 1) / 2 * (M1**2)) ** (gamma / (gamma - 1))
    rho0_rho1 = (1 + (gamma - 1) / 2 * (M1**2)) ** (1 / (gamma - 1))

    T0_T2 = (1 + ((gamma - 1) / 2) * (M2 ** 2))
    P0_P2 = (1 + (gamma - 1) / 2 * (M2**2)) ** (gamma / (gamma - 1))
    rho0_rho2 = (1 + (gamma - 1) / 2 * (M2**2)) ** (1 / (gamma - 1))

    # Print results
    print("ν(M1) =", round(nu_M1, 3), "°, ν(M2) =", round(nu_M2, 3), "°, Max θ =", round(max_theta, 3), "°")
    print("\nM2 =", round(M2, 3))
    print("Mach angles: μ1 =", round(mu1, 2), "°", "μ2 =", round(mu2, 2), "°")
    print("\nIsentropic formulas:")
    print("M1 =", round(M1, 3))
    print("T0/T1 =", round(T0_T1, 3))
    print("P0/P1 =", round(P0_P1, 3))
    print(" ρ0/ρ1 =", round(rho0_rho1, 3))
    print("\nM2 =", round(M2, 3))
    print("T0/T2 =", round(T0_T2, 3))
    print("P0/P2 =", round(P0_P2, 3))
    print("ρ0/ρ2 =", round(rho0_rho2, 3))


# Call function
prandtl_meyer(M1, theta, gamma)

