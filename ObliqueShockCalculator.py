
#-----Oblique Shock Calculator-----# 
import sys
import numpy as np
import math
from scipy.optimize import fsolve
sqrt = np.emath.sqrt

# User input for M1 and theta
M1_input = float(input("Enter the value of M1: "))
theta_input_degrees = float(input("Enter the value of theta in degrees: "))
theta_input = np.radians(theta_input_degrees)

# Optional input for T1 and P1
T1_input_str = input("Enter the value of T1 (hit Enter if not available): ")
P1_input_str = input("Enter the value of P1 (hit Enter if not available): ")

# Convert T1_input and P1_input to floats if they are not empty
T1_input = float(T1_input_str) if T1_input_str else None
P1_input = float(P1_input_str) if P1_input_str else None

#Theta Max

def theta_array(beta, M1):
    return beta - np.arctan(((5/(M1*M1*np.sin(beta)**2))+1)*np.tan(beta)/6)

beta_array = np.arange(0, (math.pi)/2, 0.1)
theta_values = theta_array(beta_array, M1_input)
theta_max = np.nanmax(theta_values)
theta_max= np.degrees(theta_max)

# beta-theta-M formula
initial_guess = np.radians(20)
gamma = 1.4

def f(beta, theta, M1):
    return theta - beta + np.arctan(((5/(M1*M1*np.sin(beta)**2))+1)*np.tan(beta)/6)


# Theta max check
mu = np.arcsin(1/M1_input)

if theta_input_degrees > theta_max:
    print("Shock is detached.")
    sys.exit()
else:
    # Solve for beta using fsolve
    beta = fsolve(f, initial_guess, args=(theta_input, M1_input))
    beta = np.degrees(beta)

    # Check if the obtained beta is physically meaningful
    if beta < 0 or beta > 90:
        print("Invalid solution for beta. Check input values.")
        sys.exit()

    # Oblique Shock Relations
    M1_N = M1_input * np.sin(np.radians(beta))
    M2_N = np.sqrt((1 + (((gamma - 1)/2) * (M1_N)**2)) / ((gamma * (M1_N**2)) - ((gamma-1)/2)))
    M2 = M2_N / (np.sin(np.radians(beta - theta_input_degrees)))
    P2_P1 = 1 + ((2*gamma)/(gamma+1)) * ((M1_N*M1_N) - 1)
    rho2_rho1 = ((gamma + 1) * M1_N * M1_N) / (2 + (gamma-1) * M1_N * M1_N)
    T2_T1 = P2_P1 * (1/rho2_rho1)
    P02_P01 = (P2_P1)*((1/(T2_T1))**(7/2))

    # Print results
    print(f"Mach angle (mu): {np.degrees(mu).item():.4f} degrees")
    print(f"Max wave angle (theta_max): {theta_max:.4f} degrees")
    print(f"Oblique shock angle (beta): {beta} degrees")
    print(f"M2: {M2.item():.4f}")
    print(f"T2/T1: {T2_T1.item():.4f}")
    print(f"P2/P1: {P2_P1.item():.4f}")
    print(f"P02/P01: {P02_P01.item():.4f}")
    
   # Optional: Print T2 and P2 if T1 and P1 are provided
    if T1_input is not None and P1_input is not None:
        T2 = T1_input * T2_T1
        P2 = P1_input * P2_P1
        print(f"T2: {T2.item():.4f}")
        print(f"P2: {P2.item():.4f}")
    





