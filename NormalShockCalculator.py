#-----Normal Shock Calculator-----#
import numpy as np
from scipy.optimize import fsolve

def shock_jump():
    """
    Calculates the downstream conditions after a normal shock wave using the jump formulas.

    Args:
        M1: Mach number upstream of the shock.

    Returns:
        M2: Mach number downstream of the shock.
        T2_T1: Temperature ratio downstream and upstream of the shock.
        P2_P1: Pressure ratio downstream and upstream of the shock.
        rho2_rho1: Density ratio downstream and upstream of the shock.
        P02_P01: Ratio of stagnation pressures downstream and upstream.
        P02_P1: Pitot Tube Formula
    """

    print('Choose the variable to input:')
    print('1. M₁ (Mach number upstream)')
    print('2. P₂/P₁ (pressure ratio)')
    print('3. ρ₂/ρ₁ (density ratio)')
    print('4. T₂/T₁ (temperature ratio)')
    print('5. P₀₂/P₀₁ (stagnation pressure ratio)')
    print('6. P₀₂/P₁ (Pitot Tube Formula)')
    print('7. M₂ (Mach number downstream)')

    choice = int(input('Enter the corresponding number: '))

    if choice == 1:
        M1 = float(input('Enter M₁: '))
    elif choice == 2:
        P2_P1 = float(input('Enter P₂/P₁: '))
        M1 = np.sqrt((6*P2_P1)+1/7)
    elif choice == 3:
        rho2_rho1 = float(input('Enter ρ₂/ρ₁: '))
        M1 = np.sqrt((-2*rho2_rho1)/(0.4*rho2_rho1-2.4))
        P2_P1 = ((7*(M1**2)-1)/6)
    elif choice == 4:
        T2_T1 = float(input('Enter T₂/T₁: '))
        def f(M1):
            return ((7*(M1**2)-1)*((5+(M1**2)))/(36*(M1**2))) - T2_T1
        M1 = fsolve(f, 2)[0]
        P2_P1 = ((7*(M1**2)-1)/6)
        rho2_rho1 = ((6*(M1**2))/(5 + M1**2))
    elif choice == 5:
        P02_P01 = float(input('Enter P₀₂/P₀₁: '))
        def f(M1):
            return (((2.4*M1**2)/((0.4*M1**2)+2))**(7/2))*((2.4/((2.8*M1**2)-0.4))**(5/2))- P02_P01
        M1 = fsolve(f, 1)[0]
        P2_P1 = ((7*(M1**2)-1)/6)
        rho2_rho1 = ((6*(M1**2))/(5*M1**2))
        T2_T1 = ((7*(M1**2)-1)*((5+(M1**2)))/(36*(M1**2)))
    elif choice == 6:
        P02_P1 = float(input('Enter P₀₂/P₁: '))
        def f(M1):
            return ((6**6)*(M1**7))/((5**(7/2))*((7*(M1**2)-1)**(5/2))) - P02_P1
        M1 = fsolve(f, 2)[0]
        P2_P1 = ((7*(M1**2)-1)/6)
        rho2_rho1 = ((6*(M1**2))/(5*M1**2))
        T2_T1 = ((7*(M1**2)-1)*((5+(M1**2)))/(36*(M1**2)))
        P02_P01 = ((P2_P1)*((1/T2_T1)**(7/2)))
    elif choice == 7:
        M2 = float(input('Enter M₂: '))
        M1 = np.sqrt((5+(M2**2))/(7*(M2**2)-1))
        P2_P1 = ((7*(M1**2)-1)/6)
        rho2_rho1 = ((6*(M1**2))/(5*M1**2))
        T2_T1 = ((7*(M1**2)-1)*((5+(M1**2)))/(36*(M1**2)))
        P02_P01 = ((P2_P1)*((1/T2_T1)**(7/2)))
        P02_P1 = ((6**6)*(M1**7))/((5**(7/2))*((7*(M1**2)-1)**(5/2)))
    else:
        print('Invalid choice. Please enter a number between 1 and 7.')
        return
    
    if 'M1' not in locals():
        M1 = np.sqrt((6*P2_P1)+1/7)
    if 'P2_P1' not in locals():
        P2_P1 = ((7*(M1**2)-1)/6)
    if 'rho2_rho1' not in locals():
        rho2_rho1 = ((7*(M1**2)-1)/5)
    if 'T2_T1' not in locals():
        T2_T1 = ((7*(M1**2)-1)*((5+(M1**2)))/(36*(M1**2)))
    if 'P02_P01' not in locals():
        P02_P01 = ((P2_P1)*((1/T2_T1)**(7/2)))
    if 'P02_P1' not in locals():
        P02_P1 = ((6**6)*(M1**7))/((5**(7/2))*((7*(M1**2)-1)**(5/2)))
    if 'M2' not in locals():
        M2 = np.sqrt((5+(M1**2))/(7*(M1**2)-1))
    
    # Print downstream conditions
    print("M₁ = {:.4f}".format(M1))
    print("P₂/P₁ = {:.4f}".format(P2_P1))
    print("ρ₂/ρ₁ = {:.4}".format(rho2_rho1))
    print("T₂/T₁ = {:.4f}".format(T2_T1))
    print("P₀₂/P₀₁ = {:.4f}".format(P02_P01))
    print("P₀₂/P₁ = {:.4f}".format(P02_P1)) 
    print("M₂ = {:.4f}".format(M2))
    
# Call the function
shock_jump()

