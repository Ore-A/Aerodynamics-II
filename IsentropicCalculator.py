import numpy as np
from scipy.optimize import fsolve

def isentropic_flow_calculator():
    # Specify the specific heat ratio (gamma)
    gamma = 1.4  # Replace with the appropriate value for your case

    # Display menu for user to choose variable
    print('Choose the variable to input:')
    print('1. P0/P (total to static pressure ratio)')
    print('2. M (Mach number)')
    print('3. T0/T (total to static temperature ratio)')
    print('4. rho0/rho (total to static density ratio)')
    print('5. A/A* (area ratio)')

    choice = int(input('Enter the corresponding number: '))

    # Get user input based on the chosen variable
    if choice == 1:
        P0_P = float(input('Enter P0/P: '))
    elif choice == 2:
        M = float(input('Enter M: '))
        P0_P = (1 + ((gamma - 1) / 2) * (M ** 2)) ** (gamma / (gamma - 1))
    elif choice == 3:
        T0_T = float(input('Enter T0/T: '))
        M = ((5 * (T0_T - 1)) ** 0.5)
        P0_P = (1 + ((gamma - 1) / 2) * (M ** 2)) ** (gamma / (gamma - 1))
    elif choice == 4:
        rho0_rho = float(input('Enter rho0/rho: '))
        T0_T = rho0_rho ** (2 / 5)
        M = ((5 * (T0_T - 1)) ** 0.5)
        P0_P = (1 + ((gamma - 1) / 2) * (M ** 2)) ** (gamma / (gamma - 1))
    elif choice == 5:
        A_Astar = float(input('Enter A/A*: '))
        
        # Define functions to find Mach numbers from area ratio
        def f_supersonic(M):
            return (1/M) * ((2/(gamma + 1)) * (1 + 0.5*(gamma - 1)*M**2))**((gamma + 1)/(2*(gamma - 1))) - A_Astar
        
        def f_subsonic(M):
            return (1/M) * ((2/(gamma + 1)) * (1 + 0.5*(gamma - 1)*M**2))**((gamma + 1)/(2*(gamma - 1))) - A_Astar

        # Solve for supersonic Mach number
        M_sup = fsolve(f_supersonic, 2)[0]
    
        # Solve for subsonic Mach number
        M_sub = fsolve(f_subsonic, 0.5)[0]
        
        P0_P_sup = (1 + ((gamma - 1) / 2) * (M_sup ** 2)) ** (gamma / (gamma - 1))  # Assuming supersonic
        P0_P_sub = (1 + ((gamma - 1) / 2) * (M_sub ** 2)) ** (gamma / (gamma - 1))  # Assuming subsonic
        T0_T_sup = (1 + ((gamma - 1) / 2) * (M_sup ** 2))
        T0_T_sub = (1 + ((gamma - 1) / 2) * (M_sub ** 2))
        rho0_rho_sup = (1 + ((gamma - 1) / 2) * (M_sup ** 2)) ** (1 / (gamma - 1))
        rho0_rho_sub = (1 + ((gamma - 1) / 2) * (M_sub ** 2)) ** (1 / (gamma - 1))

    else:
        print('Invalid choice. Please enter a number between 1 and 5.')
        return

    # Calculate other columns using isentropic flow formulas
    if choice != 5:
        if 'M' not in locals():
            M = ((2 / (gamma - 1)) * ((P0_P) ** ((gamma - 1) / gamma) - 1)) ** 0.5
        if 'P0_P' not in locals():
            P0_P = (1 + ((gamma - 1) / 2) * (M ** 2)) ** (gamma / (gamma - 1))
        if 'T0_T' not in locals():
            T0_T = (1 + ((gamma - 1) / 2) * (M ** 2))
        if 'rho0_rho' not in locals():
            rho0_rho = (1 + ((gamma - 1) / 2) * (M ** 2)) ** (1 / (gamma - 1))
        if 'A_Astar' not in locals():
            A_Astar = (1/(6**3))*((5+(M**2))**3/M)

    # Display the results
    print('Results:')
    if choice == 5:
        print(f'M_sub: {M_sub:.4f}  P0/P_sub: {P0_P_sub:.4f}  rho0/rho_sub: {rho0_rho_sub:.4f}  T0/T_sup: {T0_T_sub:.4f}  A/A*: {A_Astar:.4f}')
        print(f'M_sup: {M_sup:.4f}  P0/P_sup: {P0_P_sup:.4f}  rho0/rho_sup: {rho0_rho_sup:.4f}  T0/T_sup: {T0_T_sup:.4f}  A/A*: {A_Astar:.4f}')
    else:
        print(f'M: {M:.4f}')
        print(f'P0/P (total to static pressure ratio): {P0_P:.4f}')
        print(f'T0/T (total to static temperature ratio): {T0_T:.4f}')
        print(f'rho0/rho (total to static density ratio): {rho0_rho:.4f}')
        print(f'A/A* (area ratio): {A_Astar:.4f}')

# Call the function
isentropic_flow_calculator()




