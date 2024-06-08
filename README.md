# Aerodynamics-II

This repository contains Python scripts for various aerodynamics calculations. These scripts were used in a class to perform specific calculations related to isentropic processes, normal shocks, oblique shocks, and Prandtl-Meyer expansion waves.

## Contents

1. [IsentropicCalculator.py](#IsentropicCalculatorpy)
2. [NormalShockCalculator.py](#NormalShockCalculatorpy)
3. [ObliqueShockCalculator.py](#ObliqueShockCalculatorpy)
4. [Prandtl-MeyerCalculator.py](#Prandtl-MeyerCalculatorpy)

### IsentropicCalculator.py

This script calculates various properties of isentropic flows. Isentropic processes are idealized thermodynamic processes in which the entropy of the fluid remains constant.

#### Features:
- Calculation of temperature, pressure, and density ratios.
- Mach number computations.
- Validation of isentropic relations.

#### Usage:
```python
python IsentropicCalculator.py
```
### NormalShockCalculator.py

This script calculates properties of normal shocks. Normal shocks are abrupt, nearly vertical disruptions in the flow, typically occurring in supersonic flows.

#### Features:
- Calculation of post-shock properties such as pressure, temperature, and density.
- Determination of Mach number before and after the shock.

#### Usage:
```python
python NormalShockCalculator.py
```

### ObliqueShockCalculator.py

This script calculates properties of oblique shocks, which occur when the shock wave forms at an angle to the flow direction.

#### Features:
- Calculation of shock angles and downstream properties.
- Determination of Mach number and pressure ratios.

#### Usage:
```python
python ObliqueShockCalculator.py
```

### Prandtl-MeyerCalculator.py

This script calculates properties of Prandtl-Meyer expansion waves, which describe the expansion process in supersonic flows.

#### Features:
- Calculation of flow properties through expansion waves.
- Determination of Mach number before and after the expansion fan.

#### Usage:
```python
python Prandtl-MeyerCalculator.py
```

## Requirements

- Python 3.x
- NumPy
- SciPy

You can install the required libraries using pip:
```bash
pip install numpy scipy
```

## Author

- Oreoluwa Abejide


