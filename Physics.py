class Constants:
    c    = 299792458.    # speed of light in m/s
    e    = 1.6021766e-19 # electron charge in C
    me   = 0.510998928   # electron mass in MeV/c2
    rho0 = 2.6867774e19  # density of a gas at T = 0 K and P = 1 atm in 1/cm3
    NA   = 6.02214179e23 # Avogadro constant in 1/mol

class Xe:
    P    = 10. / 1.01325    # Xe pressure in bar
    T    = 293.15 / 273.15  # Xe temperature relative to 0 Celsius
    X0   = 8.48             # Xe radiation length in g / cm2
    mass = 131.293          # Xe nuclear mass in amu
    rho  = Constants.rho0 * P / T * mass / Constants.NA # Xe density in 1/cm3
    x0   = X0 / rho         # Xe radiation length in cm
