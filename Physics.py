class Constants:
    c    = 299792458.    # speed of light in m/s
    e    = 1.6021766e-19 # electron charge in C
    me   = 0.510998928   # electron mass in MeV/c2
    rho0 = 2.6867774e19  # density of a gas at T = 0 K and P = 1 atm in 1/cm3
    NA   = 6.02214179e23 # Avogadro constant in 1/mol

class NEXT:
    '''
        NEXT experiment properties.
    '''
    def __init__( self, P = 15., T = 20. ):
        '''
            P -> Xenon pressure in atm.
            T -> Xenon temperature in Celsius.
            '''
        self.P   = P / 1.01325        # Xe pressure in bar
        self.T   = 1. + T/273.15      # Xe temperature relative to 0 celsius
        self.X0  = 8.48               # Xe radiation length in g / cm2
        self.A   = 131.293            # Xe nuclear mass in amu
        self.rho = Constants.rho0 * self.P / self.T * self.A / Constants.NA # Xe density in 1/cm3
        self.x0  = self.X0 / self.rho # Xe radiation length in cm

def K2P( k, m = Constants.me ):
    return ( k * ( k + 2*m ) ) ** 0.5

