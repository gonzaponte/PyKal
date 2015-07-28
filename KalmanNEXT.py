import Array
from math import *
from KalmanFilter import KalmanFilter, KalmanData
from Physics import Constants, NEXT, K2P
from numpy import loadtxt
from scipy.interpolate import interp1d
from scipy.integrate import quad as integrate

def Beta( p, m = Constants.me ):
    '''
        Return the beta value for a given momentum.
    '''
    return p / sqrt( p**2 + m**2 )

def thetaMS( p, x, Q = 1 ):
    '''
        Return the angle of multiple scattering for a certain momentum p in MeV/c and a number of radiation lengths x.
    '''
    return 13.6 / ( Beta(p) * p ) * Q * sqrt( x ) * ( 1.0 + 0.038 * log( x ) )

_Eloss_data = loadtxt('/Users/Gonzalo/github/PyKalman/KalmanFilter/xe_estopping_power_NIST.dat')

_ds, _dEds = zip(*_Eloss_data)

_Eloss = interp1d( (0,) + _ds, (_dEds[0],) + _dEds, kind = 'cubic' )

class NEXTKF(KalmanFilter):
    '''
        State = ( x - dz tanx, y - dz tany, tanx, tany, E )
    '''
    xyresolution = 0.1 # cm
    Eresolution  = 0.5 # MeV
    
    xyres2 = xyresolution**2
    Eres2  = Eresolution**2
    
    def __init__( self, P = 15., T = 27. ):
        KalmanFilter.__init__( self, name = 'NEXT Kalman Filter' )
        self.gas = NEXT( P, T )
    
    def TransportMatrix( self, index0, index1 ):
        This  = self.Track.GetNode(index0)
        Other = self.Track.GetNode(index1)
        
        state  = This.Pstate.State if This.status < 2 else This.Fstate.State
        tx = state[2]
        ty = state[3]
        dz = Other.running - This.running
        ds = abs(dz) * ( tx**2 + ty**2 + 1 ) ** 0.5
        
        F = Array.Identity( self.Sdim )
        F[-1][-1] = self.Eloss( state[-1], ds )
        return F

    
    def MeasurementMatrix( self, index0, index1 ):
        dz = self.Track.GetNode(index1).running - self.Track.GetNode(index0).running
        return Array.Matrix( [ 1., 0., dz, 0., 0. ],
                             [ 0., 1., 0., dz, 0. ],
                             [ 0., 0., 0., 0., 1. ])

    def MultipleScatteringMatrix( self, index0, index1 ):
        This  = self.Track.GetNode(index0)
        Other = self.Track.GetNode(index1)
        z0    = This.running
        z02   = z0**2;
        z1    = Other.running
        L     = abs( z1 - z0 ) / self.gas.x0 * 0.1

        dx, dy, tx, ty, e = This.Pstate.State
        theta02 = thetaMS( K2P(e) , L )**2
        
        factor = theta02 * ( 1 + tx**2 + ty**2 )
        p33 =  factor * ( 1 + tx**2 )
        p44 =  factor * ( 1 + ty**2 )
        p34 =  factor *  tx * ty
        
        return Array.Matrix( [  z02 * p33,  z02 * p34, -z0 * p33, -z0 * p34,      0.    ],
                             [  z02 * p34,  z02 * p44, -z0 * p34, -z0 * p44,      0.    ],
                             [ -z0  * p33, -z0  * p34,       p33,       p34,      0.    ],
                             [ -z0  * p34, -z0  * p44,       p34,       p44,      0.    ],
                             [         0.,         0.,        0.,        0., self.Eres2 ])
    
    def NoiseMatrix( self, index ):
        V = Array.Identity(self.Mdim) * self.xyres2
        V[-1][-1] = self.Eres2
        return V

    def Eloss( self, e0, ds ):
        return _Eloss( e0 ) * ds

class KFNEXT(KalmanFilter):
    '''
        State = ( x, y, tanx, tany, E )
    '''
    xyresolution = 0.1 # cm
#    Eresolution  = 0.01 #/ 2.458 # MeV
    Eresolution  = 0.01 * 2.458**0.5 #/ 2.458 # MeV
    
    xyres2 = xyresolution**2
    Eres2  = Eresolution**2
    
    def __init__( self, P = 15., T = 27. ):
        KalmanFilter.__init__( self, name = 'NEXT Kalman Filter' )
        self.gas = NEXT( P, T )
    
    def SetMeasurements( self, runs, hits ):
        KalmanFilter.SetMeasurements( self, runs, hits )
        
        totalE = sum( hit.State[-1] for hit in hits )
        for i,hit in enumerate(hits):
            self.Track.GetNode(i).ParticleEnergy = totalE
            totalE -= hit.State[-1]
    
    def TransportMatrix( self, index0, index1 ):
        This  = self.Track.GetNode(index0)
        Other = self.Track.GetNode(index1)
        
        state  = This.Pstate.State if This.status < 2 else This.Fstate.State
        tx = state[2]
        ty = state[3]
        dz = Other.running - This.running
        ds = abs(dz) * ( tx**2 + ty**2 + 1 ) ** 0.5
        
        F = Array.Identity( self.Sdim )
        F[0][2]   = dz
        F[1][3]   = dz
        F[-1][-1] = self.Eloss( This.ParticleEnergy, ds ) / This.hit.State[-1]
        return F
    
    def MeasurementMatrix( self, index0, index1 ):
        return Array.Matrix( [ 1., 0., 0., 0., 0. ],
                             [ 0., 1., 0., 0., 0. ],
                             [ 0., 0., 0., 0., 1. ])
    
    def MultipleScatteringMatrix( self, index0, index1 ):
        This  = self.Track.GetNode(index0)
        Other = self.Track.GetNode(index1)

        dx, dy, tx, ty, e = This.Pstate.State

        z0    = This.running
        z02   = z0**2;
        z1    = Other.running
        dz    = z1 - z0
        ds    = abs(dz) * ( tx**2 + ty**2 + 1 ) ** 0.5
        L     = ds / self.gas.x0
        theta02 = thetaMS( K2P(This.ParticleEnergy) , L )**2
        
        factor = theta02 * ( 1 + tx**2 + ty**2 )
        p33 =  factor * ( 1 + tx**2 )
        p44 =  factor * ( 1 + ty**2 )
        p34 =  factor *  tx * ty
        
        state  = This.Pstate.State if This.status < 2 else This.Fstate.State
        self.Eres2 = self.Eresolution**2 / state[-1]

        return Array.Matrix( [  z02 * p33/3,  z02 * p34/3, -z0 * p33/2, -z0 * p34/2,      0.    ],
                             [  z02 * p34/3,  z02 * p44/3, -z0 * p34/2, -z0 * p44/2,      0.    ],
                             [ -z0  * p33/2, -z0  * p34/2,       p33  ,       p34  ,      0.    ],
                             [ -z0  * p34/2, -z0  * p44/2,       p34  ,       p44  ,      0.    ],
                             [         0.,         0.,        0.,           0.     , self.Eres2 ])# * 10.
    
    def NoiseMatrix( self, index ):
        V = Array.Identity(self.Mdim) * self.xyres2
        V[-1][-1] = self.Eresolution**2 / self.Track.GetNode(index).hit.State[-1]
        return V
    
    def Eloss( self, e0, ds ):
        return _Eloss( e0 ) * ds

#class KFNEXT(KalmanFilter):
#    '''
#        State = ( x, y, tanx, tany, E )
#    '''
#
#    xyresolution = 0.1 # cm
#    Eresolution  = 0.5 # MeV
#
#    xyres2 = xyresolution**2
#    Eres2 = Eresolution**2
#    
#    def __init__( self, P = 15., T = 27. ): # MeV
#      KalmanFilter.__init__( self, name = 'NEXT Kalman Filter' )
#      self.gas = NEXT( P, T )
#    
#    def TransportMatrix( self, index0, index1 ):
#      F   = Array.Identity( self.Sdim )
#      F[0][2] = F[1][3] = dz
#      Prev  = self.Track.GetNode(index0)
#      This = self.Track.GetNode(index1)
#      s   = Prev.Fstate.Vector
#      dz  = (This.running - Prev.running)**2 * ( 1 + s[2]**2  +s[3]**2 ) * thetaMS( p, L )**2
#      
#      p3p3 = dz#...
#      
#      return
#    
#    def MeasurementMatrix( self, index0, index1 ):
#      dz = self.Track.GetNode(index1).running - self.Track.GetNode(index0).running
#      dz = - dz
#      return Array.Matrix( [ 1., 0., dz, 0., 0. ],
#                          [ 0., 1., 0., dz, 0. ],
#                          [ 0., 0., 0., 0., 1. ])
#    
#    def MultipleScatteringMatrix( self, index ):
#      z0   = self.Track.GetNode(index).running
#      z02  = z0**2;
#      edep = self.this_hit[-1]
#      try:
#        L    = abs( self.Track.GetNode(index+1).running - z0 ) / self.gas.x0
#      except:
#        L    = 1e3
#      
#      x0, y0, tx, ty, p = self.Pstate.State
#      theta02 = thetaMS( p , L )**2
#      
#      factor = theta02 * ( 1 + tx**2 + ty**2 )
#      p33 =  factor * ( 1 + tx**2 )
#      p44 =  factor * ( 1 + ty**2 )
#      p34 =  factor *  tx * ty
#      
#      return Array.Matrix( [  z02 * p33,  z02 * p34, -z0 * p33, -z0 * p34 ],
#                          [  z02 * p34,  z02 * p44, -z0 * p34, -z0 * p44 ],
#                          [ -z0  * p33, -z0  * p34,       p33,       p34 ],
#                          [ -z0  * p34, -z0  * p44,       p34,       p44 ])
#    
#    def NoiseMatrix( self, index ):
#      return Array.Identity(2) * xyresolution**2
#

class TrigoKF( KalmanFilter ):
    '''
        State = ( x, y, tantheta, tanphi )
    '''
    Ndim = 4
    xyresolution = 0.1 # cm
    
    def __init__( self ):
        KalmanFilter.__init__( self, name = 'Trigo Kalman Filter' )
    
    def TransportMatrix( self, index ):
        z = self.Track.GetNode(index).running
        return Array.Matrix( [ 1., 0., z , 0. ],
                             [ 0., 1., z , 0. ],
                             [ 0., 0., 1., 0. ],
                             [ 0., 0., 0., 1. ])
    
    def MeasurementMatrix( self, index ):
        return Array.Matrix( [ 1., 0., 0., 0. ],
                             [ 0., 1., 0., 0. ])
    
    def MultipleScatteringMatrix( self, index ):
        return Array.Identity(4) * 1e-12
    
    def NoiseMatrix( self, index ):
        return Array.Identity(2) * self.xyresolution**2



if __name__ == '__main__':
    import ROOT
    R            = ROOT.TRandom3(0)
    V            = Array.Identity(2) * 0.01
    Nhits        = 500
    hits         = [ Array.Vector( sin(0.01*i) + R.Gaus(0,.1), -cos(0.1*i) + R.Gaus(0,.1) ) for i in range(Nhits) ]
    runnings     = map( float, range(Nhits) )
    measurements = [ KalmanData( h, V ) for h in hits ]
    istate       = Array.Vector( hits[0][0], hits[0][1], 0., 0. )
    istate       = Array.Vector( 0., 0., 0., 0. )
    icvmatrix    = Array.Identity(4) * 2
    
    nextkf = NEXTKF()
    nextkf = TrigoKF()
    nextkf.SetMeasurements( runnings, measurements )
    nextkf.SetInitialState( KalmanData( istate, icvmatrix ) )
    track = nextkf.Fit()
    p = track.Plot()
    p3 = track.Plot3D()
    cchi = ROOT.TCanvas()
    chi2 = ROOT.TGraph()
    chi2.SetMarkerStyle(20)
    for i in range(Nhits):
        chi2.SetPoint(i,i,track.GetNode(i).chi2)
    chi2.Draw('AP')
    raw_input()
