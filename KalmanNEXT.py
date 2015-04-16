import Array
from math import *
from KalmanFilter import KalmanFilter, KalmanData
from Physics import Constants, NEXT

NEXT = NEXT( P = 10., T = 20. )

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

class NEXTKF(KalmanFilter):
    '''
        State = ( x - dz tanx, y - dz tany, tanx, tany, p )
    '''
    Sdim = 5 # State dimensionality
    xyresolution = 0.1 # cm
    energyresolution = 0.5 # MeV
    xy2 = xyresolution**2
    Er2 = energyresolution**2
    
    def __init__( self, p0 = 2.49 ): # MeV
        KalmanFilter.__init__( self, name = 'NEXT Kalman Filter' )
        self.p  = p0
    
    def TransportMatrix( self, index ):
        F   = Array.Identity( Sdim )
        F[-1][-1] = self.Eloss( [-1], dx )

        s   = self.Track.GetNode(index-1).Fstate.Vector
        dx  = self.Track.GetNode(index).running - self.Track.GetNode(index-1).running
        dx *= ( 1 + s[2]**2  +s[3]**2 ) ** 0.5


        return
    
    def MeasurementMatrix( self, index ):
        dz = self.Track.GetNode(index).running - self.Track.GetNode(index-1).running
        dz = - dz
        return Array.Matrix( [ 1., 0., dz, 0., 0. ],
                             [ 0., 1., 0., dz, 0. ],
                             [ 0., 0., 0., 0., 1. ])
    
    def MultipleScatteringMatrix( self, index ):
        z0   = self.Track.GetNode(index).running
        z02  = z0**2;
        edep = self.this_hit[-1]
        try:
            L    = abs( self.Track.GetNode(index+1).running - z0 ) / NEXT.x0
        except:
            L    = 1e3
        
        x0, y0, tx, ty, p = self.Pstate.State
        theta02 = thetaMS( p , L )**2
        
        factor = theta02 * ( 1 + tx**2 + ty**2 )
        p33 =  factor * ( 1 + tx**2 )
        p44 =  factor * ( 1 + ty**2 )
        p34 =  factor *  tx * ty
        
        return Array.Matrix( [  z02 * p33,  z02 * p34, -z0 * p33, -z0 * p34 ],
                             [  z02 * p34,  z02 * p44, -z0 * p34, -z0 * p44 ],
                             [ -z0  * p33, -z0  * p34,       p33,       p34 ],
                             [ -z0  * p34, -z0  * p44,       p34,       p44 ])
    
    def NoiseMatrix( self, index ):
        return Array.Identity(2) * xyresolution**2

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
        return Array.Identity(4)
    
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
