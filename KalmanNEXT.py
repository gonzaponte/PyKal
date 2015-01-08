import Array
import math
from KalmanFilter import KalmanFilter, KalmanData
import Physics

def Beta( p, m = Physics.Constants.me ):
    '''
        Return the beta value for a given momentum.
    '''
    return p / math.sqrt( p**2 + m**2 )

def thetaMS( p, x, Q = 1 ):
    '''
        Return the angle of multiple scattering for a certain momentum p in MeV/c and a number of radiation lengths x.
    '''
    return 13.6 / ( Beta(p) * p ) * Q * math.sqrt( x ) * ( 1 + 0.038 * math.log( x ) )

class NEXTKF(KalmanFilter):
    '''
        State = ( x - dz tanx, y - dz tany, tanx, tany )
    '''
    Ndim = 4
    xyresolution = 0.1 # cm
    
    def __init__( self, p0 = 2.49 ): # MeV
        KalmanFilter.__init__( self, name = 'NEXT Kalman Filter' )
        self.p0 = p0
        self.p  = p0
    
    def TransportMatrix( self, index ):
        return Array.Identity( self.Ndim )
    
    def MeasurementMatrix( self, index ):
        dz = self.Track.GetNode(index).running - self.Track.GetNode(index-1).running
        dz = - dz
        return Array.Matrix( [ 1., 0., dz, 0. ],
                             [ 0., 1., 0., dz ])
    
    def MultipleScatteringMatrix( self, index ):
        z0   = self.Track.GetNode(index).running
        z02  = z0**2;
        edep = self.this_hit[-1]
        try:
            L    = abs( self.Track.GetNode(index+1).running - z0 ) / Physics.Xe.x0
        except:
            L    = 1e3
        
        p1, p2, p3, p4 = self.prev_state
        
        theta02 = thetaMS( self.p , L )**2
                            
        p3p3 =  theta02 * ( 1 + p3**2 ) * ( 1 + p3**2 + p4**2 )
        p4p4 =  theta02 * ( 1 + p4**2 ) * ( 1 + p3**2 + p4**2 )
        p3p4 =  theta02 *  p3 * p4      * ( 1 + p3**2 + p4**2 )
        
        return Array.Matrix( [  z02 * p3p3,  z02 * p3p4, -z0 * p3p3, -z0 * p3p4 ],
                             [  z02 * p3p4,  z02 * p4p4, -z0 * p3p4, -z0 * p4p4 ],
                             [ -z0  * p3p3, -z0  * p3p4,       p3p3,       p3p4 ],
                             [ -z0  * p3p4, -z0  * p4p4,       p3p4,       p4p4 ])
    
    def NoiseMatrix( self, index ):
        return Array.Identity(2) * self.xyresolution**2

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
    hits         = [ Array.Vector( math.sin(0.01*i) + R.Gaus(0,.1), -math.cos(0.1*i) + R.Gaus(0,.1) ) for i in range(Nhits) ]
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

#    import ROOT
#    
#    ### Free fall example
#    class FreeFall( KalmanFilter ):
#        #### State variables:
#        #### [ y_k, y'_k, g ] where ' denotes time derivative
#        
#        def TransportMatrix( self, index ):
#            return Array.Matrix( [ 1., 1., -0.5],
#                                 [ 0., 1., -1.0],
#                                 [ 0., 0., +1.0])
#        
#        def MeasurementMatrix( self, index ):
#            return Array.Matrix( [1.,0.,0.] )
#        
#        def MultipleScatteringMatrix( self, index ):
#            return Array.Zeros( 3, 3 )
#        
#        def NoiseMatrix( self, index ):
#            return Array.Identity( self.N )
#
#    def GenerateData( h = 100., sigma = 1., N = 100 ):
#        N = float(N)
#        R = ROOT.TRandom3(0)
#        return [ ( Array.Vector(h - i/N + R.Gaus(0,sigma)), Array.Matrix( [sigma**2] ) ) for i in range(N) ]
#
#    FreeFallExample = FreeFall()
#    
#    measurements = GenerateData()
#    initialstate = KalmanData( Array.Vector( 100., 0., 9.81 ), Array.Matrix( [1.,0.,0.],[0.,1.,0.],[0.,0.,0.] ) )
#
#    FreeFallExample.SetMeasurements( measurements )
#    FreeFallExample.SetInitialState( initialstate )
#    FreeFallExample.SetInitialGuess()
