import ialex
import Array
import math
import KalmanFMWK
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

class NEXTKalmanFilter(KalmanFMWK.KalmanFilter):
    Ndim = 4
    
    def TransportMatrix( self, index ):
        return Array.Matrix([1.,0.,z0,0.],
                            [0.,1.,0.,z0],
                            [0.,0.,1.,0.],
                            [0.,0.,0.,1.])


    def MeasurementMatrix( self, index ):
        return Array.Matrix( [1.,0.,z0,0.],
                             [0.,1.,0.,z0])
    
    def MultipleScatteringMatrix( self, index ):
        dz =self.system.Nodes[k].ZSlice.dz
        z0 =self.system.Nodes[k].ZSlice.Zi
        z02=z0*z0;
        edep = self.system.Nodes[k].ZSlice.Edep*1.
        Lr = self.system.Nodes[k].ZSlice.Lr
        L=abs(dz)/Lr
        
        stateDict=self.system.States[k-1]
        state =stateDict["F"]
        
        ak = state.V
        p1 = ak[0]
        p2 = ak[1]
        p3 = ak[2]
        p4 = ak[3]
        
        p3p3 = s2tms*(1+p3**2)*(1+p3**2+p4**2)
        p4p4 =  s2tms*(1+p4**2)*(1+p3**2+p4**2)
        p3p4 = s2tms*p3*p4*(1+p3**2+p4**2)
        
        return Array.Matrix( [  z02*p3p3,  z02*p3p4, -z0*p3p3, -z0*p3p4],
                             [  z02*p3p4,  z02*p4p4, -z0*p3p4, -z0*p4p4],
                             [ -z0 *p3p3, -z0 *p3p4,     p3p3,     p3p4],
                             [ -z0 *p3p4, -z0 *p4p4,     p3p4,     p4p4])

    def NoiseMatrix( self, index ):
        return Array.Identity( 2 ) * 0.01

class NEXTKalmanFilterBis(KalmanFMWK.KalmanFilter):
    '''
        State = ( x - dz tanx, y - dz tany, tanx, tany )
    '''
    Ndim = 4
    xyresolution = 0.1 # cm
    
    def __init__( self, p0 = 2.49 ): # MeV
        KalmanFMWK.KalmanFilter.__init__( self, name = 'NEXT Kalman Filter' )
        self.p0 = p0
        self.p  = p0
    
    def TransportMatrix( self, index ):
        return Array.Identity( self.Ndim )
    
    def MeasurementMatrix( self, index ):
        dz = self.Track.GetNode(index).running - self.Track.GetNode(index-1).running
        return Array.Matrix( [ 1., 0., dz, 0. ],
                             [ 0., 1., 0., dz ])
    
    def MultipleScatteringMatrix( self, index ):
        z0   = self.Track.GetNode(index).running
        z02  = z0**2;
        edep = self.this_hit[-1]
        L    = abs( self.next_node.running - z0 ) / Physics.Xe.x0
        
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
        return Array.Matrix( [ [ 1., 0. ],
                               [ 0., 1. ]] ) * self.xyresolution**2

class PyKal( ialex.IAlg ):
    def __init__( self, measurements = [], name = 'NEXTKalmanFilter' ):
        ialex.IAlg.__init__( self, name )
        self.measurements = copy.deepcopy(measurements)
    
    def define( self ):
        self.p0 = 2.49 #MeV
    
    def initialize( self ):
        self.KF = NEXTKalmanFilterBis( self.p0 )
        return
    
    def execute( self ):
        
        return True
    
    def finalize( self ):
        return



if __name__ == '__main__':
    import ROOT
    R            = ROOT.TRandom3(0)
    V            = Array.Matrix( [0.1**2,0.], [0.,0.1**2] )
    hits         = [ Array.Vector( R.Gaus(0,1), R.Gaus(0,1) ) for i in range(100) ]
    runnings     = map( float, range(100) )
    measurements = [ KalmanFMWK.KalmanMeasurement( h, V ) for h in hits ]
    istate       = Array.Vector( 0., 0., 0., 0. )
    icvmatrix    = Array.Identity(4)
    
    nextkf = NEXTKalmanFilterBis( 2.49 )
    nextkf.SetMeasurements( runnings, measurements )
    nextkf.SetInitialState( KalmanFMWK.KalmanMeasurement( istate, icvmatrix ) )
    track = nextkf.Fit()
    p = track.Plot()
    raw_input()

#    import ROOT
#    
#    ### Free fall example
#    class FreeFall( KalmanFMWK.KalmanFilter ):
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
#    initialstate = KalmanMeasurement( Array.Vector( 100., 0., 9.81 ), Array.Matrix( [1.,0.,0.],[0.,1.,0.],[0.,0.,0.] ) )
#
#    FreeFallExample.SetMeasurements( measurements )
#    FreeFallExample.SetInitialState( initialstate )
#    FreeFallExample.SetInitialGuess()
