import ialex
import Array
import KalmanFMWK

class NEXTKalmanFilter(KalmanFMWK.KalmanFilter):
    self.N = 4
    
    def TransportMatrix( self, index ):
        return Array.Identity( self.N )
    
    def MeasurementMatrix( self, index ):
        zmean = 0.5 * ( self.Track.GetNode(index).hit[2] + self.Track.GetNode(index+1).hit[2] )
        return Array.Matrix( [1.,0.,z0,0.],
                             [0.,1.,0.,z0],
                             [ ?, ?, ?, ?])
    
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
        return Array.Identity( self.N )



class PyKal( ialex.IAlg ):
    def __init__( self, measurements = [], name = 'NEXTKalmanFilter' ):
        ialex.IAlg.__init__( self, name )
        self.measurements = copy.deepcopy(measurements)
    
    def define( self ):
        self.p0 = 2.49 #MeV
    
    def initialize( self ):
        return
    
    def execute( self ):
        
        return True
    
    def finalize( self ):
        return



if __name__ = '__main__':
    import ROOT
    
    ### Free fall example
    class FreeFall( KalmanFMWK.KalmanFilter ):
        #### State variables:
        #### [ y_k, y'_k, g ] where ' denotes time derivative
        
        def TransportMatrix( self, index ):
            return Array.Matrix( [ 1., 1., -0.5],
                                 [ 0., 1., -1.0],
                                 [ 0., 0., +1.0])
        
        def MeasurementMatrix( self, index ):
            return Array.Matrix( [1.,0.,0.] )
        
        def MultipleScatteringMatrix( self, index ):
            return Array.Zeros( 3, 3 )
        
        def NoiseMatrix( self, index ):
            return Array.Identity( self.N )

    def GenerateData( h = 100., sigma = 1., N = 100 ):
        N = float(N)
        R = ROOT.TRandom3(0)
        return [ ( Array.Vector(h - i/N + R.Gaus(0,sigma)), Array.Matrix( [sigma**2] ) ) for i in range(N) ]

    FreeFallExample = FreeFall()
    
    measurements = GenerateData()
    initialstate = KalmanMeasurement( Array.Vector( 100., 0., 9.81 ), Array.Matrix( [1.,0.,0.],[0.,1.,0.],[0.,0.,0.] ) )

    FreeFallExample.SetMeasurements( measurements )
    FreeFallExample.SetInitialState( initialstate )
    FreeFallExample.SetInitialGuess()
