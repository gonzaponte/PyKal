'''
    Kalman filter framework.
'''
from copy import copy as _copy, deepcopy as _deepcopy
from Array import Vector as _Vector, Matrix as _Matrix, Identity as _IdentityMatrix

class KalmanData:
    '''
        Represents a state of the Kalman Filter with the associated covariance. It is a pair (vector, covariance matrix).
    '''

    def __init__( self, vector = _Vector(), CovarianceMatrix = _Matrix() ):
        ''' 
            Initialize with a state vector and its covariance matrix.
        '''
        self.State = vector
        self.CovarianceMatrix = CovarianceMatrix
    
    def __len__( self ):
        '''
            Return the dimension of the state vector.
        '''
        return len(self.State)
    
    def __bool__( self ):
        '''
            False if empty.
        '''
        return bool( len(self) )

    def __str__( self ):
        '''
            String representation for printing purposes.
        '''
        return '''
Vector = {0}
Cov =
{1}
'''.format( self.State, self.CovarianceMatrix )

class KalmanNode:
    '''
        Contains all the information about a given step. Contains the index of the step, the hit, the predicted, filtered and smoothed states and residuals as well as the chi2 and the cumulative chi2 of the step.
    '''
    
    def __init__( self,
                  step       = 0,      running = 0, hit = KalmanData(), true_hit   = KalmanData(),
                  Pstate     = KalmanData(), Fstate     = KalmanData(), Sstate     = KalmanData(),
                  Presiduals = KalmanData(), Fresiduals = KalmanData(), Sresiduals = KalmanData(),
                  chi2       = -1, cumchi2 = -1 ):

        self.step       = step
        self.running    = running
        self.hit        = hit
        self.true_hit   = true_hit
        self.Pstate     = Pstate
        self.Fstate     = Fstate
        self.Sstate     = Sstate
        self.Presiduals = Presiduals
        self.Fresiduals = Fresiduals
        self.Sresiduals = Sresiduals
        self.chi2       = chi2
        self.cumchi2    = cumchi2

    def __str__( self ):
        '''
            String representation for printing purposes.
        '''
        return '''
step number {0}
running variable value {1}
true hit:
{2}
hit:
{3}
predicted state:
{4}
filtered state:
{5}
smoothed state:
{6}
predicted residuals:
{7}
filtered residuals:
{8}
smoothed residuals:
{9}
chi2 of this node:
{10}
cumulative chi2 at this node:
{11}
'''.format( self.step, self.running, self.true_hit, self.hit,
            self.Pstate, self.Fstate, self.Sstate,
            self.Presiduals, self.Fresiduals, self.Sresiduals,
            self.chi2, self.cumchi2 )

class KalmanTrack:
    '''
        A set of KalmanNodes. Contains all possible information about the track.
    '''
    def __init__( self, nodes = list() ):
        '''
            Initialize with the list of nodes. It is also possible to add them later with the AddNode method.
        '''
        self.nodes  = _deepcopy(nodes)
        self.Nnodes = len( self.nodes )

    def AddNode( self, node ):
        '''
            Add a node.
        '''
        self.nodes.append( _copy( node ) )
        self.Nnodes += 1
    
    def GetNode( self, index ):
        '''
            Return the node.
        '''
        return self.nodes[index]
    
    def Plot( self ):
        '''
            Plot the track.
        '''
        import ROOT
        xmeasurements = ROOT.TGraph()
        ymeasurements = ROOT.TGraph()
        xfiltered     = ROOT.TGraph()
        yfiltered     = ROOT.TGraph()
        xsmoothed     = ROOT.TGraph()
        ysmoothed     = ROOT.TGraph()
        
        xmeasurements.SetMarkerStyle(20)
        ymeasurements.SetMarkerStyle(20)
        xfiltered.SetMarkerStyle(29)
        yfiltered.SetMarkerStyle(29)
        xsmoothed.SetMarkerStyle(34)
        ysmoothed.SetMarkerStyle(34)

        xmeasurements.SetMarkerColor(1)
        ymeasurements.SetMarkerColor(1)
        xfiltered.SetMarkerColor(3)
        yfiltered.SetMarkerColor(3)
        xsmoothed.SetMarkerColor(34)
        ysmoothed.SetMarkerColor(34)

        for node in self.nodes[:-1]:
            xmeasurements.SetPoint( node.step, node.running, node.hit.State[0] )
            ymeasurements.SetPoint( node.step, node.running, node.hit.State[1] )
            xfiltered    .SetPoint( node.step, node.running, node.Fstate.State[0] )
            yfiltered    .SetPoint( node.step, node.running, node.Fstate.State[1] )
            xsmoothed    .SetPoint( node.step, node.running, node.Sstate.State[0] )
            ysmoothed    .SetPoint( node.step, node.running, node.Sstate.State[1] )

        cx = ROOT.TCanvas()
        xmeasurements.Draw('AP')
        xfiltered.Draw('P')
        xsmoothed.Draw('P')

        cy = ROOT.TCanvas()
        ymeasurements.Draw('AP')
        yfiltered.Draw('P')
        ysmoothed.Draw('P')
        
        return cx, xmeasurements, xfiltered, xsmoothed, cy, ymeasurements, yfiltered, ysmoothed
    
    def Plot3D( self ):
        import ROOT
        measurements = ROOT.TGraph2D()
        filtered     = ROOT.TGraph2D()
        smoothed     = ROOT.TGraph2D()
        
        measurements.SetMarkerStyle(20)
        filtered    .SetLineWidth(2)
        smoothed    .SetLineWidth(2)
        
        measurements.SetMarkerColor(1)
        filtered    .SetLineColor(2)
        smoothed    .SetLineColor(4)
        
        for node in self.nodes[:-1]:
            measurements.SetPoint( node.step, node.running, node.hit.State[0]   , node.hit.State[1]     )
            filtered    .SetPoint( node.step, node.running, node.Fstate.State[0], node.Fstate.State[1] )
            smoothed    .SetPoint( node.step, node.running, node.Sstate.State[0], node.Sstate.State[1] )
        
        c = ROOT.TCanvas()
        measurements.Draw('AP')
        filtered    .Draw('Clinesame')
        smoothed    .Draw('Clinesame')
        
        return c, measurements, filtered, smoothed
    
    def __str__( self ):
        '''
            String representation for printing purposes.
        '''
        return '\nNumber of steps: {0} \n\n'.format( self.Nnodes ) + '\n\n'.join( map( str, self.nodes ) )

class KalmanFilter:
    '''
        Abstract implementation of the Kalman Filter.
    '''
    
    def __init__( self, name = 'KalmanFilter' ):
        '''
            Initializer. It is needed to call the SetInitialState and SetMeasurements methods mandatorily.
        '''
        self.name         = name
        self.Track        = KalmanTrack()
        self.IState       = False # To ensure all the information is present when fitting.
        self.Measurements = False # To ensure all the information is present when fitting.
        self.TrueHits     = False

    def SetInitialState( self, state = KalmanData() ):
        '''
            Sets "state" as the initial state of the Kalman filter, both as the predicted state and the filtered state.
        '''
        First            = self.Track.GetNode(0)
        First.Pstate     = state
        First.Fstate     = state
        First.Presiduals = First.hit.State * 0.
        First.Fresiduals = First.hit.CovarianceMatrix * 0.
        First.chi2       = 0.
        First.cumchi2    = 0.
        self.Sdim        = len(state)
        self.Mdim        = len(First.hit.State)
        self.IState      = True
    
    def SetTrue( self, hits ):
        '''
            Sets "hits" as the true hits of the track.
        '''
        assert self.HaveMeasurements, 'True hits cannot be setted before measurements. True hits not included.'
        
        indices = range( len(hits) )
        for i,hit in enumerate(hits):
            self.Track.GetNode(i).true_hit = hit
        self.TrueHits = True

    def SetMeasurements( self, runs, hits ):
        '''
            Sets "hits" as the measurements of the track for running variables runs.
        '''
        for i,(r,hit) in enumerate( zip( runs, hits ) ):
            self.Track.AddNode( KalmanNode( step = i, running = r, hit = hit) )
        
        self.Measurements = True # Now we can compute stuff

    def TransportMatrix( self, index ):
        '''
            To be coded in the particular implementation.
        '''
        raise NotImplementedError('TransportMatrix is not defined.')
        return
    
    def MeasurementMatrix( self, index ):
        '''
            To be coded in the particular implementation.
        '''
        raise NotImplementedError('MeasurementMatrix is not defined.')
        return
    
    def MultipleScatteringMatrix( self, index ):
        '''
            To be coded in the particular implementation.
        '''
        raise NotImplementedError('MultipleScatteringMatrix is not defined.')
        return

    def NoiseMatrix( self, index ):
        '''
            To be coded in the particular implementation.
        '''
        raise NotImplementedError('NoiseMatrix is not defined.')
        return
    
    def Filter( self, index ):
        '''
            Filter the index-th step.
        '''
    
        if not index:
            raise ValueError( 'The initial step cannot be predicted' )
        
        ### Predicting
        Last = self.Track.GetNode( index - 1 )
        This = self.Track.GetNode( index )
        Hit  = This.hit.State

        TMatrix   = self.TransportMatrix( index )
        MSMatrix  = self.MultipleScatteringMatrix( index )
        MMMatrix  = self.MeasurementMatrix( index )
        NMatrix   = self.NoiseMatrix( index )
        
        TMatrixT  = TMatrix.T()
        MMMatrixT = MMMatrix.T()
        NMatrixI  = NMatrix.Inverse()
        
        x0 = Last.Fstate.State
        C0 = Last.Fstate.CovarianceMatrix
        
        x_predicted = TMatrix ** x0
        C_predicted = TMatrix ** C0 ** TMatrixT + MSMatrix
        
        r_predicted = Hit - MMMatrix ** x_predicted
        R_predicted = NMatrix + MMMatrix ** C_predicted ** MMMatrixT
        
        This.Pstate     = KalmanData( x_predicted, C_predicted )
        This.Presiduals = KalmanData( r_predicted, R_predicted )
        
        ### Filtering
        C_filtered = ( C_predicted.Inverse() + MMMatrixT ** NMatrixI ** MMMatrix ).Inverse()
        GainMatrix = C_filtered ** MMMatrixT ** NMatrixI
        x_filtered = x0 + GainMatrix ** ( Hit - MMMatrix ** x0 )
        
        projector  = _IdentityMatrix( self.Mdim ) - MMMatrix ** GainMatrix
        r_filtered = projector ** r_predicted
        R_filtered = projector ** NMatrix
        
        chi2plus = r_filtered ** R_filtered.Inverse() ** r_filtered
        newchi2  = Last.cumchi2 + chi2plus
        
        This.Fstate     = KalmanData( x_filtered, C_filtered )
        This.Fresiduals = KalmanData( r_filtered, R_filtered )
        This.chi2       = chi2plus
        This.cumchi2    = newchi2
    
    def Smooth( self, index ):
        '''
            Smooth the index-th step.
        '''
        Next = self.Track.GetNode( index + 1 )
        This = self.Track.GetNode( index     )
        Hit  = This.hit.State
        
        MMMatrix  = self.MeasurementMatrix( index )
        NMatrix   = self.NoiseMatrix( index )
        MMMatrixT = MMMatrix.T()
        TMatrixT  = self.TransportMatrix( index ).T()
        
        x0  = This.Fstate.State
        C0  = This.Fstate.CovarianceMatrix
        x1p = Next.Pstate.State
        C1p = Next.Pstate.CovarianceMatrix
        x1s = Next.Sstate.State
        C1s = Next.Sstate.CovarianceMatrix
        
        GainMatrix = C0 ** TMatrixT   ** C1p.Inverse()
        x_smooth   = x0  + GainMatrix ** ( x1s - x1p )
        C_smooth   = C0  + GainMatrix ** ( C1s - C1p ) ** GainMatrix.T()
        
        r_smooth   = Hit - MMMatrix ** x_smooth
        R_smooth   = NMatrix - MMMatrix ** C_smooth ** MMMatrixT
    
        This.Sstate     = KalmanData( x_smooth, C_smooth )
        This.Sresiduals = KalmanData( r_smooth, R_smooth )
    
    def Fit( self ):
        assert self.IState and self.Measurements, '''
            You must set the initial state and the measurements. Use the SetInitialState and SetMeasurements methods.
            '''
        
        for i in range( 1, self.Track.Nnodes ):
            self.Filter( i )
        
        LastNode = self.Track.GetNode(-1)
        LastNode.Sstate = LastNode.Fstate
        LastNode.Sresiduals = LastNode.Fresiduals
        
        for i in reversed(range( self.Track.Nnodes - 1 )):
            self.Smooth( i )

        return self.Track
        





