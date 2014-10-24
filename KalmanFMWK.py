'''
    Kalman stuff
'''
import copy
import Array

class KalmanMeasurement:
    '''
        Represents a state of the Kalman Filter. It is a pair (vector, covariance matrix).
    '''

    def __init__( self, vector = Array.Vector(), CovarianceMatrix = Array.Matrix() ):
        ''' 
            Initialize with a state vector and its covariance matrix.
        '''
        self.Vector = vector
        self.CovarianceMatrix = CovarianceMatrix
    
    def __len__( self ):
        '''
            Return the dimension of the state vector.
        '''
        return len(self.Vector)

    def __str__( self ):
        '''
            String representation for printing purposes.
        '''
        return '''
Vector = {0}
Cov =
{1}
'''.format( self.Vector, self.CovarianceMatrix )

class KalmanNode:
    '''
        Contains all the information about a given step. Contains the index of the step, the hit, the predicted, filtered and smoothed states and residuals as well as the chi2 and the cumulative chi2 of the step.
    '''
    
    def __init__( self, step = 0, running = 0, hit = KalmanMeasurement(),
                  pred_state = KalmanMeasurement(), filt_state = KalmanMeasurement(), smooth_state = KalmanMeasurement(),
                  pred_resid = KalmanMeasurement(), filt_resid = KalmanMeasurement(), smooth_resid = KalmanMeasurement(),
                  chi2       = -1, cumchi2 = -1 ):

        self.step         = step
        self.running      = running
        self.hit          = hit
        self.pred_state   = pred_state
        self.filt_state   = filt_state
        self.smooth_state = smooth_state
        self.pred_resid   = pred_resid
        self.filt_resid   = filt_resid
        self.smooth_resid = smooth_resid
        self.chi2         = chi2
        self.cumchi2      = cumchi2

    def __str__( self ):
        '''
            String representation for printing purposes.
        '''
        return '''
step number {0}
running variable value {1}
hit:
{2}
predicted state:
{3}
filtered state:
{4}
smoothed state:
{5}
predicted residuals:
{6}
filtered residuals:
{7}
smoothed residuals:
{8}
chi2 of this node:
{9}
cumulative chi2 at this node:
{10}
'''.format( self.step, self.hit, self.running,
            self.pred_state, self.filt_state, self.smooth_state,
            self.pred_resid, self.filt_resid, self.smooth_resid,
            self.chi2, self.cumchi2 )

class KalmanTrack:
    '''
        A set of KalmanNodes. Contains all possible information about the track.
    '''
    def __init__( self, nodes = list() ):
        '''
            Initialize with the list of nodes. It is also possible to add them later with the AddNode method.
        '''
        self.nodes = copy.deepcopy(nodes)
        self.nnodes = len( self.nodes )

    def AddNode( self, node ):
        '''
            Add a node.
        '''
        self.nodes.append( copy.copy( node ) )
        self.nnodes += 1
    
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
        
        xmeasurements.SetMarkerStyle(20)
        ymeasurements.SetMarkerStyle(20)
        xfiltered.SetMarkerStyle(29)
        yfiltered.SetMarkerStyle(29)

        xmeasurements.SetMarkerColor(1)
        ymeasurements.SetMarkerColor(1)
        xfiltered.SetMarkerColor(3)
        yfiltered.SetMarkerColor(3)

        for node in self.nodes[:-1]:
            print node.step, node.hit.Vector[0], node.filt_state.Vector[0], node.hit.Vector[1], node.filt_state.Vector[1]
            xmeasurements.SetPoint( node.step, node.running, node.hit.Vector[0] )
            ymeasurements.SetPoint( node.step, node.running, node.hit.Vector[1] )
            xfiltered    .SetPoint( node.step, node.running, node.filt_state.Vector[0] )
            yfiltered    .SetPoint( node.step, node.running, node.filt_state.Vector[1] )

        cx = ROOT.TCanvas()
        xmeasurements.Draw('AP')
        xfiltered.Draw('P')

        cy = ROOT.TCanvas()
        ymeasurements.Draw('AP')
        yfiltered.Draw('P')

        return cx, xmeasurements, xfiltered, cy, ymeasurements, yfiltered
        
    def __str__( self ):
        '''
            String representation for printing purposes.
        '''
        return '\nNumber of steps: {0} \n\n'.format( self.nnodes ) + '\n\n'.join( map( str, self.nodes ) )

class KalmanFilter:
    '''
        Abstract implementation of the Kalman Filter.
    '''
    
    def __init__( self, name = 'KalmanFilter' ):
        '''
            Initializer. It is needed to call the SetInitialState and SetMeasurements methods mandatorily.
        '''
        self.name             = name
        self.Track            = KalmanTrack()
        self.HaveInitialState = False # To ensure all the information is present when fitting.
        self.HaveMeasurements = False # To ensure all the information is present when fitting.
        self.guess            = None

    def SetInitialState( self, state = KalmanMeasurement() ):
        '''
            Sets "state" as the initial state of the Kalman filter, both as the predicted state and the filtered state.
        '''
        first_node            = self.Track.GetNode(0)
        first_node.pred_state = state
        first_node.filt_state = state
        first_node.pred_resid = self.MeasurementMatrix(0) * 0.
        first_node.filt_resid = self.MeasurementMatrix(0) * 0.
        first_node.chi2       = 0.
        first_node.cumchi2    = 0.
        self.state_dim        = len(state)
        self.measurement_dim  = len(first_node.hit.Vector)
        self.HaveInitialState = True
    
    def SetInitialGuess( self, guess = None ):
        '''
            Sets "guess" as the initial guess for fitting. If no guess is given an internal algorithm is used.'
        '''
        self.guess = guess if guess else self._ComputeInitialGuess()
    
    def SetMeasurements( self, runs, hits ):
        '''
            Sets "hits" as the measurements of the track.
        '''
        indices = range(len(runs))
        for i,r,hit in zip( indices, runs, hits ):
            self.Track.AddNode( KalmanNode( step = i, running = r, hit = hit) )
        self.HaveMeasurements = True # Now we can compute stuff

    def _ComputeInitialGuess( self ):
        #To be implemented
        self.guess = None

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

    def _NewNode( self, index ):
        '''
            Shorter names are defined in order to make the code more readable.
        '''
        
        self.prev_node  = self.Track.GetNode( index - 1 )
        self.this_node  = self.Track.GetNode( index )
        self.next_node  = self.Track.GetNode( index + 1 )
        self.this_hit   = self.this_node.hit.Vector
        self.prev_state = self.prev_node.filt_state.Vector
        self.prev_cov   = self.prev_node.filt_state.CovarianceMatrix
        self.MSMatrix   = self.MultipleScatteringMatrix( index )
        self.MMMatrix   = self.MeasurementMatrix( index )
        self.MMMatrixT  = self.MMMatrix.T()
        self.TMatrix    = self.TransportMatrix( index )
        self.TMatrixT   = self.TMatrix.T()
        self.NMatrix    = self.NoiseMatrix( index )
        self.NMatrixI   = self.NMatrix.Inverse()


    def Predict( self, index ):
        '''
            Predict the value for step "index".
        '''
        if not index:
            raise ValueError( 'The initial step cannot be predicted' )

        x_predicted = self.TMatrix ** self.prev_state
        C_predicted = self.TMatrix ** self.prev_cov ** self.TMatrixT + self.MSMatrix
        r_predicted = self.this_hit - self.MMMatrix ** x_predicted
        R_predicted = self.NMatrix + self.MMMatrix ** C_predicted ** self.MMMatrixT
    
        self.this_node.pred_state = KalmanMeasurement( x_predicted, C_predicted )
        self.this_node.pred_resid = KalmanMeasurement( r_predicted, R_predicted )
        
    def Filter( self, index ):
        '''
            Filter the index-th step.
        '''
        C_filtered = ( self.this_node.pred_state.CovarianceMatrix.Inverse() + self.MMMatrixT ** self.NMatrixI ** self.MMMatrix ).Inverse()
        GainMatrix = C_filtered ** self.MMMatrixT ** self.NMatrixI
        x_filtered = self.prev_state + GainMatrix ** ( self.this_hit - self.MMMatrix ** self.prev_state )
        
        projector  = Array.Identity( self.measurement_dim ) - self.MMMatrix ** GainMatrix
        r_filtered = projector ** self.this_node.pred_resid.Vector
        R_filtered = projector ** self.NMatrix
        
        chi2plus = r_filtered ** R_filtered.Inverse() ** r_filtered
        newchi2  = self.prev_node.cumchi2 + chi2plus
        
        self.this_node.filt_state = KalmanMeasurement( x_filtered, C_filtered )
        self.this_node.filt_resid = KalmanMeasurement( r_filtered, R_filtered )
        self.this_node.chi2       = chi2plus
        self.this_node.cumchi2    = newchi2
    
    def Smooth( self, index ):
        return
        #Still not working
        GainMatrix = self.this_node.filt_state.CovarianceMatrix ** self.TMatrixT ** self.next_node.pred_state.CovarianceMatrix.Inverse()
        x_smooth   = self.this_node.filt_state.Vector + GainMatrix# ** ( ??? )
    
    def Fit( self ):
        assert self.HaveInitialState and self.HaveMeasurements, '''
            You must set the initial state and the measurements. Use the SetInitialState and SetMeasurements methods.
            '''
        if self.guess is None:
            self._ComputeInitialGuess()
        
        for i in range( 1, self.Track.nnodes - 1 ):
            self._NewNode( i )
            self.Predict( i )
            self.Filter( i )
        for i in reversed(range( 1, self.Track.nnodes - 1 )):
            self._NewNode( i )
            self.Smooth( i )

        return self.Track
        





