import ROOT
import shelve
import random
from Array import *
from KalmanFilter import KalmanFilter, KalmanData
from Physics import Constants, NEXT
from KalmanNEXT import NEXTKF, KFNEXT

ifile = shelve.open('../PyKalman/Analysis/singlebeta.shelve')
E0 = 2.5

_R = random.Random()

def single_track():
    from track import track
    x,y,z,u,u,u,e = zip(*track)
    XYEMatrix = Identity(3) * 1.0
    hits = [ KalmanData( Vector(*xye), XYEMatrix) for xye in zip(x,y,e) ]
    kf = NEXTKF()
    kf.SetMeasurements(z,hits)
    kf.SetInitialState( KalmanData( Vector(0.,0.,0.,0.,E0), Identity(5) ) )
    
    track = kf.Fit()
    p     = track.Plot()
    p3    = track.Plot3D()
    cchi  = ROOT.TCanvas()
    chi2  = ROOT.TGraph()
    chi2.SetMarkerStyle(20)
    for i in range(track.Nnodes): chi2.SetPoint(i,i,track.GetNode(i).chi2)
    chi2.Draw('AP')
    raw_input()

def Randomize( x, sigma ):
    for i in range(len(x)):
        x[i] += _R.gauss(0.,sigma)

def full():
    for t in range( ifile['N'] ):
        x,y,z,u,u,u,e = map( list, zip(*ifile[str(t)]) )
        Randomize(x,.1)
        Randomize(y,.1)
        Randomize(z,.1)
#        Randomize(e,.1)
        XYEMatrix = Diagonal( [1., 1., 100.] )
        hits = [ KalmanData( Vector(*xye), XYEMatrix) for xye in zip(x,y,e) ]
#        kf = NEXTKF()
        kf = KFNEXT()
        kf.SetMeasurements(z,hits)
        x0  = 0.
        y0  = 0.
        tx0 = ( x[1] - x[0] ) / ( z[1] - z[0] )
        ty0 = ( y[1] - y[0] ) / ( z[1] - z[0] )
        kf.SetInitialState( KalmanData( Vector(x0,y0,tx0,ty0,E0), Identity(5)*1 ) )
        
        track = kf.Fit()
#        print '\n'*10
#        print 'tnumber',t
#        print '\n'*2
#        print track

        p     = track.Plot()
        p3    = track.Plot3D()
        cchi  = ROOT.TCanvas()
        chi2  = ROOT.TGraph()
        chi2.SetMarkerStyle(20)
        for i in range(track.Nnodes): chi2.SetPoint(i,i,track.GetNode(i).chi2)
        chi2.Draw('AP')
        
        raw_input()
full()
