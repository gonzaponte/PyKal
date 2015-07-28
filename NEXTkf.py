import ROOT
import shelve
import random
from Array import *
from KalmanFilter import KalmanFilter, KalmanData
from Physics import Constants, NEXT
from KalmanNEXT import NEXTKF, KFNEXT

#ifile = shelve.open('../PyKalman/Analysis/singlebeta.shelve')
ifile = shelve.open('../PyKalman/Analysis/singlebeta_voxelized.shelve')
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

def Randomize( x, sigma, minval = -float('inf'), maxval = float('inf') ):
    
    if not isinstance(sigma, list):
        sigma = [ sigma ] * len(x)
    for i in range(len(x)):
        smear = _R.gauss(0.,sigma[i])
        if minval < x[i] + smear < maxval:
            x[i] += smear
        else:
            x[i] = minval if smear < 0 else maxval

def full():
    for t in range( ifile['N'] ):
        x,y,z,u,u,u,e = map( list, zip(*ifile[str(t)]) )
        ### mm -> cm
        x = [ xi*1e-1 for xi in x ]
        y = [ yi*1e-1 for yi in y ]
        z = [ zi*1e-1 for zi in z ]
        
        ### randomize hits according to resolution
#        Randomize(x,1e-2)
#        Randomize(y,1e-2)
#        Randomize(e,[0.01*(2.458/ei)**0.5 for ei in e], 1e-8)

        xres = 0.1**2
        eres = lambda e: 0.01*(2.458/e)**0.5
        hits = [ KalmanData( Vector(*xye), Diagonal( [ xres, xres, eres(xye[-1])] )) for xye in zip(x,y,e) ]
#        for hit in hits: print hit
#        kf = NEXTKF()
        kf = KFNEXT()
        kf.SetMeasurements(z,hits)
        x0  = x[0]
        y0  = y[0]
        tx0 = ( x[1] - x[0] ) / ( z[1] - z[0] )
        ty0 = ( y[1] - y[0] ) / ( z[1] - z[0] )
        e0  = e[0]
        kf.SetInitialState( KalmanData( Vector(x0,y0,tx0,ty0,e0), Identity(5)*10 ) )
        
        track = kf.Fit()
        
        p     = track.Plot()
        p3    = track.Plot3D()
        cchi  = ROOT.TCanvas()
        chi2  = ROOT.TGraph()
        chi2.SetMarkerStyle(20)
        for i in range(track.Nnodes): chi2.SetPoint(i,i,track.GetNode(i).chi2)
        chi2.Draw('AP')
        
        raw_input()
full()
