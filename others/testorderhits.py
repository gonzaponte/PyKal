from orderhits import *
from Array import *
from ROOT import *
from Plots import *
from array import array


def toarray( x ):
    return array( 'f', x )

def BuildLine( x0, x1, color = kBlack, style = 1, width = 2 ):
    line = TPolyLine3D( 2, *map( toarray, zip(x0,x1) ) )
    line.SetLineColor( color )
    line.SetLineWidth( width )
    line.SetLineStyle( style )
    return line

def BuildMarkers( hits ):
    p = TPolyMarker3D(len(hits),20)
    for i,hit in enumerate(hits):
        p.SetPoint(i,*hit)
    return p

def rotation_test():
    o = Vector3(0,0,0)
    u = Vector3(1.,1.,1.).Unit()
    R = RotationMatrix(u)

    x0 = Vector3(1,0,0)
    y0 = Vector3(0,1,0)
    z0 = Vector3(0,0,1)
    x1 = R ** x0
    y1 = R ** y0
    z1 = R ** z0


    h = TH3F('h','',2,-1,1,2,-1,1,2,-1,1)
    uline  = BuildLine( o, u, kViolet )
    x0line = BuildLine( o, x0, kBlack )
    y0line = BuildLine( o, y0, kRed )
    z0line = BuildLine( o, z0, kBlue )
    x1line = BuildLine( o, x1, kBlack, 2 )
    y1line = BuildLine( o, y1, kRed, 2 )
    z1line = BuildLine( o, z1, kBlue, 2 )

    h.Draw()
    uline.Draw('same')
    x0line.Draw('same')
    y0line.Draw('same')
    z0line.Draw('same')
    x1line.Draw('same')
    y1line.Draw('same')
    z1line.Draw('same')

    print x1 ** y1
    print x1 ** z1
    print z1 ** y1
    raw_input('ok')

def GetEvent(filename):
    data = open(filename).readlines()
    x,y,zi,zf,u,u,u,e,de,dx = map( list, zip(*[ map( float, line ) for line in map( str.split, data[1:] ) ] ) )
    z = [ (zii+zff)*0.5 for zii,zff in zip(zi,zf) ]
    return map( Vector, x, y, z ), de

def ordering_test():
    hits, energies = GetEvent('singlebeta_0.dat')
    h0, e0 = FindPath( hits, energies )
#    plt = Plot4D( *zip(*h0), t = range(len(h0)) )
    marks = BuildMarkers( h0 )
    line = TPolyLine3D( len(h0), *map( toarray, zip(*h0) ) )
    marks.Draw()
    line.SetLineColor(kRed)
    line.SetLineWidth(2)
    line.Draw('same')
    raw_input()

ordering_test()
