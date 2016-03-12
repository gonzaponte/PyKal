from Array import *
from ROOT import *
from array import array
from operator import itemgetter
from math import *

I3 = Identity(3)

def BuildVector( x ):
    return Vector3(*x)

def RelativeTo( hit0 ):
    return lambda hit: hit - hit0

def Rcut( cut ):
    return lambda x: abs(x[1])<cut

def ComputeProbability( E0 ):
    def x( hit ):
        return abs(Vector3( hit[0], hit[1], 0. )) / MSvariance( E0, abs(hit[2]) )**0.5
    return x

def MSvariance( E, dx, RL = 1.0 ):
    p = ( E**2 + 2*E*0.511 )**0.5
    b = p/E
    return 13.6 / (b*p) * (dx/RL)**0.5 * ( 1 + 0.038 * log(dx/RL) )

def FindPath( initial_hits, initial_energies, seed = [0,1], rcut = 3.0 ):
    seed = sorted( seed, reverse=True )
    
    initial_hits = map( BuildVector, initial_hits )
    
    ordered_hits     = list( reversed( map( initial_hits.pop, seed ) ) )
    ordered_energies = list( reversed( map( initial_energies.pop, seed ) ) )

    while initial_hits:
        E = sum(initial_energies)
        
        zdir = ordered_hits[-1] - ordered_hits[-2]
        rotation_matrix = RotationMatrix( zdir )
        
        indices, relative_hits = zip( *filter( Rcut( rcut ), enumerate( map( RelativeTo( ordered_hits[-1] ), initial_hits ) ) ) )
        rotated_hits    = map( rotation_matrix.__pow__, relative_hits )
        best_hit, index = sorted( zip( map( ComputeProbability(E), rotated_hits ), indices ) )[0]
        
        ordered_hits.append( initial_hits.pop(index) )
        ordered_energies.append( initial_energies.pop(index) )
    
    return ordered_hits, ordered_energies

def RotationMatrix( u, v = Vector3(0,0,1) ):
    '''
        Compute the rotation matrix for translating the vector v to u.
    '''
    u = u.Unit()
    w = v^u
    s = w ** w
    c = v**u
    R = Matrix([0.,-w[2],w[1]],[w[2],0.,-w[0]],[-w[1],w[0],0.])
    return I3 + R + (1.-c)/s * ( R ** R )

