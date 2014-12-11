import timeit # module with timing subroutines                                                              
import random # module to generate random numnbers                                                          
from itertools import imap,starmap,izip
from operator import mul

def v(N=3,min=-10,max=10):
    """Generates a random vector (in an array) of dimension N; the                                          
    values are integers in the range [min,max]."""
    out = []
    for k in range(N):
        out.append(random.randint(min,max))
    return out

def check(v1,v2):
    if len(v1)!=len(v2):
        raise ValueError,"the lenght of both arrays must be the same"
    pass

def d0(v1,v2):
    """                                                                                                     
    d0 is Nominal approach:                                                                                 
    multiply/add in a loop                                                                                  
    """
    check(v1,v2)
    out = 0
    for k in range(len(v1)):
        out += v1[k] * v2[k]
    return out

def d1(v1,v2):
    """                                                                                                     
    d1 uses an imap (from itertools)                                                                        
    """
    check(v1,v2)
    return sum(imap(mul,v1,v2))

def d2(v1,v2):
    """                                                                                                     
    d2 uses a conventional map                                                                              
    """
    check(v1,v2)
    return sum(map(mul,v1,v2))

def d3(v1,v2):
    """                                                                                                     
    d3 uses a starmap (itertools) to apply the mul operator on an izipped (v1,v2)                           
    """
    check(v1,v2)
    return sum(starmap(mul,izip(v1,v2)))

# generate the test vectors                                                                                 
v1 = v()
v2 = v()

if __name__ == '__main__':

    # Generate two test vectors of dimension N                                                              

    t0 = timeit.Timer("d0(v1,v2)","from dot_product import d0,v1,v2")
    t1 = timeit.Timer("d1(v1,v2)","from dot_product import d1,v1,v2")
    t2 = timeit.Timer("d2(v1,v2)","from dot_product import d2,v1,v2")
    t3 = timeit.Timer("d3(v1,v2)","from dot_product import d3,v1,v2")

    print "d0 elapsed: ", t0.timeit()
    print "d1 elapsed: ", t1.timeit()
    print "d2 elapsed: ", t2.timeit()
    print "d3 elapsed: ", t3.timeit()
