import numpy

a=numpy.array([0.1,0.2,0.3,0.4,0.5])
i=0
while i < len(a):
    print "\ncurrent: "+str(i)+" "+str(a[i])
    print a
    if i == 0:
        a=numpy.insert(a,1,0.15)
        print a
    i += 1
    