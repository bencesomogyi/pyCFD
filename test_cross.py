import numpy
import time
import math

a=numpy.array([1.,0.,0.])
b=numpy.array([0.,1.,0.])

t0=time.time()
c=numpy.cross(a,b)
t1=time.time()

t2=time.time()
x = ((a[1] * b[2]) - (a[2] * b[1]))
y = ((a[2] * b[0]) - (a[0] * b[2]))
z = ((a[0] * b[1]) - (a[1] * b[0]))
d=numpy.array([x,y,z])
t3=time.time()

print c
print t1-t0
print d
print t3-t2

print "numpy.add"
t0=time.time()
c=numpy.add(a,b)
t1=time.time()
t2=time.time()
x = a[0] + b[0]
y = a[1] + b[1]
z = a[2] + b[2]
d=numpy.array([x,y,z])
t3=time.time()

print c
print t1-t0
print d
print t3-t2

print "\nnumpy norm"
a=numpy.random.rand(3)
t4=time.time()
e=numpy.linalg.norm(a)
t5=time.time()
t6=time.time()
f=math.sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
t7=time.time()

print e
print t1-t0
print ""
print f
print t3-t2
