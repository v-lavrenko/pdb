#!/usr/bin/env python3

# setup:
# sudo apt install python3
# python3 -m pip install numpy matplotlib

from matplotlib import pyplot as plt
from numpy import random, array, sign, dot, diag, sin, cos
from numpy.linalg import svd, det
import sys

def kabsch(A,B):              # A,B [n x 3]   ... source / target point clouds
    n,w = A.shape
    H = dot(A.T,B)            # H = A' x B    ... covariance O(3*3*N)
    U,s,V = svd(H)            # H = U,s,V'    ... singular-value decomposition
    U,V = U.T,V.T
    d = sign(det(dot(V,U)))   # d = sgn |VU'| ... left/right handed coords
    D = diag([1]*(w-1) + [d]) # D = [1 1 d]
    R = dot(dot(V,D),U)       # R = V D U'  ... optimal rotation
    return R

# toy example of Kabsch

# A = 6 arbitrary points in the XY plane
A = array([[ 1, 1], 
           [-1,-1],
           [.5, 0],
           [-1,.5],
           [ 0,-1],
           [ 0, 0]])

r = 0.5 # arbitrary rotation angle (radians)

# R = rotation matrix
R = array([[ cos(r),sin(r)], 
           [-sin(r),cos(r)]])

# B = rotated(A) + random noise
B = dot(A,R.T) + 0.1 * random.rand(*A.shape)

# Q = Kabsch rotation matrix from A to B
Q = kabsch(A,B)

# C = attempt to recover B using Q
C = dot(A,Q.T)

print (R)
print (Q)
plt.scatter(B[:,0], B[:,1], marker='o', label='target  B')
plt.scatter(A[:,0], A[:,1], marker='v', label='source  A')
plt.scatter(C[:,0], C[:,1], marker='^', label='rotated A')
plt.grid()
plt.legend()
plt.show()
