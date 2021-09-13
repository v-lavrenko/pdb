#!/usr/bin/python3

from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
from numpy import sign, dot, diag, ones, sin, cos
from numpy.linalg import svd, det, norm
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
import sys

def arg(i): return sys.argv[i] if i < len(sys.argv) else ''

def argp(key): return key in sys.argv

def getarg(key,default):
    for i in range(len(sys.argv)):
        if sys.argv[i] == key: return sys.argv[i+1]
    return default

def getprm(prm,key,default=''):
    if key not in prm: return default
    beg = prm.index(key) + len(key)
    return prm[beg:].partition(',')[0]

# https://github.com/ericjang/svd3 ... 3x3 SVD in C: 2μs
# https://en.wikipedia.org/wiki/Kabsch_algorithm

def kabsch(A,B):              # A,B [n x 3]   ... source / target point clouds
    n,w = A.shape
    H = dot(A.T,B)            # H = A' x B    ... covariance O(3*3*N)
    U,s,V = svd(H)            # H = U,s,V'    ... singular-value decomposition
    U,V = U.T,V.T
    d = sign(det(dot(V,U)))   # d = sgn |VU'| ... left/right handed coords
    D = diag([1]*(w-1) + [d]) # D = [1 1 d]
    R = dot(dot(V,D),U)       # R = V D U'  ... optimal rotation
    return R

def toy_kabsch():
    A = np.array([[ 1, 1],
                  [-1,-1],
                  [.5, 0],
                  [-1,.5],
                  [ 0,-1],
                  [ 0, 0]])
    r = 0.5 # rotation angle (radians)
    R = np.array([[ cos(r),sin(r)],
                  [-sin(r),cos(r)]])
    B = dot(A,R.T) + 0.1 * np.random.rand(*A.shape)
    Q = kabsch(A,B)
    C = dot(A,Q.T)
    print (R)
    print (Q)
    plt.scatter(B[:,0], B[:,1], marker='o', label='target  B')
    plt.scatter(A[:,0], A[:,1], marker='v', label='source  A')
    plt.scatter(C[:,0], C[:,1], marker='^', label='rotated A')
    plt.grid()
    plt.legend()
    plt.show()

def read_aligned():
    A,B,a,b = [],[],[],[]
    for line in sys.stdin:
        if line[:2] == 'AB':
            X = [float(x) for x in line.split()[3:]]
            A += [X[:3]]
            B += [X[3:]]
        elif line[0] == 'A':
            a += [[float(x) for x in line.split()[3:]]]
        elif line[0] == 'B':
            b += [[float(x) for x in line.split()[3:]]]
    return np.array(A),np.array(B),np.array(a),np.array(b)

def fit(A,B): return norm(A-B,axis=1).mean().round(3)

def scatter3d(XX):
    markers = 'x^ov*'
    colors = 'red blue green magenta brown'.split()
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for i,X in enumerate(XX):
        ax.scatter(X[:,0],X[:,1],X[:,2],marker=markers[i%5],color=colors[i%5])
    plt.show()

def quiver3d(A,B):
    D = B-A
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.quiver(A[:,0],A[:,1],A[:,2],D[:,0],D[:,1],D[:,2])
    plt.show()
        
def df2xyz(AA,I=None): # dataframe of atoms -> Ix3 matrix (index -> x,y,z)
    return np.array(AA.iloc[AA.index if I is None else list(I),[3,4,5]])

# returns [(a,b)]: AA[a] is best match to BB[b] based on radius ||XYZ||
def align_radius(AA,BB,eps=0.1,far=0.2): 
    AB = [] 
    for el in set(AA.el): # el = H,C,N,CA,...
        iA = AA.index[AA.el == el] # [2,18,42] ... indices of 'CA' atoms in A
        iB = BB.index[BB.el == el] #
        if not len(iB): continue # no el in B -> nothing to align
        LA = np.array(AA.iloc[iA].L2) # [0.0, 3.8, 5.8] ... L2 of 'CA' atoms in A
        LB = np.array(BB.iloc[iB].L2) #
        D = np.array([[abs(L2a-L2b) for L2b in LB] for L2a in LA])
        for a,b in zip(*linear_sum_assignment(D)): # best match
            b2 = sorted(D[a,:])[1] if len(LB)>1 else far # 2nd best in row a
            a2 = sorted(D[:,b])[1] if len(LA)>1 else far # 2nd best in col b
            if D[a,b] <= eps and b2 >= far and a2 >= far:
                AB += [(iA[a],iB[b])]
    return AB

# returns [(a,b)]: AA[a] is best match to BB[b] based on actual distance
def align_nearest(AA,BB):
    AB = [] 
    for el in set(AA.el): # el = H,C,N,CA,...
        iA = AA.index[AA.el == el] # [2,18,42] ... indices of 'CA' atoms in A
        iB = BB.index[BB.el == el] #
        if not len(iB): continue
        XA,XB = df2xyz(AA,iA),df2xyz(BB,iB) # [x,y,z] for each 'CA' atom in A,B
        D = cdist(XA,XB)                    # A x B distance matrix
        ab = zip(*linear_sum_assignment(D)) # (a,b) ... minimal-distance alignment
        AB += [(iA[a],iB[b],D[a,b]) for a,b in ab] # convert to dataframe indices
    As,Bs,Ds = zip(*AB)
    AB += [(a,None,AA.iloc[a].L2) for a in set(AA.index) - set(As)] # unmatched in A
    AB += [(None,b,BB.iloc[b].L2) for b in set(BB.index) - set(Bs)] # unmatched in B
    return AB

def matched(AB): return [(a,b,d) for a,b,d in AB if a is not None and b is not None]
def onlyA(AB): return [(a,d) for a,b,d in AB if b is None]
def onlyB(AB): return [(b,d) for a,b,d in AB if a is None]

def dist(A,B): return np.sqrt((A.x-B.x)**2 + (A.y-B.y)**2 + (A.z-B.z)**2)

def print_aligned(AB,AA,BB):
    for a,b in AB:
        A,B = AA.iloc[a],BB.iloc[b]
        d0,d1 = abs(A.L2-B.L2),dist(A,B)
        print('AB:%s %.3f %.3f %+.3f,%+.3f,%+.3f %+.3f,%+.3f,%+.3f'.replace(' ','\t') %
              (A.el,d0,d1,A.x,A.y,A.z,B.x,B.y,B.z))

def print_unaligned(AB,AA,BB,eps=0.5):
    for a,b in AB:
        A,B = AA.iloc[a],BB.iloc[b]
        if dist(A,B) <= eps: continue
        print('AB:%s %.3f %+.3f,%+.3f,%+.3f %+.3f,%+.3f,%+.3f'.replace(' ','\t') %
              (A.el,dist(A,B),A.x,A.y,A.z,B.x,B.y,B.z))
    for a in set(set(AA.index) - set([a for a,b in AB])):
        A = AA.iloc[a]
        print('A:%s %.3f %+.3f,%+.3f,%+.3f'.replace(' ','\t') % (A.el,A.L2,A.x,A.y,A.z))
    for b in set(set(BB.index) - set([b for a,b in AB])):
        B = BB.iloc[b]        
        print('B:%s %.3f %+.3f,%+.3f,%+.3f'.replace(' ','\t') % (B.el,B.L2,B.x,B.y,B.z))

def F1(AB,A,B): return np.round(2*len(AB) / (len(A) + len(B)),2)
def FN(AB,A,B): return np.round(1 - len(AB)/len(A), 2)
def FP(AB,A,B): return np.round(1 - len(AB)/len(B), 2)

def eFN(AB,A,B,eps=0.1): # Miss rate for epsilon-insensitive matching atoms: A -> B
    matches = [(a,b,d) for a,b,d in matched(AB) if d <= eps]
    return 1 - len(matches)/len(A)
    #print(len(matches), len(A), len(B), len(AB))

def EL2(AB,A,B): # mean distance from atom in A to corresponding atom in B
    ab = [d for a,b,d in matched(AB)] 
    a = [d for a,d in onlyA(AB)] # distance to origin for unmatched atoms
    b = [d for b,d in onlyB(AB)]
    return sum(ab+a)/len(ab+a) # ignore unmatched in B
    #return sum(ab+a+b)/len(ab+a+b)
    #print('%.2f/%d %.2f/%d %.2f/%d' % (sum(ab), len(ab), sum(a), len(a), sum(b), len(b)))

def rotate_xyz(AA,R):
    AX = np.array(AA.iloc[:,[3,4,5]])
    AA.iloc[:,[3,4,5]] = dot(AX,R.T)

def center(AA):
    AX = np.array(AA.iloc[:,[3,4,5]])
    AA.iloc[:,[3,4,5]] = AX - AX.mean(axis=0)        

#         id   el  asn      x      y      z     L2	sim
#0    484D:2    N    1  0.699  1.261 -0.392  1.494    1.000
#1    484D:2   CA    3  0.000  0.000  0.000  0.000    0.612
def do_eval(TSV, eps=.1, far=.2, thr=.8):
    DD = pd.read_csv(TSV, sep='\t') # all atoms in all balls
    if not len(DD):
        sys.stderr.write('ERROR: %s is empty, exiting\n' % (TSV))
        return
    qid = DD.iloc[0].id # query id = id of 1st atom in TSV
    Q = DD[DD.id == qid].reset_index(drop=True) # Q = all atoms for qid
    print("query_id document_id Edist FN0.1 FN0.2 FN0.5 sim".replace(' ','\t'))
    for id in DD.id.unique():
        D = DD[DD.id == id].reset_index(drop=True) # document D
        QD = align_radius(Q,D,eps,far) # rough-align subset of Q,D based on L2 distances
        QX = df2xyz(Q,[q for q,d in QD])
        DX = df2xyz(D,[d for q,d in QD])
        RX = kabsch(DX,QX) # optimal rotation R: D -> Q
        rotate_xyz(D,RX) # apply rotation to all atoms in D
        A1 = align_nearest(Q,D) # alignment after rotation
        L2,FN1,FN2,FN5 = EL2(A1,Q,D),eFN(A1,Q,D,0.1),eFN(A1,Q,D,0.2),eFN(A1,Q,D,0.5)
        print ("%s %s %.3f %.3f %.3f %.3f %.3f".replace(' ','\t') %
               (qid, id, L2, FN1, FN2, FN5, D.iloc[0].sim))
        sys.stdout.flush()

def do_align(_A,_B,eps=.1,far=.5):     #         id   el  asn      x      y      z     L2
    if False:
        AA = pd.read_csv(_A, sep='\t') #0    484D:2    N    1  0.699  1.261 -0.392  1.494
        BB = pd.read_csv(_B, sep='\t') #1    484D:2   CA    3  0.000  0.000  0.000  0.000
    else:
        All = pd.read_csv(_A, sep='\t')
        qid = All.iloc[0].id
        AA = All[All.id==qid].reset_index(drop=True)
        BB = All[All.id!=qid].reset_index(drop=True)
    AB = align_radius(AA,BB,eps,far) # roughly align subset of A,B based on L2 distances
    ma,mb = zip(*AB) # matched indices for A,B
    oa,ob = set(AA.index)-set(ma),set(BB.index)-set(mb) # unmatched in A,B
    AX,BX = df2xyz(AA,ma),df2xyz(BB,mb)
    R = kabsch(AX,BX) # learn rotation
    if False:
        CX = dot(AX,R.T)
        scatter3d([AX,CX,BX])
        print_aligned(AB,AA,BB)
        print ('aligned:', fit(AX,BX), '->', fit(CX,BX), 'matched:', len(ma), 'todo:', len(oa), len(ob))
    rotate_xyz(AA,R) # apply rotation
    AB = align_nearest(AA,BB)
    ma,mb = zip(*AB) # matched indices in A,B
    oa,ob = set(AA.index)-set(ma),set(BB.index)-set(mb) # unmatched in A,B
    AX,BX = df2xyz(AA,ma),df2xyz(BB,mb)
    A0,B0 = df2xyz(AA,oa),df2xyz(BB,ob)
    scatter3d([AX,BX,A0,B0])
    #print_aligned(AB,AA,BB)
    print ('aligned:', fit(AX,BX), 'matched:', len(ma), 'todo:', len(oa), len(ob))
    print_unaligned(AB,AA,BB)
    print (eval_alignment(AB,AA,BB))
    #
    #scatter3d([C,B])
        
def do_topKnn(_A):
    AA = pd.read_csv(_A, sep='\t')
    AX = df2xyz(AA)
    D = cdist(AX,AX)
    for i in AA.index:
        R = sorted(D[i])
        m = '\t'.join(['%.3f' % R[r] for r in [1,2,5,10]])
        print ('%s\t%s' % (AA.iloc[i].el, m))
    
usage = '''
Usage:
pdb.py toy        ... toy example of Kabsch algorithm
pdb.py eval BALLS ... 
'''

def show_usage(): print (usage); sys.exit(1)
    
if   arg(1) == 'toy': toy_kabsch()
elif arg(1) == 'eval': do_eval(arg(2))
elif arg(1) == 'align': do_align(arg(2),arg(3))
elif arg(1) == 'topK': do_topKnn(arg(2))
else: show_usage()


# H = np.dot(A.T,B)

# U,s,V = svd(H)

# S = np.diag(s)

# print(U,S,V)

# err = H - np.dot(np.dot(U,S),V) # H == U x S x V

# print (np.round(err,10))

# d = det(np.dot(V.T,U.T))

# D = np.diag([1,np.sign(d)])

# Q = np.dot(np.dot(V.T,D),U.T)

# def do_eval0():
#     A,B,a,b = read_aligned()
#     R = kabsch(A,B)
#     C = dot(A,R.T)
#     print ('aligned:', fit(A,B), '->', fit(C,B))
#     print (R.round(3))
#     scatter3d([A,C,B])
#     scatter3d([A,B,a,b])
#     cfit(a,b)
#     cfit(dot(a,R.T),b)
#     #print (A.shape,B.shape,a.shape,b.shape,R.shape)

# def cfit(A,B):
#     D = cdist(A,B)
#     ri, ci = linear_sum_assignment(D)
#     print (D.round(2))
#     print (ri)
#     print (ci)
#     print (D[ri,ci].mean())

# def eval_alignment(AB,AA,BB):
#     recall = len(AB)/len(AA) # recall: what % of the query atoms were found in the doc
#     precis = len(AB)/len(BB) # precision: what % of the doc is made up of query atoms
#     F1 = 2 * recall * precis / (recall + precis)
#     return {'recall':np.round(recall,2),
#             'precision':np.round(precis,2),
#             'F1':np.round(F1,2)}
