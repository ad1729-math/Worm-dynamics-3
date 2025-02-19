#The dual dimer graph for a square lattice on a torus is considered. 
import numpy as np 
import matplotlib.pyplot as plt 
import networkx as nx
from matplotlib.widgets import Button
from matplotlib.backend_bases import MouseEvent
import matplotlib.pyplot as plt
from pyvis.network import Network
from numpy import random
import time 
import math
from collections import Counter
import csv

st=time.time()
b0=np.log(1+np.sqrt(2))/2

N,M=101,101 

def ad_x(x):
   if x==1:
     return [N,2]
   elif x==N:
     return [N-1, 1]
   else:
     return [x-1,x+1]
   
def ad_y(x):
   if x==1:
     return [M,2]
   elif x==M:
     return [M-1, 1]
   else:
     return [x-1,x+1]
   
A,Kas=[],[]

#II-Defining dual decorated lattice for torus

for y in range(1,M+1): #Here, m=0, M+1 are boundaries of S1 X R  #Checked to be completely correct

    for x in range(1,N+1):
      if y==1: 

        A.extend([[[x,y,2,1],[ad_x(x)[0],y,1,1]],[[x,y,1,1],[ad_x(x)[1],y,2,1]],[[x,y,1,2],[x,ad_y(y)[0],0,2]],[[x,y,2,2],[x,ad_y(y)[1],1,2]]])

      elif y==M:

        A.extend([[[x,y,2,1],[ad_x(x)[0],y,1,1]],[[x,y,1,1],[ad_x(x)[1],y,2,1]],[[x,y,2,2],[x,ad_y(y)[1],0,2]],[[x,y,1,2],[x,ad_y(y)[0],2,2]]])

      else:
        A.extend([[[x,y,2,1],[ad_x(x)[0],y,1,1]],[[x,y,1,1],[ad_x(x)[1],y,2,1]],[[x,y,1,2],[x,ad_y(y)[0],2,2]],[[x,y,2,2],[x,ad_y(y)[1],1,2]]])

    for x in range(1,N+1):

      for k in range(1,3):
        A.extend([[[x,y,k,1],[x,y,k,2]],[[x,y,k,2],[x,y,k,3]],[[x,y,k,3],[x,y,k,1]]])
        A.extend([[[x,y,k,2],[x,y,k,1]],[[x,y,k,1],[x,y,k,3]],[[x,y,k,3],[x,y,k,2]]])

      A.extend([[[x,y,1,3],[x,y,2,3]],[[x,y,2,3],[x,y,1,3]]])


def Enum(x,y,k,l):
  return 6*(N*(y-1)+(x-1))+3*(k-1)+l


def Rev_Enum(c):

  m1=math.ceil(c/(6*N))
  n1=math.ceil((c-(6*(m1-1))*N)/6)
  x,y=n1,m1 
  k=math.ceil((c-(6*N*(m1-1)+6*(n1-1)))/3)
  l=c-(6*N*(m1-1)+6*(n1-1))-3*(k-1)

  return [x,y,k,l]

p0=1  #Probability of nature of the bond, ferromagnetic or anti-ferromagnetic 


def ProbB(p):
  choices = [1, -1]
  return np.random.choice(choices, p=[p, 1-p])

J0=1
J_st=0 #J_st=0 means square and 1 is triangular
h=0 #External field is zero

def Prob(v,J1):
  return np.random.normal(v,J1)

def Inter(b):
  Interact, J_list=[], []

  for y in range(1,M+1):
    for x in range(1,N+1):
      j=ProbB(p0)
      Interact.append([[x,y],[ad_x(x)[1],y],j])
      J_list.append(j)
      
      j1=J_st*ProbB(p0) #Diagonal edge making square to triangular
      Interact.append([[x,y],[x,y],j1])

      if 1<=y<=M:
        j=ProbB(p0) 
        Interact.append([[x,y],[x,ad_y(y)[1]],j])
        J_list.append(j)

  j_max=np.sort(J_list)[-1]

  return Interact, j_max

Int_list=Inter(1)
Int=Int_list[0]
j0=Int_list[1]

#The weights 
def W(r1,r2,b):

  vh,vh1, vh2=np.exp(b*h), np.exp(-b*h), np.exp(-b*h/2)
  
  c1, c2=Rev_Enum(r1), Rev_Enum(r2)
  x1,y1,k1,l1=c1[0],c1[1],c1[2],c1[3]
  x2,y2,k2,l2=c2[0],c2[1],c2[2],c2[3]
  
  if y2==ad_y(y1)[0]:
    if x1==x2 and [k1,k2,l1,l2]==[1,2,2,2]:
      for e in Int:
        if e[0]==[x1,y1] and e[1]==[ad_x(x1)[1],y1]: #inner loop connenting to the 2nd loop
          v=np.exp(e[2]*b) 

    else: v=0

  elif y2==ad_y(y1)[1]:
    if x1==x2 and [k1,k2,l1,l2]==[2,1,2,2]:
      for e in Int:
        if e[0]==[x1,ad_y(y1)[1]] and e[1]==[ad_x(x1)[1],ad_y(y1)[1]]: #inner loop connenting to the 2nd loop
          v=np.exp(e[2]*b) 
          
    else: v=0

  else: 

    if y2==y1:
      if x2==ad_x(x1)[0] and [k1,k2,l1,l2]==[2,1,1,1]:
        for e in Int:
          if e[0]==[x1,y1] and e[1]==[x1,ad_y(y1)[1]]: #inner loop connenting to the 2nd loop
            v=np.exp(e[2]*b) 
      
      elif x2==ad_x(x1)[1] and [k1,k2,l1,l2]==[1,2,1,1]:

        for e in Int:
          if e[0]==[ad_x(x1)[1],y1] and e[1]==[ad_x(x1)[1],ad_y(y1)[1]]: #inner loop connenting to the 2nd loop
            v=np.exp(e[2]*b) 

      else:
        if x1==x2:

          if k1==1:
            if k2==2:
              if [l1,l2]==[3,3]:
                for e in Int:
                  if e[0]==[x1,y1] and e[1]==[x1,y1]:
                    v=np.exp(e[2]*b)
              else: v=0

            elif k2==k1:
              if [l1,l2]==[3,2] or [l1,l2]==[2,3]:
                for e in Int:
                  if e[0]==[x1,y1] and e[1]==[ad_x(x1)[1],y1]: 
                    v1=np.exp(-e[2]*b/2)
                  elif e[0]==[x1,y1] and e[1]==[x1,y1]:
                    v2=np.exp(-e[2]*b/2)
                v=v1*v2

              elif [l1,l2]==[3,1] or [l1,l2]==[1,3]:
                for e in Int:
                  if e[0]==[ad_x(x1)[1],y1] and e[1]==[ad_x(x1)[1],ad_y(y1)[1]]: 
                    v1=np.exp(-e[2]*b/2)
                  elif e[0]==[x1,y1] and e[1]==[x1,y1]:
                    v2=np.exp(-e[2]*b/2)
                v=v1*v2

              else:
                if [l1,l2]==[1,2] or [l1,l2]==[2,1]:
                  for e in Int: 
                    if e[0]==[x1,y1] and e[1]==[ad_x(x1)[1],y1]: 
                      v1=np.exp(-e[2]*b/2)
                    elif e[0]==[ad_x(x1)[1],y1] and e[1]==[ad_x(x1)[1],ad_y(y1)[1]]: 
                      v2=np.exp(-e[2]*b/2)

                  v=v1*v2
                else: v=0
                  
          elif k1==2:

            if k2==1:
              if [l1,l2]==[3,3]:
                for e in Int:
                  if e[0]==[x1,y1] and e[1]==[x1,y1]:
                    v=np.exp(e[2]*b)
              else: v=0

            elif k2==k1:
              if [l1,l2]==[3,2] or [l1,l2]==[2,3]:
                for e in Int:
                  if e[0]==[x1,ad_y(y1)[1]] and e[1]==[ad_x(x1)[1],ad_y(y1)[1]]: 
                    v1=np.exp(-e[2]*b/2)
                  elif e[0]==[x1,y1] and e[1]==[x1,y1]:
                    v2=np.exp(-e[2]*b/2)
                v=v1*v2

              elif [l1,l2]==[3,1] or [l1,l2]==[1,3]:
                for e in Int:
                  if e[0]==[x1,y1] and e[1]==[x1,ad_y(y1)[1]]: 
                    v1=np.exp(-e[2]*b/2)
                  elif e[0]==[x1,y1] and e[1]==[x1,y1]:
                    v2=np.exp(-e[2]*b/2)
                v=v1*v2

              else:
                if [l1,l2]==[1,2] or [l1,l2]==[2,1]:

                  for e in Int: 
                    if e[0]==[x1,y1] and e[1]==[x1,ad_y(y1)[1]]: 
                      v1=np.exp(-e[2]*b/2)

                    elif e[0]==[x1,ad_y(y1)[1]] and e[1]==[ad_x(x1)[1],ad_y(y1)[1]]: 
                      v2=np.exp(-e[2]*b/2)

                  v=v1*v2     

                else: v=0 

        else: v=0

    else: v=0

  return v

#Defining neighbouts of each vertex

Neighbours=[]
for y in range(1,M+1):
  for x in range(1,N+1):
    c11=Enum(x,y,1,1)
    Neighbours.append([c11, Enum(x,y,1,2),Enum(x,y,1,3),Enum(ad_x(x)[1],y,2,1)])

    c12=Enum(x,y,1,2)
    Neighbours.append([c12, Enum(x,y,1,1),Enum(x,y,1,3),Enum(x,ad_y(y)[0],2,2)])

    c13=Enum(x,y,1,3)
    Neighbours.append([c13, Enum(x,y,1,1),Enum(x,y,1,2),Enum(x,y,2,3)])

    c21=Enum(x,y,2,1)
    Neighbours.append([c21, Enum(x,y,2,2),Enum(x,y,2,3),Enum(ad_x(x)[0],y,1,1)])

    c22=Enum(x,y,2,2)
    Neighbours.append([c22, Enum(x,y,2,1),Enum(x,y,2,3),Enum(x,ad_y(y)[1],1,2)])

    c23=Enum(x,y,2,3)
    Neighbours.append([c23, Enum(x,y,2,1),Enum(x,y,2,2),Enum(x,y,1,3)])


def Neigh(c):
  L=Neighbours[c-1]
  return [L[1],L[2],L[3]]

# st1=time.time()
# r=Enum(10,16,1,2)
# for i in range(2000):
#   print(Neigh(r+i))
# et1=time.time()
# print(et1-st1)

#Initial unique domain matching
Matching_0=[]

for x in range(1,N+1):
  for y in range(1,M+1):

    cx1,cx_2=Enum(x,y,1,1),Enum(ad_x(x)[1],y,2,1)
    cx2,cx_1=Enum(x,y,2,1),Enum(ad_x(x)[0],y,1,1)

    c1,c2=Enum(x,y,1,3),Enum(x,y,2,3)

    Matching_0.extend([[cx1,cx_2],[cx_2,cx1],[c1,c2],[c2,c1]])

  for y in range(1,M+1): 
    cy1,cy_2=Enum(x,y,1,2),Enum(x,ad_y(y)[0],2,2)
    cy2,cy_1=Enum(x,y,2,2),Enum(x,ad_y(y)[1],1,2)
    
    Matching_0.extend([[cy1,cy_2],[cy2,cy_1]])

def Dimer_identifier(a, Matching_0):
  for s in Matching_0:
    if s[0]==a:
      return s[1]
      break
    elif s[1]==a:
      return s[0]
      break

#Once the detailed balance matrix is constructed, there is no need for the weights, this significantly imporves the speed.
cmax=Enum(N,M,2,3)
T=[0,1,2]
A=[] #Detailed balance matrix at each point
for c in range(1,cmax+1):
  Ac=[[0,0,0],[0,0,0],[0,0,0]]
  L_=Neigh(c)

  n1,n2,n3=L_[0],L_[1],L_[2] #Neighbours to the vertex c
  w1,w2,w3=W(c,n1,b0),W(c,n2,b0),W(c,n3,b0) #Weights corresponding to those
  W_l=[w1,w2,w3]
  Ws=np.sort(W_l) #Sorted the weights

  ws1,ws2,wsm=Ws[0],Ws[1],Ws[2]

  im=W_l.index(wsm)  #This is the neighbour with maximum weight
  T0=[]
  for t in T:
    if t!=im:
      T0.append(t)
      
  i1,i2=T0[0],T0[1]

  if wsm>ws1+ws2:
    Ac[im][im]=wsm-(ws1+ws2)
    Ac[i1][im],Ac[im][i1]=ws1,ws1
    Ac[i2][im],Ac[im][i2]=ws2,ws2
  else:
    Ac[i1][i2],Ac[i2][i1]=(ws1+ws2-wsm)/2,(ws1+ws2-wsm)/2
    Ac[i1][im],Ac[im][i1]=(wsm+ws2-ws1)/2,(wsm+ws2-ws1)/2
    Ac[i2][im],Ac[im][i2]=(ws1+wsm-ws2)/2,(ws1+wsm-ws2)/2

  A.append([Ac[0][0],Ac[0][1],Ac[0][2],Ac[1][0],Ac[1][1],Ac[1][2],Ac[2][0],Ac[2][1],Ac[2][2],w1,w2,w3])

with open('NEIG.csv', 'w', newline='') as file:
  writer = csv.writer(file)
  writer.writerows(Neighbours)

with open('Mat.csv', 'w', newline='') as file:
  writer = csv.writer(file)
  writer.writerows(Matching_0)

with open('Det_balance.csv', 'w', newline='') as file:
  writer = csv.writer(file)
  writer.writerows(A)


#Dimer worm dynamics
def Dimer_move(x1,x2): #Here, x2 is the pivot and x1 is a dimer
  A1=A[x2-1] #This list contains the detailed balance matrix as well as the edge weigths
  NB=Neighbours[x2-1] #This list contains the neighbours of the pivot x2
  NB_=[NB[1],NB[2],NB[3]]
  ind=NB_.index(x1)
  P=[A1[(ind)*3]/A1[9+ind],A1[(ind)*3+1]/A1[9+ind],A1[(ind)*3+2]/A1[9+ind]] #The transition matrix 

  r=np.random.rand()

  if r<=P[0]:
    return NB_[0]
  elif P[0]<r<=P[0]+P[1]:
    return NB_[1]
  else:
    return NB_[2]

# r=500
# En=10000
# N_=Neigh(r)
# st=time.time()
# for i in range(En):
#   print(Dimer_move(N_[0],r))
# et=time.time()

# print(et-st)

#Performing dimer worm dynamics 

# runs=100
# x10,x20=1201,Dimer_identifier(1201,Matching_0) #x1 is the worm head and x2 is the initial pivot

# X10,Y10=Rev_Enum(x10)[0],Rev_Enum(x10)[1]
# X20,Y20=Rev_Enum(x20)[0],Rev_Enum(x20)[1]
  
# Ensemble=100
# D=[]

# st=time.time()

# Length, Distance=[],[]
# for en in range(Ensemble):
#   length=0
#   x1,x2=x10,x20

#   Matching=Matching_0.copy()

#   for i in range(runs):
#     c1=Dimer_move(x1,x2) #New dimer worm head 
    
#     Matching=[x for x in Matching if x not in [[x1,x2], [x2,x1]]]

#     if c1==x10:
#       c2=x2

#     else:
#       c2=Dimer_identifier(c1,Matching) #New dimer pivot 
#       Matching.extend([[x2,c1],[c1,x2]])

#       x1,x2=c1,c2 

#       length+=1

#   X1,Y1=Rev_Enum(x1)[0],Rev_Enum(x1)[1]
#   X2,Y2=Rev_Enum(x2)[0],Rev_Enum(x2)[1]

#   d=abs(X1-X10)+abs(X2-Y20)

#   # print(d)

#   # Length.append(length) #length of the dimer loop
#   D.append(d) #Distance of the final worm head
  
# et=time.time()

# print(et-st)

# L=np.sort(Length)
# D0=np.sort(D)

# print(D)

# def count_elements(sorted_list, N):
#     # Count occurrences of each element in the list
#     element_count = Counter(sorted_list)
    
#     # Create a list to hold counts from 1 to N
#     result = []
#     for i in range(1, N + 1):
#         result.append([i, element_count.get(i, 0)])  # Get the count, default to 0 if not found
    
#     return result

# Count_length=count_elements(L, runs)
# Count_dist=count_elements(D, runs)

# Len, C_len, C_dist=[],[],[]
# for i in range(runs):
#   Len.append(Count_length[i][0])
#   C_len.append(Count_length[i][1])
#   C_dist.append(Count_dist[i][1]/Ensemble)

# print(C_len)

