#!/usr/bin/env python
# coding: utf-8
#You need Sage

#Code provided by Rachelle-Heim Boissier and Yann Rotella
#Everything is free of use, re-use and change
#Don't hesitate to contact us if you find flaws
#Or more importantly if you find ways to improve the attack !!!

w = 4
l = 2
r = 30
c = 70
b = 100
from brial import *
B = declare_ring([Block('m', r), Block('s', c), Block('h',1)],globals())

def creation_tab_3dim(p,n,m):
    tab = [0]*p
    for i in range (p):
        tab[i] = [0]*n
    for i in range (p):
        for j in range (n):
            tab[i][j] = [0]*m
    return tab

##Create a vector containing Ti's polynomials for further evaluation
def creation_tab_perm(A):
    tabA = [0]*b
    for y in range(5):
        for x in range(5):
            for z in range(w):
                tabA[w*(5*y +x) +z] = A[x][y][z]
    return tabA

def eval_poly(A,tab):
    M = A
    for i in range(b):
        #if not (tab[i] == B.gens()[200]):
        M = M.subs({B.gens()[i] : tab[i]})
    return M

T = creation_tab_3dim(5,5,w)

#theta
for x in range (5):
    for y in range (5):
        for z in range(w):
            T[x][y][z] = B.gens()[w*(5*y + x) + z]
            for p in range (5):
                T[x][y][z] = T[x][y][z] + B.gens()[w*(5*p + (x-1)%5) + z]
                T[x][y][z] = T[x][y][z] + B.gens()[w*(5*p + (x+1)%5) + (z-1)%w]

tabT = creation_tab_perm(T)


#Permutation $\rho$
R = creation_tab_3dim(5,5,w)
for z in range(w):
        R[0][0][z] = B.gens()[z]
x = 1;
y = 0;
temp = 0; 
for t in range(24):
    for z in range(w):
        R[x][y][z] = B.gens()[   w*(5*y + x) + ((z-(int((t+1)*(t+2)/2)))%w) ]
    temp = y;
    y = (2*x + 3*y)%5;
    x = temp;

#Permutation $\rho \circ \theta $
RT = creation_tab_3dim(5,5,w)
for y in range(5):
    for x in range(5):
        for z in range(w):
            RT[x][y][z] = substitute_variables(B, tabT, R[x][y][z])

tabRT = creation_tab_perm(RT)

#Permutation $\pi$
P = creation_tab_3dim(5,5,w)
for x in range (5):
    for y in range (5):
        for z in range(w):
            P[x][y][z] = B.gens()[ w*(5*x + ((x+3*y)%5)) + z ]

#Permutation $\pi \circ \rho \circ \theta $ 
Lin = creation_tab_3dim(5,5,w)
for y in range(5):
    for x in range(5):
        for z in range(w):
            Lin[x][y][z] = substitute_variables(B, tabRT, P[x][y][z])

tabLin = creation_tab_perm(Lin)

#Permutation $\chi$
K = creation_tab_3dim(5,5,w)
for x in range(5):
    for y in range(5):
        for z in range(w):
            K[x][y][z] = B.gens()[w*(5*y + x) +z] 
            K[x][y][z] += ( B.gens()[w*(5*y + ((x+1)%5) ) +z] +1)*(B.gens()[w*(5*y +((x+2)%5)) +z])

#Permutation $\chi \circ \pi \circ \rho \circ \theta $ 
NotConst = creation_tab_3dim(5,5,w)
for y in range(5):
    for x in range(5):
        for z in range(w):
            NotConst[x][y][z] = substitute_variables(B, tabLin, K[x][y][z])


#After $\chi$
tabLin2 = creation_tab_perm(Lin)

##print('Slice 0 : ')
##print('z = 0, x = 0, y = 0,1,2,3,4\n')
z = 0
x = 0
for y in range(0,5): 
    tabLin2[w*(5*y + x) +z] = 0*B.gens()[1]
    
##print('\nSlice 1 :')
##print('z = 1, x = 1, y = 3,4\n')
z = 1
x = 1
for y in range(3,5): 
    tabLin2[w*(5*y + x) +z] = 0*B.gens()[1]
    
##print('\nSlice 2 :')
##print('z = 2, x = 2, y = 0,1,2,3,4\n')
z = 2
x = 2
for y in range(0,5): 
    tabLin2[w*(5*y + x) +z] = 0*B.gens()[1]

##print('\nSlice 3 :')
##print('z = 3, x = 3, y = 0,1,2,3,4\n')
z = 3
x = 3
for y in range(0,5): 
    tabLin2[w*(5*y + x) +z] = 0*B.gens()[1]

AfterAlloc = creation_tab_3dim(5,5,w)
for y in range(5):
    for x in range(5):
        for z in range(w):
            AfterAlloc[x][y][z] = substitute_variables(B, tabLin2, K[x][y][z])


#Linear system building
LinSys = [0]*30

i = 0

##print('Slice 0 : ')
##print('z = 0, x = 0, y = 0,1,2,3,4\n')
z = 0
x = 0
for y in range(0,5): 
    LinSys[i] = tabLin[w*(5*y + x) +z]
    i = i+1
    
##print('\nSlice 1 :')
##print('z = 1, x = 1, y = 3,4\n')
z = 1
x = 1
for y in range(3,5): 
    LinSys[i] = tabLin[w*(5*y + x) +z]
    i = i+1
    
    
##print('\nSlice 2 :')
##print('z = 2, x = 2, y = 0,1,2,3,4\n')
z = 2
x = 2
for y in range(0,5): 
    LinSys[i] = tabLin[w*(5*y + x) +z]
    i = i+1
    

##print('\nSlice 3 :')
##print('z = 3, x = 3, y = 0,1,2,3,4\n')
z = 3
x = 3
for y in range(0,5): 
    LinSys[i] = tabLin[w*(5*y + x) +z]
    i = i+1

##print('Slice 0 : ')
##print('z = 0, x = 4, y = 0,2,3\n')
z = 0
x = 4
LinSys[i] = AfterAlloc[x][0][z] + AfterAlloc[x][2][z]
i += 1
LinSys[i] = AfterAlloc[x][2][z] + AfterAlloc[x][3][z]
i = i+1

##print('\n')
    
##print('z = 0, x = 3, y = 1,2,4\n')
z = 0
x = 3
LinSys[i] = AfterAlloc[x][1][z] + AfterAlloc[x][2][z]
i = i+1
LinSys[i] = AfterAlloc[x][2][z] + AfterAlloc[x][4][z]
i = i+1 

##print('Slice 3 : ')
##print('z = 3, x = 1 , y = 0,2,4\n')
z = 3
x = 1
LinSys[i] = AfterAlloc[x][0][z] + AfterAlloc[x][2][z]
i += 1
LinSys[i] = AfterAlloc[x][2][z] + AfterAlloc[x][4][z]
i += 1

##print('z = 3, x = 2, y = 0,1,3\n')
z = 3
x = 2
LinSys[i] = AfterAlloc[x][0][z] + AfterAlloc[x][1][z]
i = i+1
LinSys[i] = AfterAlloc[x][1][z] + AfterAlloc[x][3][z]
i = i+1


##print('Slice 2 : ')
##print('z = 2, x = 1, y = 0,2,4\n')
z = 2
x = 1
LinSys[i] = AfterAlloc[x][0][z] + AfterAlloc[x][2][z]
i += 1
LinSys[i] = AfterAlloc[x][2][z] + AfterAlloc[x][4][z]
i += 1

##print('z = 2, x = 0, y = 1,3,4\n')
z = 2
x = 0
LinSys[i] = AfterAlloc[x][1][z] + AfterAlloc[x][3][z]
i = i+1
LinSys[i] = AfterAlloc[x][3][z] + AfterAlloc[x][4][z]
i = i+1 

##print('Slice 1 : ')


##print('z = 1, x = 0, y = 3,4\n')
z = 1
x = 0
LinSys[i] = AfterAlloc[x][3][z] + AfterAlloc[x][4][z]

LinSysM = [0]*30
LinSysS = [0]*30
tabM = []
for i in range(30):
    tabM.append(B.gens()[i])
for i in range(70):
    tabM.append(B.gens()[i+30]*0)


for i in range(30):
    LinSysM[i] = substitute_variables(B, tabM, LinSys[i])
LinSysM


# In[26]:


for i in range(30):
    LinSysS[i] = LinSys[i] + LinSysM[i]
LinSysS


#Gaussian 

#number of lines n, number of columns (variables) m
n = 30
m = 30
tabSysM = [0]*n
for i in range(n):
    tabSysM[i] = [0]*m
tabSysS = [0]*n
for i in range(n):
    tabSysS[i] = [0]*70

MatriceM = matrix(GF(2),n,m)
for i in range(n):
    for j in range(m):
        if B.gens()[j] in LinSysM[i]:
            MatriceM[i,j] = 1
print(MatriceM.rank())
print(MatriceM)
MatriceS = matrix(GF(2),30,70)
for i in range(30):
    for j in range(70):
        if B.gens()[j+30] in LinSysS[i]:
            MatriceS[i,j] = 1

A = MatriceM.extended_echelon_form()
U = MatriceM.echelon_form()
C = A.submatrix(0,30)
MatriceM == C**(-1)*U
NewMatriceS = C*MatriceS
print(NewMatriceS)
print('\n')
print(U)

for i in range(21,30):
    print('i =',i, '   ',U[i])

finalTab = matrix(GF(2),33,70)
for i in range(22):
    for j in range(70):
        finalTab[i,j] = NewMatriceS[i,j]

for j in range(70):
    finalTab[24,j] = NewMatriceS[24-2,j]

for i in range(26,33):
    for j in range(70):
        finalTab[i,j] = NewMatriceS[i-3,j]
print(finalTab)

fichier = open("finalTab.txt", 'w')
fichier.write(str(finalTab))
fichier.close()