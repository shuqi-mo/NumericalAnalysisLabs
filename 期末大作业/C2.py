import numpy as np
import time

start = time.time()

def LUDec(A):  
    n = len(A)
    L = np.zeros(shape=(n,n))
    U = np.zeros(shape=(n,n))  
    for base in range(n-1):
        for i in range(base+1,n):
            L[i,base]=A[i,base]/A[base,base]
            A[i]=A[i]-L[i,base]*A[base]
    for i in range(n):
        L[i,i]=1
    U=np.array(A)
    return L,U

def solve(L,U,b):
    rows = len(b)
    y = np.zeros(rows)
    y[0] = b[0]/L[0,0]
    for k in range(1,rows):
        y[k] = (b[k] - np.sum(L[k,:k]*y[:k]))/L[k,k]     
    x = np.zeros(rows)
    k = rows-1
    x[k] = y[k]/U[k,k] 
    for k in range(rows-2,-1,-1):
        x[k] = (y[k] - np.sum(x[k+1:]*U[k,k+1:]))/U[k,k]    
    return x
    
A = np.array( [[ 20.,   2.,   3.,   0.],
               [  1.,   8.,   1.,   1.],
               [  2.,  -3.,  15.,   0.],
               [  1.,   0.,   0.,   1.]],dtype='float')
b1 = np.array([[ 1., 0., 0., 0.]],dtype='float').T
b2 = np.array([[ 0., 1., 0., 0.]],dtype='float').T
b3 = np.array([[ 0., 0., 1., 0.]],dtype='float').T
b4 = np.array([[ 0., 0., 0., 1.]],dtype='float').T
inverse = np.linalg.inv(A)
La,Ua = LUDec(A)

end1 = time.time()

print("(a)")
print("LU分解得到的L为:")
print(La)
print("LU分解得到的U为:")
print(Ua)
print("(b)")
print("解得X1,X2,X3,X4的值分别为:")
x1 = solve(La,Ua,b1)     
print("X1=",x1)
x2 = solve(La,Ua,b2)      
print("X2=",x2)
x3 = solve(La,Ua,b3)      
print("X3=",x3)
x4 = solve(La,Ua,b4)      
print("X4=",x4)
print("矩阵A的逆为:")
print(inverse)

end2 = time.time()

print("\n程序执行时间(不考虑I/O)为：{}".format(end1-start))
print("\n程序执行时间(考虑I/O)为：{}".format(end2-start))
