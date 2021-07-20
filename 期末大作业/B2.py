import numpy as np
from scipy import integrate
import time

start = time.time()

def f(x):
    return np.cos(x) / (1 + np.sin(x)**3)

a = 0
b = 1
res_nc = []
res_gl = []

for n in range(5,21):
    h = (b - a) / n
    coefficient = np.zeros(n+1)
    for k in range(0,n+1):
        coefficient[k] = n * np.math.factorial(k) * np.math.factorial(n-k)
        if (n-k) % 2 != 0:
            coefficient[k] = -coefficient[k]
        order = []
        for j in range(n+1):
            if j != k:
                order.append(j)
        poly = np.poly1d(order, r=True, variable=["x"])
        poly_afterintegral = np.array(np.polyint(poly))
        result = 0
        for j in range(n+1):
            result += poly_afterintegral[j] * n ** (n-j+1)
        coefficient[k] = result / coefficient[k]
    ans_nc = 0
    for k in range(n+1):
        ans_nc += coefficient[k] * f(a+k*h)
    ans_nc *= b - a
    res_nc.append(ans_nc)

T = 20
L = []
L.append(np.poly1d([1], r=False, variable=["x"]))
L.append(np.poly1d([1,0], r=False, variable=["x"]))
for i in range(1,T+1):
    L.append(np.polysub((2*i+1)/(i+1)*np.polymul(L[1],L[i]),i/(i+1)*L[i-1]))

for n in range(5,21):
    Ak = []
    xk = np.roots(L[n+1])
    xk = np.sort(xk)
    for i in range(n+1):
        denominator = 1
        numerator = []
        for j in range(n+1):
            if j != i:
                denominator *= xk[i] -xk[j]
                numerator.append(xk[j])
        poly_numerator = np.poly1d(numerator, r=True, variable=["x"])
        l = poly_numerator / denominator
        l_afterintegral = np.array(np.polyint(l))
        ans = 0
        for j in range(n+2):
            if j % 2 != 0:
                continue
            else:
                ans += 2 * l_afterintegral[j]
        Ak.append(abs(ans))
    ans_gl = 0
    for i in range(n+1):
        ans_gl += Ak[i]*f((b-a)/2*xk[i]+(a+b)/2)
    ans_gl *= (b-a)/2
    res_gl.append(ans_gl)

result, err = integrate.quad(f,a,b)
n = []
quad = []
for i in range(5,21):
    n.append(i)
    quad.append(result)

end1 = time.time()

print("(a) Newton-Cotes")
for i in range(16):
    print("n = {}, the result is {}".format(i+5, res_nc[i]))
print("(b) Gauss-Legendre")
for i in range(16):
    print("n = {}, the result is {}".format(i+5, res_gl[i]))
print("(c) quad()")
print("Using quad(), the result is {}".format(result))

end2 = time.time()
print("\n程序执行时间(不考虑I/O)为：{}".format(end1-start))
print("\n程序执行时间(考虑I/O)为：{}".format(end2-start))