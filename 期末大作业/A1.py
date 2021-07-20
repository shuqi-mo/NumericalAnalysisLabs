import numpy as np
import time

start = time.time()

def lagrange(x):
    y_lr = 0
    for i in range(n):
        denominator = 1.0
        numerator = 1.0
        for j in range(n):
            if i == j:
                continue
            denominator *= x_sample[i] - x_sample[j]
            numerator *= x - x_sample[j]
        base = numerator / denominator
        y_lr += base * y_sample[i]
    return y_lr
    
def newton(x):
    y_newton = diff[1][0]
    delta_x = x - x_sample[0]
    for i in range(1,n):
        y_newton += diff[i+1][i] * delta_x
        delta_x *= x - x_sample[i]
    return y_newton

start1 = time.time()
x_sample = np.array([0.4,0.55,0.65,0.8,0.9,1.05])
y_sample = np.array([0.41075,0.57825,0.69675,0.88811,1.02652,1.25382])
n = x_sample.shape[0]
# 均差表
diff = np.zeros((n+1,n))
for i in range(n):
    diff[0][i] = x_sample[i]
    diff[1][i] = y_sample[i]
for diff_index in range(1,n):
    for i in range(diff_index,n):
        diff[diff_index+1][i] = (diff[diff_index][i]-diff[diff_index][i-1])/(diff[0][i]-diff[0][i-diff_index])
end1 = time.time()

print("(a)")
print("拉格朗日插值法的结果：")
print("f(0.42) = {}".format(lagrange(0.42)))
print("f(0.596) = {}".format(lagrange(0.596)))
print("f(1.0) = {}".format(lagrange(1.0)))
print("牛顿插值法的结果：")
print("f(0.42) = {}".format(newton(0.42)))
print("f(0.596) = {}".format(newton(0.596)))
print("f(1.0) = {}".format(newton(1.0)))

start2 = time.time()
y_sample = np.array([0.4,0.55,0.65,0.8,0.9,1.05])
x_sample = np.array([0.41075,0.57825,0.69675,0.88811,1.02652,1.25382])
n = x_sample.shape[0]
# 均差表
diff = np.zeros((n+1,n))
for i in range(n):
    diff[0][i] = x_sample[i]
    diff[1][i] = y_sample[i]
for diff_index in range(1,n):
    for i in range(diff_index,n):
        diff[diff_index+1][i] = (diff[diff_index][i]-diff[diff_index][i-1])/(diff[0][i]-diff[0][i-diff_index])
end2 = time.time()

print("(b)")
print("拉格朗日插值法的结果：")
print("f(z1) = 0.5, z1 = {}".format(lagrange(0.5)))
print("f(z2) = 0.75, z2 = {}".format(lagrange(0.75)))
print("f(z3) = 1.0, z3 = {}".format(lagrange(1.0)))
print("牛顿插值法的结果：")
print("f(z1) = 0.5, z1 = {}".format(newton(0.5)))
print("f(z2) = 0.75, z2 = {}".format(newton(0.75)))
print("f(z3) = 1.0, z3 = {}".format(newton(1.0)))
end = time.time()

print("\n程序执行时间(不考虑I/O)为：{}".format(end2-start2+end1-start1))
print("\n程序执行时间(考虑I/O)为：{}".format(end-start))