import numpy as np
import sys
sys.path.append('../test')

import build_CS20NG_example_data
import matplotlib.pyplot as plt


data=build_CS20NG_example_data.build_cs20ng_test_data(256)

b = data['b']
b = np.array(b)
print(b.shape)
pix_idx = data['pix_idx']
pix_idx = np.array(pix_idx, dtype=np.int32)

print(pix_idx.shape)

n = data['nrow']
m = data['ncol']

import l1cPy

x,p = l1cPy.l1qc_dct_py(n, m, b, pix_idx, epsilon=0.01, cgmaxiter=200, verbose=2)

print(x[1:20])
print(p)


x = np.reshape(x, (n, m))

plt.figure(1)
plt.imshow(x)

plt.show()

# 2114
