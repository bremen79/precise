import numpy as np
import matplotlib.pyplot as plt

from precise import precise_co96


n=10000;
delta=0.05;

# i.i.d. uniform samples in [0,1]
vector = np.random.uniform(size=n)
print(vector)

lcblist, ucblist = precise_co96(vector, delta)
print("Upper and lower bounds")
print(lcblist, ucblist)

plt.plot(np.arange(1, n + 1), lcblist, label='PRECiSE-CO96', color='red')
plt.plot(np.arange(1, n + 1), ucblist, color='red')
plt.ylim(0, 1)
plt.xlim(1, n)
plt.xscale('log')
plt.xlabel('Number of samples (log scale)')
plt.ylabel('Confidence Sequences')
plt.title(r'Confidence sequences for a uniform r.v. in [0,1], $\delta=0.05$')
plt.legend()
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.show()
