import numpy as np

from precise import precise_co96


n=10000;
delta=0.05;

# i.i.d. uniform samples in [0,1]
vector = np.random.uniform(size=n)
print(vector)

lcblist, ucblist = precise_co96(vector, delta)
print("Upper and lower bounds")
print(lcblist, ucblist)
