import numpy as np

from precise import precise_co96


#vector = np.random.binomial(n=1, p=0.5, size=100)
vector = np.random.uniform(size=10)
print(vector)

lcblist, ucblist = precise_co96(vector, 0.05)
print("Upper and lower bounds")
print(lcblist, ucblist)
