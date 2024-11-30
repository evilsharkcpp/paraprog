import numpy as np
import scipy as sp
from scipy.sparse import csr_matrix

rng = np.random.default_rng()
rvs = sp.stats.poisson(25, loc=10).rvs
n = np.random.randint(1, 500)
S = sp.sparse.random(n, n, density=0.25, random_state=rng, data_rvs=rvs)
t1 = csr_matrix(S)
S = sp.sparse.random(n, n, density=0.25, random_state=rng, data_rvs=rvs)
t2 = csr_matrix(S)

o1 = open("left.txt", "w")

rows, cols = t1.nonzero()
o1.write(f"{n} {n}\n")
np.savetxt(o1, t1.indptr.astype(int), fmt="%i", newline=" ")
o1.write("\n")
np.savetxt(o1, t1.data, newline=" ")
o1.write("\n")
np.savetxt(o1, cols.astype(int), fmt="%i", newline=" ")

o1 = open("right.txt", "w")

rows, cols = t2.nonzero()
o1.write(f"{n} {n}\n")
np.savetxt(o1, t2.indptr.astype(int), fmt="%i", newline=" ")
o1.write("\n")
np.savetxt(o1, t2.data, newline=" ")
o1.write("\n")
np.savetxt(o1, cols.astype(int), fmt="%i", newline=" ")

add = np.dot(t1, t2)
add = csr_matrix(add)

o1 = open("add.txt", "w")

rows, cols = add.nonzero()
o1.write(f"{n} {n}\n")
np.savetxt(o1, add.indptr.astype(int), fmt="%i", newline=" ")
o1.write("\n")
np.savetxt(o1, add.data, newline=" ")
o1.write("\n")
np.savetxt(o1, cols.astype(int), fmt="%i", newline=" ")

sum = np.sum(t1, t2)
o1 = open("sum.txt", "w")

rows, cols = sum.nonzero()
o1.write(f"{n} {n}\n")
np.savetxt(o1, sum.indptr.astype(int), fmt="%i", newline=" ")
o1.write("\n")
np.savetxt(o1, sum.data, newline=" ")
o1.write("\n")
np.savetxt(o1, cols.astype(int), fmt="%i", newline=" ")
