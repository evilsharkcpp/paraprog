import numpy as np
import scipy as sp
from scipy.sparse import csr_matrix

rng = np.random.default_rng()
rvs = sp.stats.poisson(25, loc=10).rvs
n = np.random.randint(2, 3)
S = sp.sparse.random(n, n, density=0.75, random_state=rng, data_rvs=rvs)
t1 = csr_matrix(S)
S = sp.sparse.random(n, n, density=0.75, random_state=rng, data_rvs=rvs)
t2 = csr_matrix(S)

o1 = open("left.txt", "w")

o1.write(f"{n} {n}\n")
np.savetxt(o1, t1.indptr.astype(int), fmt="%i", newline=" ")
o1.write("\n")
np.savetxt(o1, t1.data, newline=" ")
o1.write("\n")
np.savetxt(o1, t1.indices.astype(int), fmt="%i", newline=" ")
o1.close()
o1 = open("right.txt", "w")

o1.write(f"{n} {n}\n")
np.savetxt(o1, t2.indptr.astype(int), fmt="%i", newline=" ")
o1.write("\n")
np.savetxt(o1, t2.data, newline=" ")
o1.write("\n")
np.savetxt(o1, t2.indices.astype(int), fmt="%i", newline=" ")
o1.close()

add = t1.dot(t2)

o1 = open("mult.txt", "w")

o1.write(f"{n} {n}\n")
np.savetxt(o1, add.indptr.astype(int), fmt="%i", newline=" ")
o1.write("\n")
np.savetxt(o1, add.data, newline=" ")
o1.write("\n")
np.savetxt(o1, add.indices.astype(int), fmt="%i", newline=" ")
o1.close()

print(add.tolil())

sum = t1 + t2
o1 = open("sum.txt", "w")

o1.write(f"{n} {n}\n")
np.savetxt(o1, sum.indptr.astype(int), fmt="%i", newline=" ")
o1.write("\n")
np.savetxt(o1, sum.data, newline=" ")
o1.write("\n")
np.savetxt(o1, sum.indices.astype(int), fmt="%i", newline=" ")
o1.close()
