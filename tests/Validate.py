import sys
import numpy as np
import scipy as sp
from scipy.sparse import csr_matrix


def main():
    ref = sys.argv[1]
    data = sys.argv[2]
    f1 = open(ref)
    f2 = open(data)
    sizes1 = f1.readline()
    sizes2 = f2.readline()
    rows1, rows2 = np.array(
        list(f1.readline().strip().split(" ")), dtype=int
    ), np.array(list(f2.readline().strip().split(" ")), dtype=int)
    data1, data2 = np.array(
        list(f1.readline().strip().split(" ")), dtype=float
    ), np.array(list(f2.readline().strip().split(" ")), dtype=float)
    cols1, cols2 = np.array(
        list(f1.readline().strip().split(" ")), dtype=int
    ), np.array(list(f2.readline().strip().split(" ")), dtype=int)
    f1.close()
    f2.close()
    A1 = csr_matrix((data1, cols1, rows1))
    A2 = csr_matrix((data2, cols2, rows2))
    assert (A1 != A2).nnz == 0


if __name__ == "__main__":
    main()
