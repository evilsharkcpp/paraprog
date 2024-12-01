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
    rows1, rows2 = np.array(f1.readline()), np.array(f2.readline())
    data1, data2 = np.array(f1.readline()), np.array(f2.readline())
    cols1, cols2 = np.array(f1.readline()), np.array(f2.readline())
    f1.close()
    f2.close()
    A1 = csr_matrix((data1), (rows1, cols1))
    A2 = csr_matrix((data2), (rows2, cols2))
    assert np.allclose(A1, A2)


if __name__ == "__main__":
    main()
