import numpy as np
import matplotlib.pyplot as plt

def is_psd(matrix):
    eigenvalues = np.linalg.eigvalsh(matrix)
    return np.all(eigenvalues >= 0)

def toda(E, ex, size):
    g = 2
    x = np.zeros(2*size - 1)
    x[0] = 1
    x[1] = ex 
    for i in range(2, 2*size - 1):
        x[i] = (2*E*(i - 1)*x[i - 1] - g*(i - 3/2)*x[i - 2])/g/(i - 1/2)
    return np.array([[x[i + j] for j in range(size)] for i in range(size)])

def toda_xp(E, ex, size):
    xlen = max(2*size - 1, 4)
    xp = np.zeros((xlen, 2*size - 1))
    xp[0][0] = 1
    xp[1][0] = ex 
    for i in range(2, xlen):
        xp[i][0] = (2*E*(i - 1)*xp[i - 1][0] + (- 2*i + 3)*xp[i - 2][0])/(2*i - 1)
    for j in range(2, 2*size - 1):
        for i in range(1, 3):
            xp[i][j] = E*xp[i][j - 2] - xp[i + 1][j - 2] - xp[i - 1][j - 2]
        xp[0][j] = ((3 + j)*xp[2][j] - 2*E*xp[1][j])/(j - 1)
        for i in range(3, xlen):
            xp[i][j] = (2*E*(i - 1)*xp[i - 1][j] + (j - 2*i + 3)*xp[i - 2][j])/(2*i + j - 1);
    mat = np.zeros((size*size, size*size))
    for i in range(0, size*size):
        for j in range(i, size*size):
            b = i % size
            a = i // size
            d = j % size
            c = j // size
            mat[i][j] = xp[a + c][b + d]
            mat[j][i] = xp[a + c][b + d]
    return mat, xp

def harmonics(E, size):
    k = 1
    x = np.zeros(2*size - 1)
    x[0] = 1
    x[1] = 0
    for i in range(2, 2*size - 1):
        x[i] = (i - 1)*E*x[i - 2]/i/k
    return np.array([[x[i + j] for j in range(size)] for i in range(size)])

def coulomb(E, size):
    a = 1
    x = np.zeros(2*size - 1)
    x[0] = 1
    x[1] = -1/a - 3*a/4/E 
    for i in range(2, 2*size - 1):
        x[i] = (2*i*x[i - 2] + a*(-1 - 2*n)*x[i - 1])/2/E/(i + 1)
    return np.array([[x[i + j] for j in range(size)] for i in range(size)])

def anharmonics(E, xsq, size):
    g = 1
    x = np.zeros(2*size - 1)
    x[0] = 1
    x[2] = xsq
    for i in range(4, 2*size - 1):
        x[i] = ((i - 2)*x[i - 2] + (i - 3)*x[i - 4]*E)/g/(i - 1)
    return np.array([[x[i + j] for j in range(size)] for i in range(size)])

def anharmonics_xp(E, xsq, size):
    xp = np.zeros((max(2*size - 1, 8), 2*size - 1))
    xp[0][0] = 1
    xp[2][0] = xsq
    for i in range(4, max(2*size - 1, 8)):
        xp[i][0] = ((i - 2)*xp[i - 2][0] + (i - 3)*xp[i - 4][0]*E)/(i - 1)
    for j in range(2, 2*size - 1):
        for i in range(0, 4):
            xp[i][j] = E*xp[i][j - 2] + xp[i + 2][j - 2] - xp[i+4][j - 2]
        for i in range(4, 2*size - 1):
            xp[i][j] = ((i + j - 2)*xp[i - 2][j] + E*(i - 3)*xp[i - 4][j])/(i + 2*j - 1)
    mat = np.zeros((size*size, size*size))
    for i in range(0, size*size):
        for j in range(i, size*size):
            b = i % size
            a = i // size
            d = j % size
            c = j // size
            mat[i][j] = xp[a + c][b + d]
            mat[j][i] = xp[a + c][b + d]
    return mat

def harmonics_xp(E, xp_1_1, size):
    xp = np.zeros((2*size - 1, 2*size  - 1))
    xp[0][0] = 1
    xp[1][1] = xp_1_1
    for j in range(0, 2*size - 3):
        for i in range(2, 2*size - 1):
            xp[i][j] = E*(i - 1)*xp[i - 2][j]/(i + j)
            xp[i - 2][j + 2] = E*xp[i - 2][j] - xp[i][j]
    mat = np.zeros((size*size, size*size))
    for i in range(0, size*size):
        for j in range(0, size*size):
            b = i % size
            a = i // size
            d = j % size
            c = j // size
            mat[i][j] = xp[a + c][b + d]
    return mat

