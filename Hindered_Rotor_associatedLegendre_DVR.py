"""
Hamiltonian
H = Bj^2 + V(cos(theta))

Basis
\ket{j, m}

ポテンシャルV(cos(theta))に対してDVRを適用．解析解(eigenvalue: analytical solution by spherical harmonics basis)と解が一致することを確認できる．
x(= cos(theta))行列を対角化することでグリッドおよびスペクトラル基底<->グリッド基底のユニタリー行列を計算する．

DVRは基本的にライブラリを使った方が手っ取り早いが
・Legendreの陪関数についてのGauss求積ライブラリが見つからなかったこと
・他の言語での実装のこと
を考えて自分で計算するプログラムを実装した．

"""
import numpy as np
from sympy.physics.wigner import wigner_3j

#matrix size
l = 5
dim = 0
for i in range (0, l + 1):
    dim += 2*i + 1
dimension = dim

#potential paramter
V1 = 10
l_pot = 3

#molecular paramter
B = 10


#検算のためにLegendre多項式を書いているが任意関数について計算可能  
def potential_function(x):  #x = cos(theta)
    #V = x 　　#l_pot = 1
    #V = (1/2)*(3*x**2 - 1)  #l_pot = 2
    V = (1/2)*(5*x**3 - 3*x)  #l_pot = 3
    return V

"""球面調和関数基底のDVR計算ここから"""
#球面調和関数のDVR計算のためのクラス
class MatrixElementSH:
    def __init__(self, lrow, mrow, lcolumn, mcolumn):
        self.lrow = lrow
        self.mrow = mrow
        self.lcolumn = lcolumn
        self.mcolumn = mcolumn
    
    def symbol(self, p, q):
        symbol1 = wigner_3j( self.lrow, p, self.lcolumn,
                            -self.mrow, q, self.mcolumn)
        
        symbol2 = wigner_3j( self.lrow,  p, self.lcolumn,
                                     0, 0,             0)
        symbol = symbol1*symbol2
        return symbol
    
    def Vr(self, p, q):
        Vr = ((-1)**(self.mrow))*np.sqrt((2*self.lcolumn + 1)*(2*self.lrow + 1))*self.symbol(p, q)
        return Vr
    
    def x(self):
        Vr1 = self.Vr(1, 0)
        return Vr1

def DVR_spherical_harmonics():
    x = np.array([[+ MatrixElementSH(lrow, mrow, lcolumn, mcolumn).x()
                   for lcolumn in range (0, l + 1)
                   for mcolumn in range (-lcolumn, lcolumn + 1)
                   ]
            for lrow in range (0, l + 1)
            for mrow in range (-lrow, lrow + 1)
            ])
    x = x.astype(np.float64)
    Xi, Xivec = np.linalg.eigh(x)
    return Xi, Xivec
Xi = DVR_spherical_harmonics()[0]
Xivec = DVR_spherical_harmonics()[1]
"""球面調和関数基底のDVR計算ここまで"""


#Hamiltonian行列要素を計算するためのクラス
class MatrixElement:
    def __init__(self, lrow, mrow, lcolumn, mcolumn):
        self.lrow = lrow
        self.mrow = mrow
        self.lcolumn = lcolumn
        self.mcolumn = mcolumn
        
    def row(self):
        row = 0
        for i in range (0, self.lrow):
            row += 2*i + 1
        return int(row + self.lrow + self.mrow)
    
    def column(self):
        column = 0
        for i in range (0, self.lcolumn):
            column += 2*i + 1
        return int(column + self.lcolumn + self.mcolumn)
    
    def symbol(self, p, q):
        symbol1 = wigner_3j( self.lrow, p, self.lcolumn,
                            -self.mrow, q, self.mcolumn)
        
        symbol2 = wigner_3j( self.lrow,  p, self.lcolumn,
                                     0, 0,             0)
        symbol = symbol1*symbol2
        return symbol
    
    def Tr(self):
        if self.lrow == self.lcolumn and self.mrow == self.mcolumn:
            Tr = self.lrow*(self.lrow + 1)
        else: 
            Tr = 0
        return Tr
    
    def Vr(self, p, q):
        Vr = ((-1)**(self.mrow))*np.sqrt((2*self.lcolumn + 1)*(2*self.lrow + 1))*self.symbol(p, q)
        return Vr
    
    def Vr1_SH(self):
        Vr1 = self.Vr(l_pot, 0)
        return Vr1
    
    def Vr1_DVR(self, Ri, Ri_vec):
        Vr1 = 0
        for i in range (0, dimension):
            Vr1 += potential_function(Ri[i])*Ri_vec[self.row(), i]*Ri_vec[self.column(), i]
        if self.mrow == self.mcolumn:
            V =  Vr1
        else:
            V = 0
        return V
    
    
def Hamiltonian_SH():
    result = np.array([[+ MatrixElement(lrow, mrow, lcolumn, mcolumn).Tr() * B
                        + MatrixElement(lrow, mrow, lcolumn, mcolumn).Vr1_SH() * V1
                        for lcolumn in range (0, l + 1)
                        for mcolumn in range (-lcolumn, lcolumn + 1)
                        ]
                    for lrow in range (0, l + 1)
                    for mrow in range (-lrow, lrow + 1)
                    ])
    result = result.astype(np.float64)
    return result

def Hamiltonian_DVR():
    result = np.array([[+ MatrixElement(lrow, mrow, lcolumn, mcolumn).Tr() * B
                        + MatrixElement(lrow, mrow, lcolumn, mcolumn).Vr1_DVR(Xi, Xivec) * V1
                        for lcolumn in range (0, l + 1)
                        for mcolumn in range (-lcolumn, lcolumn + 1)
                        ]
                    for lrow in range (0, l + 1)
                    for mrow in range (-lrow, lrow + 1)
                    ])
    result = result.astype(np.float64)
    return result

def diagonalization(H):
    eig_val,eig_vec = np.linalg.eigh(H)
    return eig_val, eig_vec

def main():
    print("dimension: ", dimension)
    print("l_potential: ", l_pot, "\n")

    print("grid: np.polynomial.legendre.leggauss() \n", np.round(np.polynomial.legendre.leggauss(l + 1)[0], 3))
    print("grid: DVR_spherical_harmonics() \n", np.round(Xi, 3), "\n")
    
    value_SH, vector_SH = diagonalization(Hamiltonian_SH())  #解析解の行列要素によるHamiltonian行列
    print("eigenvalue: analytical solution by spherical harmonics basis: \n", np.round(value_SH, 5))
    value_DVR, vector_DVR = diagonalization(Hamiltonian_DVR())  #DVR行列要素によるHamiltonian行列
    print("eigenvalue: numerical solution by associated Legendre-DVR: \n", np.round(value_DVR, 5))
    return

main()