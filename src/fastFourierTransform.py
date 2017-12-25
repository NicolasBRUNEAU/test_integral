#!/usr/bin/python
# coding: utf-8

import numpy as np

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

def printMat(mat):
    print('==========================================')
    for ligne in mat:
        for ele in ligne:
            print(ele)
        print("\n")

def plot_corr(mfc, n):
    #Représentation graphique de la grille.
    i=0
    while(i<n):

        #Reprsentation plot
        x = np.arange(-n/2,n/2)
        y = np.arange(-n/2,n/2)
        x, y = np.meshgrid(x, y)
        # print(x)
        # print(y)

        #pour z = 0
        cor = mfc[:,:,i]
        cor[cor<0]=0
        '''
        fig = plt.figure()
        ax = Axes3D(fig)
        plt.xlim(-n/2,n/2)
        plt.ylim(-n/2,n/2)
        plt.xlabel("X")
        plt.ylabel("Y")
        ax.set_zlim(0,10000)

        #ax.scatter(x, y, z)
        surf = ax.plot_surface(x, y, cor, cmap=cm.coolwarm,linewidth=0, antialiased=False)
        #plt.show()
        fig.savefig("results/correlation_"+str(i-n/2)+".png")
        plt.close(fig)
        
        '''
        fig = plt.figure()
        plt.imshow(cor)
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.colorbar()
        fig.savefig("results/correlation_"+str(i-n/2)+".png")
        plt.close(fig)

        i+=1

def conj_fourier(grid):
    """
    Return le conjugue d'une matrice de fourier
    """
    #return np.conj(np.fft.rfftn(grid, s=(128,128,128), axes=(0,1,2)))
    return np.conj(np.fft.fftn(grid))


def fast_fourier_transform(mfa, fb, n, plot=False):
    """
    Calcul la matrice de corrélation à partir des 2 matrices fa et fb dans taille N*N*N.
    """
    #mfa -> conjugé complexe de la FFT de la grille proteine recepteur
    # mfa = np.conj(np.fft.fftn(fa))
    #mfb -> FFT de la grille proteine ligand
    #mfb = np.fft.rfftn(fb, s=(n,n,n), axes=(0,1,2))
    mfb = np.fft.fftn(fb)

    #print("FFT LIGAND")
    #printMat(mfb)

    tmp = np.matmul(mfa,mfb)
    #print("MATMUL")
    #printMat(tmp)


    #mfc -> inverse de la matrice de corrélation
    mfc = np.fft.ifftn(tmp).real
    #print("CORRELATION")
    #printMat(mfc)

    mfc = np.roll(mfc, n/2, axis=0)
        
    if(plot):
        plot_corr(mfc, n)

    return mfc

if __name__ == '__main__':

    m1 = np.array([
        [[0,0,0], [0,1,0], [0,0,0]],
        [[0,1,0], [1,-0.5,1], [0,1,0]],
        [[0,0,0], [0,1,0], [0,0,0]],
        ])

    mfa1 = np.conj(np.fft.fftn(m1))

    m2 = np.array([
        [[0,0,0], [0,1,0], [0,0,0]],
        [[0,1,0], [1,0.5,1], [0,1,0]],
        [[0,0,0], [0,1,0], [0,0,0]],
        ])

    fast_fourier_transform(mfa1, m2, 3, plot=True)


m1 = np.array([
    [[0,0,0], [0,1,0], [0,0,0]],
    [[0,1,0], [1,-0.5,1], [0,1,0]],
    [[0,0,0], [0,3,0], [0,0,0]],
    ])