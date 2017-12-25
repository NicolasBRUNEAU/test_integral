#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def plotMat(test):
    fig = plt.figure()
    ax = Axes3D(fig)
    plt.xlim(0,32)
    plt.ylim(0,32)
    ax.set_zlim(0,32)
    

    for xs in range(len(test[:,0])):
        for ys in range(len(test[:,1])):
            for zs in range(len(test[:,2])):
                if(test[xs, ys, zs] != 0):
                    if(test[xs, ys, zs] == 1):
                        c = 'b'
                    else:
                        c = 'r'
                    ax.scatter(xs,ys,zs,c=c)
    plt.show()

if __name__ == '__main__':
    test = np.zeros((3, 3, 3))
    test[1, 1, 1] = 1
    test[1, 2, 1] = 1
    test[1, 1, 2] = 1
    test[1, 2, 2] = 15
    print(test)
    fig = plt.figure()
    ax = Axes3D(fig)

    for xs in range(len(test[:,0])):
        for ys in range(len(test[:,1])):
            for zs in range(len(test[:,2])):
                '''
                if(test[xs, ys, zs] != 0):
                    if(test[xs, ys, zs] == 1):
                        c = 'b'
                    else:
                        c = 'r'
                '''
                if(test[xs, ys, zs] != 0 and test[xs, ys, zs] != 1):
                    c = 'r'
                    ax.scatter(xs,
                                ys,
                                zs,
                                c=c)
    plt.show()