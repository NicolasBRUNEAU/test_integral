#!/usr/bin/python3
# coding: utf-8


import logging
import numpy as np
import copy
import math
from fastFourierTransform import *
from color import *

def dist(coord1, coord2):
    """
        Calcul distances between two points
    """
    # alternatives
    #     numpy.sqrt(numpy.sum((coord[numAtom]-fftGrid[numNode])**2))
    return np.linalg.norm(coord1 - coord2)


"""
Class molecule
"""

class Molecule(object):
    """
    Molecule def class

    It's attributes are :
        
    """
    def __init__(self, npdb):
        self.npdb = npdb
        self.atoms = []
        self.coord = []
        self.fftGrid = []
        self.dim = 0
        self.cdm = [0,0,0] #centre de masse, (x y z)
        self.read_pdb(npdb)
        self.calc_cdm()
        self.dim_max()
        

    def read_pdb(self, npdb):
        logging.debug("Reading pdb")
        self.npdb = npdb #split path et .pdb, juste garder nom
        listatoms = []
        coordinates = []
        with open(npdb, 'r') as fin:
            for line in fin:
                if(line.startswith('ATOM') and line[11:16].strip()[0] != 'H'):
                    listatoms.append((line[:6], line[11:30], line[54:]))
                    x = line[30:38].strip()
                    y = line[38:46].strip()
                    z = line[46:54].strip()
                    coordinates.append([float(x), float(y), float(z)])
        self.atoms = tuple(listatoms)
        self.coord = np.array(coordinates)
        

    def write_pdb(self, npdb, i = 0, chain = 64):
        logging.debug("Writting pdb")
        chaintmp = ''
        chaintmpPrec = ''
        with open(npdb, 'a') as fout: #ajouter le path
            for j, atom in enumerate(self.atoms):
                chaintmp = atom[1][10:11]
                if(chaintmp != chaintmpPrec):
                    chain+=1
                    chaintmpPrec = chaintmp
                    
                fout.write(atom[0] + '{:5d}'.format(i+1) + atom[1][0:10] + chr(chain) + atom[1][11:18] + '{:8.3f}{:8.3f}{:8.3f}'.format(self.coord[j][0], self.coord[j][1], self.coord[j][2]) + atom[2])
                i += 1
            fout.write("\nTER\n")
            return i, chain

    def init_fftgrid(self, internValue, size, r, t):
        """
            Permet de créer la grille pour le calcul de FFT.
        """
        
        self.fftGrid = np.zeros((size, size, size), dtype=float)

        # valeurs de r, t, rho et delta en fonction de l'article "Rigid Body Protein Docking by Fast Fourier Transform"

        # on va donc parcourir les coordonnées de la structure (ligand ou récepteur) et pour chaque atome on va regarder

        boundary = r + t
        for atom in self.coord:
            lim = self.limits(atom, boundary)

            for i in range(lim[0][0], lim[0][1]):
                for j in range(lim[1][0], lim[1][1]):
                    for k in range(lim[2][0], lim[2][1]):

                        d = dist(atom, [i, j, k])
                        if(d < r):
                            self.fftGrid[i, j, k] = internValue
                        elif(d < boundary):
                            if(self.fftGrid[i, j, k] != internValue):
                                self.fftGrid[i, j, k] = 1

    def rotate(self, ligand, alpha, beta, gamma):
        """
        initial author : enseignant meet-u
        rewrite : groupe 2 : implementation in class
        purpose: rotation of the atom with coords (x, y, z) according to the angles alpha, beta, gamma
        input:   alpha, beta, gamma, the angles of the rotation
        change : the final coordinates in the initial referential for each atom
        """
    
        for i, atom in enumerate(ligand.coord):
            # centering according to the center of rotation
            x1i = atom[0] - self.cdm[0]
            y1i = atom[1] - self.cdm[1]
            z1i = atom[2] - self.cdm[2]

            # computing cos and sin of each angle
            c_a = math.cos(alpha)
            s_a = math.sin(alpha)

            c_b = math.cos(beta)
            s_b = math.sin(beta)

            c_g = math.cos(gamma)
            s_g = math.sin(gamma)

            # applying rotation
            x3i = (c_a*c_b*c_g-s_a*s_g)*x1i + (-c_a*c_b*s_g-s_a*c_g)*y1i + c_a*s_b*z1i 
            y3i = (c_a*s_g+s_a*c_b*c_g)*x1i + (-s_a*c_b*s_g+c_a*c_g)*y1i + s_a*s_b*z1i
            z3i = -s_b*c_g*x1i + s_b*s_g*y1i + c_b*z1i 

            # back to the input referential
            self.coord[i, 0] = x3i + self.cdm[0]
            self.coord[i, 1] = y3i + self.cdm[1]
            self.coord[i, 2] = z3i + self.cdm[2]

    def translate(self, vect):
        """
            Permet de translater la structure
        """
        self.coord = self.coord + vect



    def calc_cdm(self):
        """
            Donne les coordonnées du centre de la structure (ligand ou récepteur).
        """
        for i in range(self.coord.shape[1]):
            self.cdm[i] = sum(self.coord[:,i]) / self.coord.shape[0]

    def dim_max(self):
        dimProt = [0, 0, 0]
        for i in range(3):
            dimProt[i] = max(self.coord[:,i]) - min(self.coord[:,i])
        self.dim = int(max(dimProt))

    def limits(self, atom, boundary):
        """
        evite le débordement grille
        """
        lim = []
        for i in range(3):
            inf = max(0, int(math.floor(atom[i] - boundary)))
            sup = min(self.fftGrid.shape[0], int(math.ceil(atom[i] + boundary)))
            lim.append([inf, sup])

        return(lim)




