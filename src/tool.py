#!/usr/bin/python
# coding: utf-8

import os
from molecule import *
import copy as cp
import sys
import argparse


def printMat1(mat):
    print('==========================================')
    for ligne in mat:
            for ele in ligne:
                print(ele)
            print("\n")

def search_max_corr(grid):
    tmp_coord = [0, 0, 0]
    tmp_value = -1e400
    gSize = grid.shape[0]
    #print(gSize)
    for x in range(gSize):
        for y in range(gSize):
            for z in range(gSize):
                if grid[x, y, z] > tmp_value:
                    tmp_value = grid[x, y, z]
                    #print (tmp_value)
                    tmp_coord = [(x - gSize / 2), (y - gSize / 2), (z - gSize / 2)]
                    #print (tmp_coord)
    return np.array(tmp_coord)

def calculRotate_CreateGrid_FFT(alpha, delta, dAngle, gridSize, r, t, conjFFTrecepteur, ligand, recepteur):
    """
    Fonction lancée par chaque appel de nouveau Process.
    La fonction parcours deux différents espaces conformationnels en Y (Beta) et Z (Gamma).
    Chaque nouvelle valeur de Beta et Gamma fait :
        - Une rotation du ligand par rapport a la position initiale en alpha (envoye en parametre), Beta et Gamma
        - Remplissage de la grille de Fourrier du ligand ayant subit la rotation
        - Calcul d'une FFT entre la nouvelle grille du ligand et le conjugue de la grille du recepteur (en parametre)
        - Recherche du vecteur de translation de la FFT
        - Ecriture du nouveau PDB avec le recepteur + ligand dont les coordonnées ont été vectorisé par le veceteur de translation

    Les parametres :
    * alpha : valeur de alpha issue de la parallélisation du processus
    * delta : valeur de remplissage de la grille de l'interieur du ligand
    * dAngle : pas de variation de la valeur de rotation du ligand (en degree)
    * gridSize : taille de la grille
    * conjFFTrecepteur : conjugue de la matrice de Fourrier du recepteur
    * ligand : coordonnées des atomes du ligand
    * recepteur : coordonnées des atomes du recepteur
    """
    beta = 0
    while beta <= 360:
        gamma = 0
        while gamma <= 180:
            
            ligandTmp = cp.deepcopy(ligand)

            # rotation du ligand autour du recepteur
            ligandTmp.rotate(ligand, alpha, beta, gamma)

            # recalcule de la nouvelle grille du ligand rotationne
            ligandTmp.init_fftgrid(delta, gridSize, r, t)

            # calcul nouvel fft
            matCor = fast_fourier_transform(conjFFTrecepteur, ligandTmp.fftGrid, gridSize, plot=False)
            
            vectTranslation = search_max_corr(matCor)
            #print(vectTranslation)

            ligandTmp.translate(vectTranslation)

            # test ecriture fichier
            fileName = str(alpha) + "_" + str(beta) + "_" + str(gamma) + ".pdb"
            write_pdb_complex(ligand, recepteur, fileName)

            #print('parentid: ' + str(os.getppid()), '  childid: ' + str(os.getpid()), 'alpha: ' + str(alpha), '  beta: ' + str(beta), '  gamma: ' + str(gamma))
            gamma += dAngle
            
        
        beta += dAngle
        

def write_pdb_complex(ligand, recepteur, name_pdb):
    """
    Fonction d'écriture du nouveau pdb généré
    Les numeros d'atomes vont à la suite, pas de remise à zéro des num atomes pour le ligand
    """
    # creer/ecraser fichier
    path = 'results/PDB/' + name_pdb
    if(os.path.isfile(path)):
         os.remove(path)
    i, chain = recepteur.write_pdb(path)
    ligand.write_pdb(path, i, chain)

def gestion_parametres():

    #Dictionnaire de paramètres (Valeurs par défaut)
    parametres = {"Receptor": None, "Ligand": None, "size": 128, "t": 1.1, "r": 1.1, "rho": -5, "delta": 1, "pas": 12}

    '''
    Gestion des paramètres
    8 paramètres :
    -h : aide
    -d : Seuil de druggabilité des poches de type float entre 0 et 1. Seules les poches ayant un score 
    de druggabilité supérieur seront considérées.
    -s : Indique le chemin vers le repertoire de sortie souhaité.
    -c : Indique le chemin vers le fichier pdb de la poche cible. Seules les poches
    ayant un recouvrement avec la cible seront considérées, chemin/fichier.pdb, par defaut = -1 -> pas de cible
    -o : Seuil de recouvrement des poches de type float entre 0 et 100 %. Seules les recouvrements supérieurs au seuil 
    seront considérées.
    -e : Méthode de détection des poches (1: prox,2: fpocket,3: both), seule la méthode 2 est fonctionelle.
    -p : Indique le chemin vers le fichier trajectoire (PDB ou XTC)
    -t : Indique le chemin vers le fichier topologie (TPR) si la trajectoire est au format XTC.
    '''

    if len(sys.argv)>17:
        sys.exit("erreur : trop de parametres (voir aide -h)") 

    parser = argparse.ArgumentParser()
    parser.add_argument("--rec", help="Chemin du fichier pdb du récepteur.")
    parser.add_argument("--lig", help="Chemin du fichier pdb du ligand")
    parser.add_argument("-s", type=int, help="Taille de la grille (en puissance de 2). Defaut: 128")
    parser.add_argument("-r", type=float, help="Rayon définissant les noeuds à l'intérieur de la protéine. Defaut: 1.1")
    parser.add_argument("-t", type=float, help="Rayon définissant les noeuds à la surface de la protéine. Defaut: 1.1")
    parser.add_argument("-p", type=float, help="Rho, valeur associée au noeuds à l'intérieur du récepteur (négatif dans la méthode). Defaut: -5")
    parser.add_argument("-d", type=float, help="Delta, valeur associée au noeuds à l'intérieur du ligand (entre 0 et 1 dans la méthode). Defaut: 1")
    parser.add_argument("--pas", type=int, help="Pas de rotation en degré. Defaut: 12")
   
    #Ajouter les autres paramètres de pockdrug ?

    args = parser.parse_args()

    

    '''Verifier le paramètre -s '''
    if args.s!=None:
        if args.s in [32,64,128,256,512]:
            parametres["size"] = int(args.s)
        else:
            sys.exit("erreur: valeur de taille de grille")

    '''Verifier le paramètre -R '''
    if args.rec!=None:
        if os.path.exists(args.rec):
            parametres["Receptor"] = os.path.abspath(args.rec)
        else:
            sys.exit("erreur: chemin récepteur inexistant")
    else:
        sys.exit("erreur: chemin récepteur obligatoire")

    '''Verifier le paramètre -L '''
    if args.lig!=None:
        if os.path.exists(args.lig):
            parametres["Ligand"] = os.path.abspath(args.lig)
        else:
            sys.exit("erreur: chemin ligand inexistant")
    else:
        sys.exit("erreur: chemin ligand obligatoire")

    if args.t!=None:
        if args.t > 0 and args.t < 3:
            parametres["t"] = float(args.t)
        else:
            sys.exit("erreur: paramètre t doit être 0 < t < 3")

    if args.r!=None:
        if args.r > 0 and args.r < 3:
            parametres["r"] = float(args.r)
        else:
            sys.exit("erreur: paramètre r doit être 0 < r < 3")

    if args.p!=None:
        if args.p < 0 :
            parametres["rho"] = float(args.p)
        else:
            sys.exit("erreur: paramètre p doit être < 0")

    if args.d!=None:
        if args.d > 0 and args.d <= 1 :
            parametres["delta"] = float(args.d)
        else:
            sys.exit("erreur: paramètre d doit être 0 < d < 1")

    if args.pas!=None:
        if args.pas > 0 and args.pas < 360:
            parametres["pas"] = int(args.pas)
        else:
            sys.exit("erreur: paramètre pas doit être 0 < pas < 360")




    print "Recepteur :", parametres["Receptor"]
    print "Ligand :", parametres["Ligand"]
    print "Taille grille :", parametres["size"]
    print "t :", parametres["t"]
    print "r :", parametres["r"]
    print "rho :", parametres["rho"]
    print "delta :", parametres["delta"]
    print "Pas de rotation :", parametres["pas"]

    return parametres


