#!/usr/bin/python
# coding: utf-8

import multiprocessing as mp
import os
import time
from fastFourierTransform import *
from molecule import *
from tool import *
from color import plotMat
import copy as cp

########################################
#MAIN

def main():

    parametres = gestion_parametres()

    '''
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Modifier les paramètres pas de rotation, t, r
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    '''

    tmps1 = time.time()
    jobs = []
    alpha = 0

    # import des fichiers pdb
    
    # structure initale garde pour appliquer rotation sur nouveaux structure
    ligand = Molecule(parametres["Ligand"])

    recepteur = Molecule(parametres["Receptor"])


    ligand.translate(np.array([parametres["size"],parametres["size"],parametres["size"]]) / 2 - ligand.cdm)
    recepteur.translate(np.array([parametres["size"],parametres["size"],parametres["size"]]) / 2 - recepteur.cdm)
    ligand.calc_cdm()
    recepteur.calc_cdm()

    recepteur.init_fftgrid(parametres["rho"], parametres["size"], parametres["r"], parametres["t"])
    ligand.init_fftgrid(parametres["delta"], parametres["size"], parametres["r"], parametres["t"])
    #plotMat(ligand.fftGrid)
    #plotMat(recepteur.fftGrid)
    #exit()

    # calcul du taille de la grille en fonction de la taille du ligand et du récepteur
    # création de la grille pour le ligand et le récepteur
    # Une grille pour un recepteur
    # n grille du ligand autour du recepteur


    # exploration du domaine conformationnel autour du recepteur
    conjFFTrecepteur = conj_fourier(recepteur.fftGrid)

    #calcul du nombre de confs attendu et du temps d'exécution
    nbConfsAttendus = ((360/parametres["pas"])+1)**2*((180/parametres["pas"])+1)
    tempsAttendu = (nbConfsAttendus * 0.88)
    print("Nb conformations attendus: " + str(nbConfsAttendus))
    print("Temps d'exécution attendu: " + str(tempsAttendu / 60) + " min soit " + str(tempsAttendu) + " s")
    
    # Création d'un tableau de Process pour calcul partagé
    
    while alpha <= 360:
        recv_end, send_end = mp.Pipe(False)
        thread = mp.Process(target=calculRotate_CreateGrid_FFT,
                            args=(alpha,
                                parametres["delta"],
                                parametres["pas"],
                                parametres["size"],
                                parametres["r"],
                                parametres["t"],
                                conjFFTrecepteur,
                                ligand,
                                recepteur,
                                )
                            )
        # ajoute le thread à la liste de travail
        jobs.append(thread)
        # lancer le thread
        thread.start()
        # incrementation de alpha
        alpha += parametres["pas"]



    # attendre fin processus fils pour continuer execution du main
    for th in jobs:
        th.join()



    tmps2 = time.time()-tmps1
    print("Temps d'execution = %f" %tmps2)
  


if __name__ == "__main__":
    main()