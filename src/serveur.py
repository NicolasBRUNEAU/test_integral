#!/usr/bin/env python3.5
# -*- coding: utf-8 -*-


import sys
import os
from flask import Flask, render_template, redirect, request, url_for, flash
from werkzeug.exceptions import default_exceptions
from werkzeug import secure_filename
import datetime


# UPLOAD_FOLDER = '/home/sdv/m2bi/nbruneau/Documents/M2BI/meet_u/2017-2018_Equipe2/data'
UPLOAD_FOLDER = '/home/nicolas/M2BI/Meetu/test_integral/results'
ALLOWED_EXTENSIONS = set(['txt', 'pdb'])
directory = os.getcwd()
server = Flask(__name__,
                    static_folder=directory+'/src/static',
                    template_folder=directory+'/src/templates')
server.secret_key = 'clE_Ma9!k'
server.config['TRAP_HTTP_EXCEPTIONS'] = True
server.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER


@server.route('/')
def accueil():
    return(render_template('index.html',
           directory=directory),
          200)

@server.route('/run', methods=['GET', 'POST'])
def run():
    if request.method == 'POST':
        
        # Traitement du fichier recpeteur
        f_receptor = request.files['receptor']
        filename_r = secure_filename(f_receptor.filename)
        date = datetime.datetime.now()
        # test fichier recepteur
        if (filename_r != ''):
            # traitement fichier ligand
            f_ligand = request.files['ligand']
            filename_l = secure_filename(f_ligand.filename)

            # test fichier ligand
            if (filename_l != ''):
                # création dossier results pour nouveau conformation
                id_dir = str(date.year) + str(date.month) + str(date.day) + str(date.hour) + str(date.minute) + str(date.second) + filename_r + filename_l
                os.mkdir(os.path.join(server.config['UPLOAD_FOLDER'],
                                        id_dir))

                # import fichier pdb recepteur
                f_receptor.save(os.path.join(server.config['UPLOAD_FOLDER'],
                                            id_dir,
                                            filename_r))
                # import fichier pdb ligand
                f_ligand.save(os.path.join(server.config['UPLOAD_FOLDER'],
                                            id_dir,
                                            filename_l))

                # reussite des test et run avec les deux pdb upload
                return(redirect(url_for('results', index=id_dir)))
                return(render_template('results.html/',
                                        resultat=id_dir),
                        200)
            else:
                # erreur fichier ligand
                flash(u'Erreur avec le fichier PDB du ligand', 'erreur')
                return(render_template('upload.html'),
                    200)
        else:
            # erreur fichier recepteur
            flash(u'Erreur avec le fichier PDB du récepteur', 'erreur')
            return(render_template('upload.html'),
                    200)
    else:
        # cas par défaut de la page
        return(render_template('upload.html'),
                200)

# 2017122516835tes1.pdbtes2.pdb
@server.route('/results/<index>', methods=['GET', 'POST'])
def results(index):
    return(render_template('results.html',
                            file=index),
            200)

@server.route('/example')
def example():
    return(render_template('example.html',
           page_name='Example',
           div_container='example.html',
           directory=directory),
          200)

@server.route('/contact')
def contact():
    return(render_template('contact.html',
           page_name='Contact',
           div_container='contact.html',
           directory=directory),
          200)

@server.route('/download')
def download():
    return(render_template('download.html',
           page_name='Download',
           div_container='download.html',
           directory=directory),
          200)

# 200 : succes de la requete
@server.errorhandler(400)
@server.errorhandler(401)
@server.errorhandler(403)
@server.errorhandler(404)
@server.errorhandler(405)
@server.errorhandler(500)
@server.errorhandler(501)
@server.errorhandler(502)
@server.errorhandler(503)
@server.errorhandler(504)
@server.errorhandler(505)
def page_not_found(error):
    dico_erreur={
        400:'Échec de l\'analyse HTTP.',
        401:'Le pseudo ou le mot de passe n\'est pas correct !',
        403:'Requête interdite !',
        404:'La page n\'existe pas ou plus !',
        405:'Méthode non autorisée.',
        500:'Erreur interne au serveur ou serveur saturé.',
        501:'Le serveur ne supporte pas le service demandé.',
        502:'Mauvaise passerelle.',
        503:' Service indisponible.',
        504:'Trop de temps à la réponse.',
        505:'Version HTTP non supportée.'
    }
    return(render_template('error.html',
                                erreur=error.code,
                                specification=dico_erreur[error.code]),
            error.code)


if __name__ == '__main__':
    server.run(debug=True, port=5000)
