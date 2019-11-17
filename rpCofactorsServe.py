#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Backend rpCofactors server

"""

import os
import json
import libsbml
from datetime import datetime
from flask import Flask, request, jsonify, send_file, abort
from flask_restful import Resource, Api
import io
import tarfile
import csv
import sys
import glob
import tempfile

sys.path.insert(0, '/home/')
import rpCofactors
import rpCache
import rpSBML

#######################################################
############## REST ###################################
#######################################################

app = Flask(__name__)
api = Api(app)

rpcache = rpCache.rpCache()

def stamp(data, status=1):
    appinfo = {'app': 'rpCofactors', 'version': '1.0', 
               'author': 'Melchior du Lac',
               'organization': 'BRS',
               'time': datetime.now().isoformat(), 
               'status': status}
    out = appinfo.copy()
    out['data'] = data
    return out


## REST App.
#
#
class RestApp(Resource):
    def post(self):
        return jsonify(stamp(None))
    def get(self):
        return jsonify(stamp(None))


## Run a single 
#
#
def runSingleSBML(rpcofactors, member_name, rpsbml_string, pathway_id, compartment_id):
    #open one of the rp SBML files
    rpsbml = rpSBML.rpSBML(member_name, libsbml.readSBMLFromString(rpsbml_string))
    if rpcofactors.addCofactors(rpsbml, compartment_id, pathway_id):
        return libsbml.writeSBMLToString(rpsbml.document).encode('utf-8')
    else:
        return ''


##
#
#
def runCofactors_mem(rpcofactors, inputTar, outputTar, pathway_id='rp_pathway', compartment_id='MNXC3'):
    #loop through all of them and run FBA on them
    with tarfile.open(fileobj=outputTar, mode='w:xz') as tf:
        with tarfile.open(fileobj=inputTar, mode='r:xz') as in_tf:
            for member in in_tf.getmembers():
                if not member.name=='':
                    data = singleCofactors(rpcofactors,
                            member.name,
                            in_tf.extractfile(member).read().decode('utf-8'),
                            pathway_id,
                            compartment_id)
                    if not data=='':
                        fiOut = io.BytesIO(data)
                        info = tarfile.TarInfo(member.name)
                        info.size = len(data)
                        tf.addfile(tarinfo=info, fileobj=fiOut)


## run using HDD 3X less than the above function
#
#
def runCofactors_hdd(rpcofactors, inputTar, outputTar, pathway_id='rp_pathway', compartment_id='MNXC3'):
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        with tempfile.TemporaryDirectory() as tmpInputFolder:
            tar = tarfile.open(fileobj=inputTar, mode='r:xz')
            tar.extractall(path=tmpInputFolder)
            tar.close()
            for sbml_path in glob.glob(tmpInputFolder+'/*'):
                fileName = sbml_path.split('/')[-1].replace('.sbml.xml', '')
                rpsbml = rpSBML.rpSBML(fileName)
                rpsbml.readSBML(sbml_path)
                rpcofactors.addCofactors(rpsbml, compartment_id, pathway_id)
                rpsbml.writeSBML(tmpOutputFolder)
                rpsbml = None
            with tarfile.open(fileobj=outputTar, mode='w:xz') as ot:
                for sbml_path in glob.glob(tmpOutputFolder+'/*'):
                    fileName = str(sbml_path.split('/')[-1].replace('.sbml', ''))
                    fileName += '.sbml.xml'
                    info = tarfile.TarInfo(fileName)
                    info.size = os.path.getsize(sbml_path)
                    ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))


## REST Query
#
# REST interface that generates the Design.
# Avoid returning numpy or pandas object in
# order to keep the client lighter.
class RestQuery(Resource):
    def post(self):
        inputTar = request.files['inputTar']
        params = json.load(request.files['data'])
        #pass the cache parameters to the rpCofactors object
        outputTar = io.BytesIO()
        rpcofactors = rpCofactors.rpCofactors()
        rpcofactors.deprecatedMNXM_mnxm = rpcache.deprecatedMNXM_mnxm
        rpcofactors.deprecatedMNXR_mnxr = rpcache.deprecatedMNXR_mnxr
        rpcofactors.mnxm_strc = rpcache.mnxm_strc
        rpcofactors.full_reactions = rpcache.full_reactions
        rpcofactors.chemXref = rpcache.chemXref
        rpcofactors.rr_reactions = rpcache.rr_reactions
        ######## HDD #######
        runCofactors_hdd(rpcofactors, inputTar, outputTar, params['pathway_id'], params['compartment_id'])
        ######## MEM #######
        #runCofactors_mem(rpcofactors, inputTar, outputTar, params['pathway_id'], params['compartment_id'])
        ###### IMPORTANT ######
        outputTar.seek(0)
        #######################
        return send_file(outputTar, as_attachment=True, attachment_filename='rpCofactors.tar', mimetype='application/x-tar')


api.add_resource(RestApp, '/REST')
api.add_resource(RestQuery, '/REST/Query')


if __name__== "__main__":
    app.run(host="0.0.0.0", port=8996, debug=True, threaded=True)
