#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Backend rpCofactors server

"""

import os
import shutil
import json
import libsbml
from datetime import datetime
from flask import Flask, request, jsonify, send_file, abort
from flask_restful import Resource, Api
import io
import tarfile
import csv
import sys
<<<<<<< HEAD

sys.path.insert(0, '/home/')
import rpCofactors
import rpSBML


########################### Processify each #################
######## code taken from 


import inspect
import traceback

from functools import wraps
from multiprocessing import Process, Queue


class Sentinel:
    pass


def processify(func):
    '''Decorator to run a function as a process.
    Be sure that every argument and the return value
    is *pickable*.
    The created process is joined, so the code does not
    run in parallel.
    '''

    def process_generator_func(q, *args, **kwargs):
        result = None
        error = None
        it = iter(func())
        while error is None and result != Sentinel:
            try:
                result = next(it)
                error = None
            except StopIteration:
                result = Sentinel
                error = None
            except Exception:
                ex_type, ex_value, tb = sys.exc_info()
                error = ex_type, ex_value, ''.join(traceback.format_tb(tb))
                result = None
            q.put((result, error))

    def process_func(q, *args, **kwargs):
        try:
            result = func(*args, **kwargs)
        except Exception:
            ex_type, ex_value, tb = sys.exc_info()
            error = ex_type, ex_value, ''.join(traceback.format_tb(tb))
            result = None
        else:
            error = None

        q.put((result, error))

    def wrap_func(*args, **kwargs):
        # register original function with different name
        # in sys.modules so it is pickable
        process_func.__name__ = func.__name__ + 'processify_func'
        setattr(sys.modules[__name__], process_func.__name__, process_func)

        q = Queue()
        p = Process(target=process_func, args=[q] + list(args), kwargs=kwargs)
        p.start()
        result, error = q.get()
        p.join()

        if error:
            ex_type, ex_value, tb_str = error
            message = '%s (in subprocess)\n%s' % (str(ex_value), tb_str)
            raise ex_type(message)

        return result

    def wrap_generator_func(*args, **kwargs):
        # register original function with different name
        # in sys.modules so it is pickable
        process_generator_func.__name__ = func.__name__ + 'processify_generator_func'
        setattr(sys.modules[__name__], process_generator_func.__name__, process_generator_func)

        q = Queue()
        p = Process(target=process_generator_func, args=[q] + list(args), kwargs=kwargs)
        p.start()

        result = None
        error = None
        while error is None:
            result, error = q.get()
            if result == Sentinel:
                break
            yield result
        p.join()

        if error:
            ex_type, ex_value, tb_str = error
            message = '%s (in subprocess)\n%s' % (str(ex_value), tb_str)
            raise ex_type(message)

    @wraps(func)
    def wrapper(*args, **kwargs):
        if inspect.isgeneratorfunction(func):
            return wrap_generator_func(*args, **kwargs)
        else:
            return wrap_func(*args, **kwargs)
    return wrapper


###############################################
###############################################
###############################################

############## run all using processify ####

@processify
def runSingleSBML(rpcofactors, member_name, rpsbml_string, path_id, compartment_id):
    #open one of the rp SBML files
    rpsbml = rpSBML.rpSBML(member_name, libsbml.readSBMLFromString(rpsbml_string))
    #rpcofactors = rpRanker.rpCofactors()
    if rpcofactors.addCofactors(rpsbml, compartment_id, path_id):
        return libsbml.writeSBMLToString(rpsbml.document).encode('utf-8')
    else:
        return ''


def runAllSBML(inputTar, outputTar, path_id, compartment_id):
    #loop through all of them and run FBA on them
    with tarfile.open(fileobj=outputTar, mode='w:xz') as tf:
        with tarfile.open(fileobj=inputTar, mode='r:xz') as in_tf:
            for member in in_tf.getmembers():
                if not member.name=='':
                    data = runSingleSBML(rpcofactors,
                            member.name,
                            in_tf.extractfile(member).read().decode("utf-8"),
                            path_id,
                            compartment_id)
                    if not data=='':
                        fiOut = io.BytesIO(data)
                        info = tarfile.TarInfo(member.name)
                        info.size = len(data)
                        tf.addfile(tarinfo=info, fileobj=fiOut)


=======
import glob
import random
import string

sys.path.insert(0, '/home/')
import rpCofactors
import rpCache
import rpSBML

>>>>>>> master_rest
#######################################################
############## REST ###################################
#######################################################

app = Flask(__name__)
api = Api(app)
<<<<<<< HEAD
#dataFolder = os.path.join( os.path.dirname(__file__),  'data' )


#TODO: test that it works well
#declare the rpReader globally to avoid reading the pickle at every instance
rpcofactors = rpCofactors.rpCofactors()


=======

rpcache = rpCache.rpCache()


## Stamp of rpCofactors
#
#
>>>>>>> master_rest
def stamp(data, status=1):
    appinfo = {'app': 'rpCofactors', 'version': '1.0', 
               'author': 'Melchior du Lac',
               'organization': 'BRS',
               'time': datetime.now().isoformat(), 
               'status': status}
    out = appinfo.copy()
    out['data'] = data
    return out


<<<<<<< HEAD
class RestApp(Resource):
    """ REST App."""
=======
## REST App.
#
#
class RestApp(Resource):
>>>>>>> master_rest
    def post(self):
        return jsonify(stamp(None))
    def get(self):
        return jsonify(stamp(None))


<<<<<<< HEAD
class RestQuery(Resource):
    """ REST interface that generates the Design.
        Avoid returning numpy or pandas object in
        order to keep the client lighter.
    """
    def post(self):
        inputTar = request.files['inputTar']
        params = json.load(request.files['data'])
        #pass the files to the rpReader
        outputTar = io.BytesIO()
        runAllSBML(inputTar, outputTar, params['path_id'], params['compartment_id'])
=======
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
    if not os.path.exists(os.getcwd()+'/tmp'):
        os.mkdir(os.getcwd()+'/tmp')
    tmpInputFolder = os.getcwd()+'/tmp/'+''.join(random.choice(string.ascii_lowercase) for i in range(15))
    tmpOutputFolder = os.getcwd()+'/tmp/'+''.join(random.choice(string.ascii_lowercase) for i in range(15))
    os.mkdir(tmpInputFolder)
    os.mkdir(tmpOutputFolder)
    tar = tarfile.open(fileobj=inputTar, mode='r:xz')
    tar.extractall(path=tmpInputFolder)
    tar.close()
    for sbml_path in glob.glob(tmpInputFolder+'/*'):
        fileName = sbml_path.split('/')[-1].replace('.sbml', '')
        rpsbml = rpSBML.rpSBML(fileName)
        rpsbml.readSBML(sbml_path)
        rpcofactors.addCofactors(rpsbml, compartment_id, pathway_id)
        rpsbml.writeSBML(tmpOutputFolder)
        rpsbml = None
    with tarfile.open(fileobj=outputTar, mode='w:xz') as ot:
        for sbml_path in glob.glob(tmpOutputFolder+'/*'):
            fileName = str(sbml_path.split('/')[-1].replace('.sbml', ''))
            info = tarfile.TarInfo(fileName)
            info.size = os.path.getsize(sbml_path)
            ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))
    shutil.rmtree(tmpInputFolder)
    shutil.rmtree(tmpOutputFolder)


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
>>>>>>> master_rest
        ###### IMPORTANT ######
        outputTar.seek(0)
        #######################
        return send_file(outputTar, as_attachment=True, attachment_filename='rpCofactors.tar', mimetype='application/x-tar')


api.add_resource(RestApp, '/REST')
api.add_resource(RestQuery, '/REST/Query')


if __name__== "__main__":
<<<<<<< HEAD
    #debug = os.getenv('USER') == 'mdulac'
    app.run(host="0.0.0.0", port=8996, debug=True, threaded=True)

=======
    app.run(host="0.0.0.0", port=8996, debug=True, threaded=True)
>>>>>>> master_rest
