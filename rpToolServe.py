#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Backend rpCofactors server

"""

import os
import json
import libsbml
import io
import tarfile
import csv
import sys
import glob
import tempfile
import shutil

sys.path.insert(0, '/home/')
import rpTool as rpCofactors
import rpToolCache
import rpSBML


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
            tar = tarfile.open(inputTar, mode='r')
            tar.extractall(path=tmpInputFolder)
            tar.close()
            if len(glob.glob(tmpInputFolder+'/*'))==0:
                logging.error('Input file is empty')
                return False
            for sbml_path in glob.glob(tmpInputFolder+'/*'):
                fileName = sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '')
                rpsbml = rpSBML.rpSBML(fileName)
                rpsbml.readSBML(sbml_path)
                rpcofactors.addCofactors(rpsbml, compartment_id, pathway_id)
                rpsbml.writeSBML(tmpOutputFolder)
                rpsbml = None
            if len(glob.glob(tmpOutputFolder+'/*'))==0:
                logging.error('rpCofactors has not produced any results')
                return False
            with tarfile.open(outputTar, mode='w:xz') as ot:
                for sbml_path in glob.glob(tmpOutputFolder+'/*'):
                    fileName = str(sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', ''))
                    fileName += '.rpsbml.xml'
                    info = tarfile.TarInfo(fileName)
                    info.size = os.path.getsize(sbml_path)
                    ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))
    return True


##
#
#
def main(inputTar,
         outputTar,
         pathway_id,
         compartment_id):
    rpcache = rpToolCache.rpToolCache()
    rpcofactors = rpCofactors.rpCofactors()
    rpcofactors.deprecatedMNXM_mnxm = rpcache.deprecatedMNXM_mnxm
    rpcofactors.deprecatedMNXR_mnxr = rpcache.deprecatedMNXR_mnxr
    rpcofactors.mnxm_strc = rpcache.mnxm_strc
    rpcofactors.full_reactions = rpcache.full_reactions
    rpcofactors.chemXref = rpcache.chemXref
    rpcofactors.rr_reactions = rpcache.rr_reactions
    #pass the files to the rpReader
    #outputTar_bytes = io.BytesIO()
    ######## HDD #######
    runCofactors_hdd(rpcofactors, inputTar, outputTar, pathway_id, compartment_id)
    ######## MEM #######
    #runCofactors_mem(rpcofactors, inputTar, outputTar, params['pathway_id'], params['compartment_id'])
    ########## IMPORTANT #####
    #outputTar_bytes.seek(0)
    ##########################
    #with open(outputTar, 'wb') as f:
    #    shutil.copyfileobj(outputTar_bytes, f, length=131072)
