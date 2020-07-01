#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Backend rpCofactors server

"""

import os
import json
import libsbml
import copy
import logging
import io
import tarfile
import csv
import sys
import glob
import tempfile
import shutil

sys.path.insert(0, '/home/')
import rpTool as rpCofactors
import rpCache
import rpSBML
import tool_rpUnicity

logging.basicConfig(
    #level=logging.DEBUG,
    #level=logging.WARNING,
    level=logging.ERROR,
    format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S',
)

#logging.getLogger().setLevel(logging.DEBUG)

## Run a single
#
#
def runSingleSBML(rpcofactors, member_name, rpsbml_string, pathway_id, compartment_id, pubchem_search):
    #open one of the rp SBML files
    rpsbml = rpSBML.rpSBML(member_name, libsbml.readSBMLFromString(rpsbml_string))
    if rpcofactors.addCofactors(rpsbml, compartment_id, pathway_id, pubchem_search):
        return libsbml.writeSBMLToString(rpsbml.document).encode('utf-8')
    else:
        return ''


##
#
#
def runCofactors_mem(rpcofactors, inputTar, outputTar, pathway_id='rp_pathway', compartment_id='MNXC3', pubchem_search=False):
    #loop through all of them and run FBA on them
    with tarfile.open(fileobj=outputTar, mode='w:gz') as tf:
        with tarfile.open(fileobj=inputTar, mode='r') as in_tf:
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
def runCofactors_hdd(rpcofactors, inputTar, outputTar, pathway_id='rp_pathway', compartment_id='MNXC3', pubchem_search=False):
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
                logging.debug('============= '+str(fileName)+' ============')
                rpsbml = rpSBML.rpSBML(fileName)
                rpsbml.readSBML(sbml_path)
                rpcofactors.addCofactors(rpsbml, compartment_id, pathway_id)
                rpsbml.writeSBML(tmpOutputFolder)
                rpsbml = None
            if len(glob.glob(tmpOutputFolder+'/*'))==0:
                logging.error('rpCofactors has not produced any results')
                return False
            with tarfile.open(outputTar, mode='w:gz') as ot:
                for sbml_path in glob.glob(tmpOutputFolder+'/*'):
                    fileName = str(sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', ''))
                    fileName += '.sbml.xml'
                    info = tarfile.TarInfo(fileName)
                    info.size = os.path.getsize(sbml_path)
                    ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))
    return True


##
#
#
def main(inputTar,
         outputTar,
         pathway_id='rp_pathway',
         compartment_id='MNXC3',
         pubchem_search=False):
    rpcache = rpCache.rpCache()
    rpcofactors = rpCofactors.rpCofactors()
    rpcofactors.rr_full_reactions = rpcache.getFullReactions()
    rpcofactors.deprecatedCID_cid = rpcache.getDeprecatedCID()
    rpcofactors.deprecatedRID_rid = rpcache.getDeprecatedRID()
    rpcofactors.cid_strc = rpcache.getCIDstrc()
    rpcofactors.inchikey_cid = rpcache.getInchiKeyCID()
    rpcofactors.rr_reactions = rpcache.getRRreactions()
    rpcofactors.cid_xref = rpcache.getCIDxref()
    rpcofactors.xref_comp, rpcofactors.comp_xref = rpcache.getCompXref()
    rpcofactors.chebi_cid = rpcache.getChebiCID()
    rpcofactors.cid_name = rpcache.getCIDname()
    #outputTar_bytes = io.BytesIO()
    ######## HDD #######
    with tempfile.TemporaryDirectory() as tmpResFolder:
        runCofactors_hdd(rpcofactors, inputTar, tmpResFolder+'/tmpRes.tar', pathway_id, compartment_id, pubchem_search)
        #shutil.copyfile(tmpResFolder+'/tmpRes.tar', outputTar)
        #shutil.copyfile(tmpResFolder+'/tmpRes.tar', 'tmpRes.tar')
        ######## MEM #######
        #runCofactors_mem(rpcofactors, inputTar, outputTar, params['pathway_id'], params['compartment_id'])
        ####### rpUnicity ######
        tool_rpUnicity.deduplicate(tmpResFolder+'/tmpRes.tar', outputTar)

##
#
#
def main_extrules(inputTar,
                  outputTar,
                  rxn_recipes,
                  rules_rall_tsv,
                  compounds_tsv,
                  pathway_id='rp_pathway',
                  compartment_id='MNXC3',
                  pubchem_search=False):
    rpcache = rpCache.rpCache()
    rpcofactors = rpCofactors.rpCofactors()
    #### parse the input files and merge with cache ####
    ''' if you want to merge
    rpcache.retroRulesFullReac(rxn_recipes)
    new_full_reactions = copy.deepcopy(rpcache.rr_full_reactions)
    rpcache.rr_full_reactions = {**rpcache.getFullReactions(), **new_full_reactions}
    rpcofactors.rr_full_reactions = rpcache.rr_full_reactions
    #reaction rules
    rpcache.retroReactions(rules_rall_tsv)
    new_rr_reactions = copy.deepcopy(rpcache.rr_reactions)
    rpreader.rr_reactions = {**rpcache.getRRreactions(), **new_rr_reactions}
    '''
    #if you want to overwite
    rpcache.retroRulesFullReac(rxn_recipes)
    rpcofactors.rr_full_reactions = rpcache.rr_full_reactions
    rpcache.retroRulesStrc(compounds_tsv)
    new_cid_strc = copy.deepcopy(rpcache.cid_strc)
    rpcache.cid_strc = {**rpcache.getCIDstrc(), **new_cid_strc}
    rpcache._inchikeyCID()
    rpcofactors.cid_strc = rpcache.cid_strc
    rpcofactors.inchikey_cid = rpcache.inchikey_cid
    #reaction rules
    rpcache.retroReactions(rules_rall_tsv)
    rpcofactors.rr_reactions = rpcache.rr_reactions
    ##
    rpcofactors.deprecatedCID_cid = rpcache.getDeprecatedCID()
    rpcofactors.deprecatedRID_rid = rpcache.getDeprecatedRID()
    rpcofactors.cid_xref = rpcache.getCIDxref()
    rpcofactors.xref_comp, rpcofactors.comp_xref = rpcache.getCompXref()
    rpcofactors.chebi_cid = rpcache.getChebiCID()
    rpcofactors.cid_name = rpcache.getCIDname()
    #outputTar_bytes = io.BytesIO()
    ######## HDD #######
    with tempfile.TemporaryDirectory() as tmpResFolder:
        runCofactors_hdd(rpcofactors, inputTar, tmpResFolder+'/tmpRes.tar', pathway_id, compartment_id, pubchem_search)
        #shutil.copyfile(tmpResFolder+'/tmpRes.tar', outputTar)
        #shutil.copyfile(tmpResFolder+'/tmpRes.tar', 'tmpRes.tar')
        ######## MEM #######
        #runCofactors_mem(rpcofactors, inputTar, outputTar, params['pathway_id'], params['compartment_id'])
        ####### rpUnicity ######
        tool_rpUnicity.deduplicate(tmpResFolder+'/tmpRes.tar', outputTar)

