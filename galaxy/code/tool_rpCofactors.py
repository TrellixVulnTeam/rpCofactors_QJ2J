#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Galaxy script to query rpCofactors REST service

"""
import argparse
import tarfile
import tempfile
import shutil
import copy
import glob
import os
import logging
import sys
sys.path.insert(0, '/home/')

import rpToolServe

##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Python wrapper to add cofactors to generate rpSBML collection')
    parser.add_argument('-input', type=str)
    parser.add_argument('-output', type=str)
    parser.add_argument('-input_format', type=str)
    parser.add_argument('-rxn_recipes', type=str, default='None')
    parser.add_argument('-rules_rall', type=str, default='None')
    parser.add_argument('-compounds', type=str, default='None')
    parser.add_argument('-pathway_id', type=str, default='rp_pathway')
    parser.add_argument('-compartment_id', type=str, default='MNXC3')
    parser.add_argument('-pubchem_search', type=str, default='False')
    params = parser.parse_args()
    if params.pubchem_search=='True' or params.pubchem_search=='T' or params.pubchem_search=='true' or params.pubchem_search=='t':
        pubchem_search = True
    elif params.pubchem_search=='False' or params.pubchem_search=='F' or params.pubchem_search=='false' or params.pubchem_search=='f':
        pubchem_search = False
    else:
        logging.error('Cannot interpret pubchem_search input: '+str(params.pubchem_search))
        exit(1)
    if params.input_format=='tar':
        if (params.rxn_recipes=='None' or params.rxn_recipes==None) and (params.compounds=='None' or params.compounds==None) and (params.rules_rall=='None' or params.rules_rall==None):
            rpToolServe.main(params.input,
                             params.output,
                             params.pathway_id,
                             params.compartment_id,
                             pubchem_search)
        else:
            if os.path.exists(params.rxn_recipes):
                rpToolServe.main_extrules(params.input,
                                          params.output,
                                          params.rxn_recipes,
                                          params.rules_rall,
                                          params.compounds,
                                          params.pathway_id,
                                          params.compartment_id,
                                          pubchem_search)
    elif params.input_format=='sbml':
        #make the tar.xz 
        with tempfile.TemporaryDirectory() as tmpOutputFolder:
            inputTar = tmpOutputFolder+'/tmp_input.tar'
            outputTar = tmpOutputFolder+'/tmp_output.tar'
            with tarfile.open(inputTar, mode='w:gz') as tf:
                info = tarfile.TarInfo('single.sbml.xml') #need to change the name since galaxy creates .dat files
                info.size = os.path.getsize(params.input)
                tf.addfile(tarinfo=info, fileobj=open(params.input, 'rb'))
            if (params.rxn_recipes=='None' or params.rxn_recipes==None) and (params.compounds=='None' or params.compounds==None) and (params.rules_rall=='None' or params.rules_rall==None):
                rpToolServe.main(inputTar,
                                 outputTar,
                                 params.pathway_id,
                                 params.compartment_id,
                                 pubchem_search)
            else:
                if os.path.exists(params.rxn_recipes):
                    rpToolServe.main_extrules(inputTar,
                                              outputTar,
                                              params.rxn_recipes,
                                              params.rules_rall,
                                              params.compounds,
                                              params.pathway_id,
                                              params.compartment_id,
                                              pubchem_search)
            with tarfile.open(outputTar) as outTar:
                outTar.extractall(tmpOutputFolder)
            out_file = glob.glob(tmpOutputFolder+'/*.sbml.xml')
            if len(out_file)>1:
                logging.warning('There are more than one output file...')
            shutil.copy(out_file[0], params.output)
    else:
        logging.error('Cannot identify the input/output format: '+str(params.input_format))
