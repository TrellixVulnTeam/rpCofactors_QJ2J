#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Galaxy script to query rpCofactors REST service

"""
import argparse
import tempfile
import rpToolServe
import tarfile

##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Python wrapper to add cofactors to generate rpSBML collection')
    parser.add_argument('-inputTar', type=str)
    parser.add_argument('-pathway_id', type=str)
    parser.add_argument('-compartment_id', type=str)
    parser.add_argument('-outputTar', type=str)
    parser.add_argument('-sbml', type=str)
    params = parser.parse_args()
    if params.sbml=='None' or params.sbml==None or params.sbml=='':
        if params.inputTar=='None' or params.inputTar==None or params.inputTar=='':
            logging.error('Cannot have no SBML and no TAR input')
            exit(0)
        rpToolServe.main(params.inputTar,
                         params.outputTar,
                         params.pathway_id,
                         params.compartment_id)
    else:
        #make the tar.xz 
        with tempfile.TemporaryDirectory() as tmpOutputFolder:
            inputTar = tmpOutputFolder+'/tmp_input.tar.xz'
            with tarfile.open(inputTar, mode='w:xz') as tf:
                tf.add(params.sbml)
            rpToolServe.main(inputTar,
                             params.outputTar,
                             params.pathway_id,
                             params.compartment_id)
