#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Extract the sink from an SBML into RP2 friendly format

"""
import argparse
import tempfile
import os
import logging
import shutil
import docker


##
#
#
def main(inputfile,
         input_format,
         rxn_recipes,
         rules_rall,
         compounds,
         output,
         pathway_id='rp_pathway',
         compartment_id='MNXC3',
         pubchem_search='False'):
    docker_client = docker.from_env()
    image_str = 'brsynth/rpcofactors-standalone:extrules'
    try:
        image = docker_client.images.get(image_str)
    except docker.errors.ImageNotFound:
        logging.warning('Could not find the image, trying to pull it')
        try:
            docker_client.images.pull(image_str)
            image = docker_client.images.get(image_str)
        except docker.errors.ImageNotFound:
            logging.error('Cannot pull image: '+str(image_str))
            exit(1)
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        shutil.copy(inputfile, tmpOutputFolder+'/input.dat')
        if os.path.exists(rules_rall) and os.path.exists(compounds) and os.path.exists(rxn_recipes):
            shutil.copy(rules_rall, tmpOutputFolder+'/rules_rall.tsv')
            shutil.copy(compounds, tmpOutputFolder+'/compounds.tsv')
            shutil.copy(rxn_recipes, tmpOutputFolder+'/rxn_recipes.tsv')
            command = ['/home/tool_rpCofactors.py',
                       '-input',
                       '/home/tmp_output/input.dat',
                       '-output',
                       '/home/tmp_output/output.dat',
                       '-rxn_recipes',
                       '/home/tmp_output/rxn_recipes.tsv',
                       '-rules_rall',
                       '/home/tmp_output/rules_rall.tsv',
                       '-compounds',
                       '/home/tmp_output/compounds.tsv',
                       '-input_format',
                       str(input_format),
                       '-pathway_id',
                       str(pathway_id),
                       '-compartment_id',
                       str(compartment_id),
                       '-pubchem_search',
                       str(pubchem_search)]
        else:
            command = ['/home/tool_rpCofactors.py',
                       '-input',
                       '/home/tmp_output/input.dat',
                       '-output',
                       '/home/tmp_output/output.dat',
                       '-rxn_recipes',
                       'None',
                       '-rules_rall',
                       'None',
                       '-compounds',
                       'None',
                       '-input_format',
                       str(input_format),
                       '-pathway_id',
                       str(pathway_id),
                       '-compartment_id',
                       str(compartment_id),
                       '-pubchem_search',
                       str(pubchem_search)]
        container = docker_client.containers.run(image_str, 
                                                 command, 
                                                 detach=True, 
                                                 stderr=True,
                                                 volumes={tmpOutputFolder+'/': {'bind': '/home/tmp_output', 'mode': 'rw'}})
        container.wait()
        err = container.logs(stdout=True, stderr=True)
        err_str = err.decode('utf-8')
        if not 'ERROR' in err_str:
            shutil.copy(tmpOutputFolder+'/output.dat', output)
            logging.info('\n'+err_str)
        else:
            print(err_str)
        container.remove()


##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Add the missing cofactors to the monocomponent reactions to the SBML outputs of rpReader')
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
    main(params.input,
         params.input_format,
         params.rxn_recipes,
         params.rules_rall,
         params.compounds,
         params.output,
         params.pathway_id,
         params.compartment_id,
         params.pubchem_search)
