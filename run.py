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
def main(inputfile, input_format, output, pathway_id='rp_pathway', compartment_id='MNXC3'):
    docker_client = docker.from_env()
    image_str = 'brsynth/rpcofactors-standalone:dev'
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
        command = ['/home/tool_rpCofactors.py',
                   '-input',
                   '/home/tmp_output/input.dat',
                   '-output',
                   '/home/tmp_output/output.dat',
                   '-input_format',
                   input_format,
                   '-pathway_id',
                   pathway_id,
                   '-compartment_id',
                   compartment_id]
        container = docker_client.containers.run(image_str, 
                                                 command, 
                                                 detach=True, 
                                                 stderr=True,
                                                 volumes={tmpOutputFolder+'/': {'bind': '/home/tmp_output', 'mode': 'rw'}})
        container.wait()
        err = container.logs(stdout=False, stderr=True)
        err_str = err.decode('utf-8')
        print(err_str)
        if not 'ERROR' in err_str:
            shutil.copy(tmpOutputFolder+'/output.dat', output)
        container.remove()


##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Add the missing cofactors to the monocomponent reactions to the SBML outputs of rpReader')
    parser.add_argument('-input', type=str)
    parser.add_argument('-output', type=str)
    parser.add_argument('-input_format', type=str)
    parser.add_argument('-pathway_id', type=str, default='rp_pathway')
    parser.add_argument('-compartment_id', type=str, default='MNXC3')
    params = parser.parse_args()
    main(params.input, params.input_format, params.output, params.pathway_id, params.compartment_id)
