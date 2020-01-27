#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Galaxy script to query rpCofactors REST service

"""
import json
import argparse
import requests

def rpCofactorsUpload(inputTar,
        pathway_id,
        compartment_id,
        server_url,
        outputTar):
    # Post request
    data = {'pathway_id': pathway_id, 'compartment_id': compartment_id}
    files = {'inputTar': open(inputTar, 'rb'),
             'data': ('data.json', json.dumps(data))}
    r = requests.post(server_url+'/Query', files=files)
    r.raise_for_status()
    with open(outputTar, 'wb') as ot:
        ot.write(r.content)

##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Python wrapper to add cofactors to generate rpSBML collection')
    parser.add_argument('-inputTar', type=str)
    parser.add_argument('-pathway_id', type=str)
    parser.add_argument('-compartment_id', type=str)
    parser.add_argument('-outputTar', type=str)
    parser.add_argument('-server_url', type=str)
    params = parser.parse_args()
    rpCofactorsUpload(params.inputTar,
                      params.pathway_id,
                      params.compartment_id,
                      params.server_url,
                      params.outputTar) 
