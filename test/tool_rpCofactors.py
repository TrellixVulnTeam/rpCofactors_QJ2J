#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Galaxy script to query rpCofactors REST service

"""
import requests
import argparse
import json
import glob
import tarfile

## Use
#
#
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

## run using HDD 3X less than the above function
#
#
def runCofactors_hdd(inputTar, outputTar, pathId='rp_pathway', compartment_id='MNXC3'):
    rpcofactors = rpCofactors.rpCofactors()
    if not os.path.exists(os.getcwd()+'/tmp'):
        os.mkdir(os.getcwd()+'/tmp')
    tmpInputFolder = os.getcwd()+'/tmp/'+''.join(random.choice(string.ascii_lowercase) for i in range(15))
    tmpOutputFolder = os.getcwd()+'/tmp/'+''.join(random.choice(string.ascii_lowercase) for i in range(15))
    os.mkdir(tmpInputFolder)
    os.mkdir(tmpOutputFolder)
    tar = tarfile.open(inputTar, 'r:xz')
    tar.extractall(path=tmpInputFolder)
    tar.close()
    for sbml_path in glob.glob(tmpInputFolder+'/*'):
        fileName = sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '')
        rpsbml = rpSBML.rpSBML(fileName)
        rpsbml.readSBML(sbml_path)
        rpcofactors.addCofactors(rpsbml, compartment_id, pathId)
        rpsbml.writeSBML(tmpOutputFolder)
        rpsbml = None
    with tarfile.open(outputTar, mode='w:xz') as ot:
        for sbml_path in glob.glob(tmpOutputFolder+'/*'):
            fileName = str(sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', ''))
            info = tarfile.TarInfo(fileName)
            info.size = os.path.getsize(sbml_path)
            ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))
    shutil.rmtree(tmpInputFolder)
    shutil.rmtree(tmpOutputFolder)


##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Python wrapper to add cofactors to generate rpSBML collection')
    parser.add_argument('-inputTar', type=str)
    parser.add_argument('-pathway_id', type=str)
    parser.add_argument('-compartment_id', type=str)
    parser.add_argument('-server_url', type=str)
    parser.add_argument('-outputTar', type=str)
    params = parser.parse_args()
    rpCofactorsUpload(params.inputTar,
            params.pathway_id,
            params.compartment_id,
            params.server_url,
            params.outputTar)
