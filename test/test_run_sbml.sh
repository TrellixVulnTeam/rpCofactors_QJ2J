#!/bin/sh

rm -f test_output.tar
docker run -d -p 8888:8888 --name test_rpCofactors brsynth/rpcofactors-rest
sleep 10
python tool_rpCofactors.py -sbml test_rpsbml.xml -outputTar test_output.tar -pathway_id rp_pathway -compartment_id MNXC3 -server_url http://0.0.0.0:8888/REST
docker kill test_rpCofactors
docker rm test_rpCofactors
