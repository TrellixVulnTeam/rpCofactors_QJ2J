#!/bin/sh

rm -f test_output.tar
docker run -d -p 5000:5000 --name test_rpCofactors brsynth/rpcofactors
sleep 10
python tool_rpCofactors.py -inputTar test_input.tar -outputTar test_output.tar -pathway_id rp_pathway -compartment_id MNXC3 -server_url http://0.0.0.0:5000/REST
docker kill test_rpCofactors
docker rm test_rpCofactors
