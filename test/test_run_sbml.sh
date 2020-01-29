#!/bin/bash

docker run -v ${PWD}/inside_run_sbml.sh:/home/inside_run_sbml.sh -v ${PWD}/tool_rpCofactors.py:/home/tool_rpCofactors.py -v ${PWD}/test_rpsbml.xml:/home/test_rpsbml.xml -v ${PWD}/results/:/home/results/ --rm brsynth/rpcofactors /bin/sh /home/inside_run_sbml.sh

cp results/test_output.tar .
