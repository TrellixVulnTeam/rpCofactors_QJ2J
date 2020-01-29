#!/bin/bash

docker run -v ${PWD}/inside_run_tar.sh:/home/inside_run_tar.sh -v ${PWD}/tool_rpCofactors.py:/home/tool_rpCofactors.py -v ${PWD}/test_input.tar:/home/test_input.tar -v ${PWD}/results/:/home/results/ --rm brsynth/rpcofactors /bin/sh /home/inside_run_tar.sh

cp results/test_output.tar .
