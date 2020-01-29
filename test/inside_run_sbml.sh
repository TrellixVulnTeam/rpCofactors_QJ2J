#!/bin/bash

python tool_rpCofactors.py -sbml test_rpsbml.xml -outputTar test_output.tar -pathway_id rp_pathway -compartment_id MNXC3
mv test_output.tar results/
