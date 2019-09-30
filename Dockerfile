FROM ibisba/rpsbml

RUN conda install -y -c rdkit rdkit && \
    mkdir input_cache && \
    wget https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_xref.tsv -P /home/input_cache/ && \
    wget https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_prop.tsv -P /home/input_cache/ && \
    wget https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_xref.tsv -P /home/input_cache/

#need to replace this with urls
COPY required_data/input_cache/rules_rall.tsv /home/input_cache/
COPY required_data/input_cache/compounds.tsv /home/input_cache/rr_compounds.tsv 
COPY required_data/input_cache/rxn_recipes.tsv /home/input_cache/

COPY rpCofactors.py /home/

#build the cache -- not sure if needed
RUN python rpCofactors.py




#wget rules_rall_TODO -P /home/input_cache/ && \
#wget rr_compounds_TODO -P /home/input_cache/ && \
#wget rxn_recipes_TODO -P /home/input_cache/
