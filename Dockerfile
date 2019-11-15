FROM brsynth/rpbase

RUN apt-get install --quiet --yes --no-install-recommends \
	libxext6  \
    	libxrender-dev  && \
    conda install -y -c rdkit rdkit && \
    conda install -c conda-forge flask-restful && \
    mkdir input_cache && \
    wget https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_xref.tsv -P /home/input_cache/ && \
    wget https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_prop.tsv -P /home/input_cache/ && \
    wget https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_xref.tsv -P /home/input_cache/

#get the rules
RUN wget https://retrorules.org/dl/preparsed/rr02/rp3/hs -O /home/rules_rall_rp3.tar.gz && \
    tar xf /home/rules_rall_rp3.tar.gz -C /home/ && \
    mv /home/retrorules_rr02_rp3_hs/retrorules_rr02_flat_all.tsv /home/input_cache/rules_rall.tsv && \
    rm -r /home/retrorules_rr02_rp3_hs && \
    rm /home/rules_rall_rp3.tar.gz

#get rr_compounds.tsv and rxn_recipes
RUN wget https://retrorules.org/dl/this/is/not/a/secret/path/rr02 -O /home/rr02_more_data.tar.gz && \
    tar xf /home/rr02_more_data.tar.gz -C /home/ && \
    mv /home/rr02_more_data/compounds.tsv /home/input_cache/rr_compounds.tsv && \
    mv /home/rr02_more_data/rxn_recipes.tsv /home/input_cache/ && \
    rm -r /home/rr02_more_data && \
    rm /home/rr02_more_data.tar.gz

COPY rpCofactors.py /home/
COPY rpCache.py /home/
COPY rpCofactorsServe.py /home/

RUN python /home/rpCache.py

ENTRYPOINT ["python"]
CMD ["/home/rpCofactorsServe.py"]

# Open server port
EXPOSE 8996
