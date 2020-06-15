FROM brsynth/rpcache:newrules

#RUN wget https://raw.githubusercontent.com/Galaxy-SynBioCAD/rpUnicity/standalone/code/rpUnicity.py .
RUN git clone https://github.com/Galaxy-SynBioCAD/rpUnicity.git -b rest-dev
RUN mv rpUnicity/code/tool_rpUnicity.py .
RUN rm -r rpUnicity
COPY rpTool.py /home/
COPY rpToolServe.py /home/
COPY tool_rpCofactors.py /home/
