FROM brsynth/rpcache:dev

COPY rpToolCache.py /home/
RUN python rpToolCache.py

RUN wget https://raw.githubusercontent.com/Galaxy-SynBioCAD/rpUnicity/standalone/code/rpUnicity.py .
COPY rpTool.py /home/
COPY rpToolServe.py /home/
COPY tool_rpCofactors.py /home/
