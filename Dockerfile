FROM brsynth/rpcache

COPY rpToolCache.py /home/
RUN python rpToolCache.py

COPY rpTool.py /home/
COPY rpToolServe.py /home/
COPY tool_rpCofactors.py /home/
