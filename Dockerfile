FROM brsynth/rpcache

COPY rpTool.py /home/
COPY rpToolServe.py /home/
COPY rpToolCache.py /home/
COPY tool_rpCofactors.py /home/

RUN python rpToolCache.py
