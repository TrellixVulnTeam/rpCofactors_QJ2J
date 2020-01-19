#FROM brsynth/rpcache
FROM brsynth/rpbase

COPY rpTool.py /home/
COPY rpToolServe.py /home/
COPY rpToolCache.py /home/

RUN python rpToolCache.py
