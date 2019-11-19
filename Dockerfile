FROM brsynth/rpcache

COPY rpCofactors.py /home/
COPY rpCofactorsServe.py /home/

ENTRYPOINT ["python"]
CMD ["/home/rpCofactorsServe.py"]

# Open server port
EXPOSE 8996
