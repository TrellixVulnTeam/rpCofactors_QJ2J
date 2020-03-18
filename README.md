# rpCofactors

* Docker Image: [brsynth/rpcofactors-standalone](https://hub.docker.com/r/brsynth/rpcofactors-standalone)

Completes monocomponent reaction output by RetroPath2.0 with the appropriate cofactors. Creates sub-paths when multiple reaction rules are associated with a single reaction. Input may be a single SBML file or a collection within a tar.xz archive

## Input

Required:
* **-input**: (string) Path to the input file. Can be either a single SBML file or a collection as a tar.xz archive file
* **-input_format**: (string) 

Advanced Options:
* **-pathway_id**: (string, default: rp_pathway) ID of the heterologous pathway
* **-compartment_id**: (string, default: MNXC3 (i.e. cytoplasm)) ID of the SBML compartment where the heterologous pathway will be expressed in (default: MNXC3 (i.e. cytoplasm))
* **server_url**: (string, default: http://0.0.0.0:8888/REST) IP address of the REST service

## Output

* **-output**: (string) Path to the output file

## Dependencies

* Base docker image: [brsynth/rpBase](https://hub.docker.com/r/brsynth/rpbase)
* Cache docker image: [brsynth/rpCache](https://hub.docker.com/r/brsynth/rpcache)

## Installing

To build the image using the Dockerfile, use the following command:

```
docker build -t brsynth/rpcofactors-rest:dev .
```

To run the service in the localhost:

```
docker run -p 8888:8888 brsynth/rpcofactors-rest:dev
```

### Running the test

Untar the test.tar.xz file and run the following command:

```
python tool_rpCofactors.py -input test_rpReader.tar -input_format tar -output test_rpCofactors.tar
```

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

v0.1

## Authors

* **Melchior du Lac**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Thomas Duigou
* Joan Hérisson

## How to cite rpCofactors?
