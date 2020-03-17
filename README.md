# rpCofactors

Because the reactions calculated by RetroPath2.0 are monocomponent reactions (only one product with possibly multiple substrates), to perform a full analysis, the reactions need to be completed with the missing species. This tool performs that action this by referring to the reaction description from wich it was generated. The following tool is a REST service that takes as input a tar.xz of rpSBML files (output by rpReader).

## Information Flow

### Input

Required information:
* **rpSBML input**: Either tar.xz input collection of rpSBML files or a single rpSBML file.

Advanced options:
* **Name of the heterologous pathway**: (default: rp_pathway) The SBML groups ID (defined in rpReader) that points to the heterologous reactions and chemical species.
* **SBML compartment ID**: (default: MNXC3) Compartment ID to add the new chemical species. The default is the cytoplasm.
* **REST IP address**: The IP addrress of the REST service

### Output

* **rpCofactors**: A tar.xz archive containing a list of rpSBML files.

## Installing

To build the image using the Dockerfile, use the following command:

```
docker build -t brsynth/rpcofactors-rest .
```

To run the service in the localhost:

```
docker run -p 8885:8888 brsynth/rpcofactors-rest
```

### Prerequisites

* [Docker](https://docs.docker.com/v17.09/engine/installation/)
* [libSBML](http://sbml.org/Software/libSBML)
* [RDkit](https://www.rdkit.org)

## Contributing

TODO

## Versioning

Version 0.1

## Authors

* **Melchior du Lac**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Thomas Duigou
* Joan Hérisson

### How to cite rpCofactors?

TODO
