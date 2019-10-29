# rpCofactors REST

REST service of rpCofactors. Takes as input a tar.xz of rpSBML files and outputs the same with cofactors added to the rpSBML.

### Prerequisites


### Installing

```
docker build -t brsynth/rpcofactors -f Dockerfile .
```

Run the service

```
docker run --network host -p 8996:8996 brsynth/rpcofactors
```

## Running the tests

TODO

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With


## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

TODO

## Authors

* **Melchior du Lac** 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Thomas Duigou
* Joan Hérisson
