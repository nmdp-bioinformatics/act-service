# swagger-client
The Allele Calling  service provides an API for converting raw sequence data to GFE and HLA allele calls.

This Python package is automatically generated by the [Swagger Codegen](https://github.com/swagger-api/swagger-codegen) project:

- API version: 0.0.3
- Package version: 1.0.0
- Build package: io.swagger.codegen.languages.PythonClientCodegen

## Requirements.

Python 2.7 and 3.4+

## Installation & Usage
### pip install

If the python package is hosted on Github, you can install directly from Github

```sh
pip install git+https://github.com//.git
```
(you may need to run `pip` with root permission: `sudo pip install git+https://github.com//.git`)

Then import the package:
```python
import swagger_client 
```

### Setuptools

Install via [Setuptools](http://pypi.python.org/pypi/setuptools).

```sh
python setup.py install --user
```
(or `sudo python setup.py install` to install the package for all users)

Then import the package:
```python
import swagger_client
```

## Getting Started

Please follow the [installation procedure](#installation--usage) and then run the following:

```python
from __future__ import print_function
import time
import swagger_client
from swagger_client.rest import ApiException
from pprint import pprint
# create an instance of the API class
api_instance = swagger_client.DefaultApi()
locus = 'locus_example' # str | Valid HLA locus
sequence = 'sequence_example' # str | Consensus sequence (optional)
gfe = 'gfe_example' # str | GFE Notation (optional)
neo4j_url = 'neo4j_url_example' # str | URL for the neo4j graph (optional)
user = 'user_example' # str | Username for the neo4j graph (optional)
password = 'password_example' # str | Password for the neo4j graph (optional)
gfe_url = 'gfe_url_example' # str | URL for the gfe-service (optional)
verbose = true # bool | Flag for running service in verbose (optional)
persist = true # bool | Flag for persisting the data in the GFE DB (optional)

try:
    api_response = api_instance.act_get(locus, sequence=sequence, gfe=gfe, neo4j_url=neo4j_url, user=user, password=password, gfe_url=gfe_url, verbose=verbose, persist=persist)
    pprint(api_response)
except ApiException as e:
    print("Exception when calling DefaultApi->act_get: %s\n" % e)

```

## Documentation for API Endpoints

All URIs are relative to *https://localhost*

Class | Method | HTTP request | Description
------------ | ------------- | ------------- | -------------
*DefaultApi* | [**act_get**](docs/DefaultApi.md#act_get) | **GET** /act | 
*DefaultApi* | [**actformat_get**](docs/DefaultApi.md#actformat_get) | **GET** /act_format | 
*DefaultApi* | [**ars_get**](docs/DefaultApi.md#ars_get) | **GET** /ars | 
*DefaultApi* | [**feature_get**](docs/DefaultApi.md#feature_get) | **GET** /feature_search | 
*DefaultApi* | [**gfe_get**](docs/DefaultApi.md#gfe_get) | **GET** /gfe | 


## Documentation For Models

 - [AlleleCall](docs/AlleleCall.md)
 - [ArsCall](docs/ArsCall.md)
 - [Error](docs/Error.md)
 - [Feature](docs/Feature.md)
 - [FeatureCall](docs/FeatureCall.md)
 - [GfeCall](docs/GfeCall.md)
 - [GfeTyping](docs/GfeTyping.md)
 - [Typing](docs/Typing.md)


## Documentation For Authorization

 All endpoints do not require authorization.


## Author

mhalagan@nmdp.org
