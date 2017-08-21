# swagger_client.DefaultApi

All URIs are relative to *https://localhost*

Method | HTTP request | Description
------------- | ------------- | -------------
[**act_get**](DefaultApi.md#act_get) | **GET** /act | 
[**actformat_get**](DefaultApi.md#actformat_get) | **GET** /act_format | 
[**ars_get**](DefaultApi.md#ars_get) | **GET** /ars | 
[**feature_get**](DefaultApi.md#feature_get) | **GET** /feature_search | 
[**gfe_get**](DefaultApi.md#gfe_get) | **GET** /gfe | 


# **act_get**
> AlleleCall act_get(locus, sequence=sequence, gfe=gfe, neo4j_url=neo4j_url, user=user, password=password, gfe_url=gfe_url, verbose=verbose, persist=persist)



Get HLA and GFE from consensus sequence or GFE notation

### Example 
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

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **locus** | **str**| Valid HLA locus | 
 **sequence** | **str**| Consensus sequence | [optional] 
 **gfe** | **str**| GFE Notation | [optional] 
 **neo4j_url** | **str**| URL for the neo4j graph | [optional] 
 **user** | **str**| Username for the neo4j graph | [optional] 
 **password** | **str**| Password for the neo4j graph | [optional] 
 **gfe_url** | **str**| URL for the gfe-service | [optional] 
 **verbose** | **bool**| Flag for running service in verbose | [optional] 
 **persist** | **bool**| Flag for persisting the data in the GFE DB | [optional] 

### Return type

[**AlleleCall**](AlleleCall.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: application/json

[[Back to top]](#) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to Model list]](../README.md#documentation-for-models) [[Back to README]](../README.md)

# **actformat_get**
> AlleleCall actformat_get(locus, sequence=sequence, gfe=gfe, format_type=format_type, neo4j_url=neo4j_url, user=user, password=password, gfe_url=gfe_url, verbose=verbose, persist=persist)



Get HLA and GFE from consensus sequence or GFE notation

### Example 
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
format_type = 'format_type_example' # str | Specify the data format that will be returned (optional)
neo4j_url = 'neo4j_url_example' # str | URL for the neo4j graph (optional)
user = 'user_example' # str | Username for the neo4j graph (optional)
password = 'password_example' # str | Password for the neo4j graph (optional)
gfe_url = 'gfe_url_example' # str | URL for the gfe-service (optional)
verbose = true # bool | Flag for running service in verbose (optional)
persist = true # bool | Flag for persisting the data in the GFE DB (optional)

try: 
    api_response = api_instance.actformat_get(locus, sequence=sequence, gfe=gfe, format_type=format_type, neo4j_url=neo4j_url, user=user, password=password, gfe_url=gfe_url, verbose=verbose, persist=persist)
    pprint(api_response)
except ApiException as e:
    print("Exception when calling DefaultApi->actformat_get: %s\n" % e)
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **locus** | **str**| Valid HLA locus | 
 **sequence** | **str**| Consensus sequence | [optional] 
 **gfe** | **str**| GFE Notation | [optional] 
 **format_type** | **str**| Specify the data format that will be returned | [optional] 
 **neo4j_url** | **str**| URL for the neo4j graph | [optional] 
 **user** | **str**| Username for the neo4j graph | [optional] 
 **password** | **str**| Password for the neo4j graph | [optional] 
 **gfe_url** | **str**| URL for the gfe-service | [optional] 
 **verbose** | **bool**| Flag for running service in verbose | [optional] 
 **persist** | **bool**| Flag for persisting the data in the GFE DB | [optional] 

### Return type

[**AlleleCall**](AlleleCall.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: text/plain

[[Back to top]](#) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to Model list]](../README.md#documentation-for-models) [[Back to README]](../README.md)

# **ars_get**
> ArsCall ars_get(allele, group, neo4j_url=neo4j_url, user=user, password=password, gfe_url=gfe_url, verbose=verbose)



Get ARS group associated with a GFE notation or HLA allele

### Example 
```python
from __future__ import print_function
import time
import swagger_client
from swagger_client.rest import ApiException
from pprint import pprint

# create an instance of the API class
api_instance = swagger_client.DefaultApi()
allele = 'allele_example' # str | HLA allele or GFE Notation
group = 'group_example' # str | ARS Group Type
neo4j_url = 'neo4j_url_example' # str | URL for the neo4j graph (optional)
user = 'user_example' # str | Username for the neo4j graph (optional)
password = 'password_example' # str | Password for the neo4j graph (optional)
gfe_url = 'gfe_url_example' # str | URL for the gfe-service (optional)
verbose = true # bool | Flag for running service in verbose (optional)

try: 
    api_response = api_instance.ars_get(allele, group, neo4j_url=neo4j_url, user=user, password=password, gfe_url=gfe_url, verbose=verbose)
    pprint(api_response)
except ApiException as e:
    print("Exception when calling DefaultApi->ars_get: %s\n" % e)
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **allele** | **str**| HLA allele or GFE Notation | 
 **group** | **str**| ARS Group Type | 
 **neo4j_url** | **str**| URL for the neo4j graph | [optional] 
 **user** | **str**| Username for the neo4j graph | [optional] 
 **password** | **str**| Password for the neo4j graph | [optional] 
 **gfe_url** | **str**| URL for the gfe-service | [optional] 
 **verbose** | **bool**| Flag for running service in verbose | [optional] 

### Return type

[**ArsCall**](ArsCall.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

[[Back to top]](#) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to Model list]](../README.md#documentation-for-models) [[Back to README]](../README.md)

# **feature_get**
> FeatureCall feature_get(hla, feature, neo4j_url=neo4j_url, user=user, password=password, gfe_url=gfe_url, verbose=verbose)



GFE notation and HLA alleles associated with an HLA allele or alleles

### Example 
```python
from __future__ import print_function
import time
import swagger_client
from swagger_client.rest import ApiException
from pprint import pprint

# create an instance of the API class
api_instance = swagger_client.DefaultApi()
hla = ['hla_example'] # list[str] | HLA Allele
feature = ['feature_example'] # list[str] | HLA feature
neo4j_url = 'neo4j_url_example' # str | URL for the neo4j graph (optional)
user = 'user_example' # str | Username for the neo4j graph (optional)
password = 'password_example' # str | Password for the neo4j graph (optional)
gfe_url = 'gfe_url_example' # str | URL for the gfe-service (optional)
verbose = true # bool | Flag for running service in verbose (optional)

try: 
    api_response = api_instance.feature_get(hla, feature, neo4j_url=neo4j_url, user=user, password=password, gfe_url=gfe_url, verbose=verbose)
    pprint(api_response)
except ApiException as e:
    print("Exception when calling DefaultApi->feature_get: %s\n" % e)
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **hla** | [**list[str]**](str.md)| HLA Allele | 
 **feature** | [**list[str]**](str.md)| HLA feature | 
 **neo4j_url** | **str**| URL for the neo4j graph | [optional] 
 **user** | **str**| Username for the neo4j graph | [optional] 
 **password** | **str**| Password for the neo4j graph | [optional] 
 **gfe_url** | **str**| URL for the gfe-service | [optional] 
 **verbose** | **bool**| Flag for running service in verbose | [optional] 

### Return type

[**FeatureCall**](FeatureCall.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

[[Back to top]](#) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to Model list]](../README.md#documentation-for-models) [[Back to README]](../README.md)

# **gfe_get**
> GfeCall gfe_get(hla, neo4j_url=neo4j_url, user=user, password=password, gfe_url=gfe_url, verbose=verbose)



Get GFE notation associated with an HLA allele

### Example 
```python
from __future__ import print_function
import time
import swagger_client
from swagger_client.rest import ApiException
from pprint import pprint

# create an instance of the API class
api_instance = swagger_client.DefaultApi()
hla = 'hla_example' # str | HLA allele
neo4j_url = 'neo4j_url_example' # str | URL for the neo4j graph (optional)
user = 'user_example' # str | Username for the neo4j graph (optional)
password = 'password_example' # str | Password for the neo4j graph (optional)
gfe_url = 'gfe_url_example' # str | URL for the gfe-service (optional)
verbose = true # bool | Flag for running service in verbose (optional)

try: 
    api_response = api_instance.gfe_get(hla, neo4j_url=neo4j_url, user=user, password=password, gfe_url=gfe_url, verbose=verbose)
    pprint(api_response)
except ApiException as e:
    print("Exception when calling DefaultApi->gfe_get: %s\n" % e)
```

### Parameters

Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **hla** | **str**| HLA allele | 
 **neo4j_url** | **str**| URL for the neo4j graph | [optional] 
 **user** | **str**| Username for the neo4j graph | [optional] 
 **password** | **str**| Password for the neo4j graph | [optional] 
 **gfe_url** | **str**| URL for the gfe-service | [optional] 
 **verbose** | **bool**| Flag for running service in verbose | [optional] 

### Return type

[**GfeCall**](GfeCall.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: Not defined

[[Back to top]](#) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to Model list]](../README.md#documentation-for-models) [[Back to README]](../README.md)

