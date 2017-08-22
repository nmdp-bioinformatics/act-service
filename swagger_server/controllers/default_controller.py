import connexion
from swagger_server.models.allele_call import AlleleCall
from swagger_server.models.ars_call import ArsCall
from swagger_server.models.error import Error
from swagger_server.models.sequence import Sequence
from datetime import date, datetime
from typing import List, Dict
from six import iteritems
from ..util import deserialize_date, deserialize_datetime
from py2neo import Graph
from gfe_db.act import Act
import os

from io import StringIO
from Bio import SeqIO


neo4jpass = ''
if os.getenv("NEO4JPASS"):
    neo4jpass = os.getenv("NEO4JPASS")

neo4juser = ''
if os.getenv("NEO4JUSER"):
    neo4juser = os.getenv("NEO4JUSER")

neo4jurl = "http://localhost:7474"
if os.getenv("NEO4JURL"):
    neo4jurl = os.getenv("NEO4JURL")

gfeurl = "http://gfe.b12x.org"
if os.getenv("GFEURL"):
    gfeurl = os.getenv("GFEURL")


def actformat_get(locus, format_type, sequence=None, neo4j_url=neo4jurl, user=neo4juser, password=neo4jpass, gfe_url=gfeurl, gfe=None, verbose=None, persist=None):
    """
    act_get
    Get HLA and GFE from consensus sequence
    :param locus: Valid HLA locus
    :type locus: str
    :param sequence: Consensus sequence
    :type sequence: str
    :param neo4j_url: URL for the neo4j graph
    :type neo4j_url: str
    :param user: Username for the neo4j graph
    :type user: str
    :param pass: Password for the neo4j graph
    :type pass: str
    :param gfe_url: URL for the gfe-service
    :type gfe_url: str
    :param verbose: Flag for running service in verbose
    :type verbose: bool

    :rtype: AlleleCall
    """
    graph = Graph(neo4j_url, user=user, password=password, bolt=False)
    typer = Act(graph, hostname=gfeurl, user=neo4juser, persist=persist)
    allele_call = typer.type_hla(locus, sequence, gfe)

    if isinstance(allele_call, Error):
        return allele_call, 404
    else:
        imgt_formatted = typer.typing_to_bioseq(allele_call)
        imgt_fh = StringIO()
        SeqIO.write(imgt_formatted, imgt_fh, format_type)
        imgt_data = imgt_fh.getvalue()
        return imgt_data, 200, {'content-type': 'text/plain' }


def act_get(locus, sequence=None, neo4j_url=neo4jurl, user=neo4juser, password=neo4jpass, gfe_url=gfeurl, gfe=None, verbose=None, persist=None):
    """
    act_get
    Get HLA and GFE from consensus sequence
    :param locus: Valid HLA locus
    :type locus: str
    :param sequence: Consensus sequence
    :type sequence: str
    :param neo4j_url: URL for the neo4j graph
    :type neo4j_url: str
    :param user: Username for the neo4j graph
    :type user: str
    :param pass: Password for the neo4j graph
    :type pass: str
    :param gfe_url: URL for the gfe-service
    :type gfe_url: str
    :param verbose: Flag for running service in verbose
    :type verbose: bool

    :rtype: AlleleCall
    """
    graph = Graph(neo4j_url, user=user, password=password, bolt=False)
    typer = Act(graph, hostname=gfeurl, user=neo4juser, persist=persist)
    allele_call = typer.type_hla(locus, sequence, gfe)
    if isinstance(allele_call, Error):
        return allele_call, 404
    else:
        return allele_call


def ars_get(allele, group, neo4j_url=neo4jurl, user=neo4juser, password=neo4jpass, gfe_url=gfeurl, verbose=None):
    """
    ars_get
    Get ARS group associated with a GFE notation or HLA allele
    :param allele: HLA allele or GFE Notation
    :type allele: str
    :param neo4j_url: URL for the neo4j graph
    :type neo4j_url: str
    :param user: Username for the neo4j graph
    :type user: str
    :param pass: Password for the neo4j graph
    :type pass: str
    :param gfe_url: URL for the gfe-service
    :type gfe_url: str
    :param verbose: Flag for running service in verbose
    :type verbose: bool

    :rtype: ArsCall
    """
    graph = Graph(neo4j_url, user=user, password=password, bolt=False)
    typer = Act(graph, hostname=gfeurl, user=neo4juser, persist=None)
    ars_call = typer.ars_redux(group, allele)
    if isinstance(ars_call, Error):
        return ars_call, 404
    else:
        return ars_call


def gfe_get(hla, neo4j_url=neo4jurl, user=neo4juser, password=neo4jpass, gfe_url=gfeurl, verbose=None):
    """
    gfe_get
    Get GFE notation associated with an HLA allele
    :param hla: HLA allele
    :type hla: str
    :param neo4j_url: URL for the neo4j graph
    :type neo4j_url: str
    :param user: Username for the neo4j graph
    :type user: str
    :param pass: Password for the neo4j graph
    :type pass: str
    :param gfe_url: URL for the gfe-service
    :type gfe_url: str
    :param verbose: Flag for running service in verbose
    :type verbose: bool

    :rtype: AlleleCall
    """
    graph = Graph(neo4j_url, user=user, password=password, bolt=False)
    typer = Act(graph, hostname=gfe_url, user=neo4juser, persist=None)
    gfe_call = typer.get_gfe_call(hla)
    if isinstance(gfe_call, Error):
        return gfe_call, 404
    else:
        return gfe_call


def feature_get(hla, feature, neo4j_url=neo4jurl, user=neo4juser, password=neo4jpass, gfe_url=gfeurl, verbose=None):
    """
    feature_get
    GFE notation and HLA alleles associated with an HLA allele or alleles
    :param hla: HLA Allele
    :type hla: List[str]
    :param feature: HLA feature
    :type feature: List[str]
    :param neo4j_url: URL for the neo4j graph
    :type neo4j_url: str
    :param user: Username for the neo4j graph
    :type user: str
    :param password: Password for the neo4j graph
    :type password: str
    :param gfe_url: URL for the gfe-service
    :type gfe_url: str
    :param verbose: Flag for running service in verbose
    :type verbose: bool

    :rtype: InlineResponse2003
    """
    graph = Graph(neo4j_url, user=user, password=password, bolt=False)
    typer = Act(graph, hostname=gfeurl, user=neo4juser, persist=None)
    feature_call = typer.search_features(hla, feature)
    if isinstance(feature_call, Error):
        return feature_call, 404
    else:
        return feature_call

