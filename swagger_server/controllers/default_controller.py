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

neo4jurl = "http://neo4j.b12x.org:80"
if os.getenv("NEO4JURL"):
    neo4jurl = os.getenv("NEO4JURL")

gfeurl = "http://localhost:3000"
if os.getenv("GFEURL"):
    gfeurl = os.getenv("GFEURL")


def actformat_get(locus, sequence=None, neo4j_url=neo4jurl, user=neo4juser, password=neo4jpass, gfe_url=gfeurl, format_type=None, gfe=None, verbose=None, persist=None):
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
    typer = Act(graph, hostname=gfeurl)
    allele_call = typer.type_hla(locus, sequence, gfe)

    if isinstance(allele_call, Error):
        return allele_call, 404
    else:
        imgt_formatted = typer.typing_to_bioseq(allele_call, sequence)
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
    typer = Act(graph, hostname=gfeurl)
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
    typer = Act(graph, hostname=gfeurl)
    ars_call = typer.ars_redux(group, allele)
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
    typer = Act(graph, hostname=gfe_url)
    gfe_output = typer.get_gfe(hla)
    allele_call = AlleleCall(gfe=hla, hla=gfe_output, version='0.0.1')
    return allele_call


def feature_get(hla, feature, neo4j_url=None, user=None, password=None, gfe_url=None, verbose=None):
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
    return 'do some magic!'


def seqsrch_get(locus, start=None, end=None, hla=None, feature=None, neo4j_url=None, user=None, password=None, gfe_url=None, verbose=None):
    """
    seqsrch_get
    Sequence data associated from a specific feature or HLA allele
    :param locus: HLA locus
    :type locus: str
    :param start: Starting point of sequence
    :type start: int
    :param end: End point of sequence
    :type end: int
    :param hla: HLA Allele
    :type hla: List[str]
    :param feature: HLA feature rank
    :type feature: str
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

    :rtype: InlineResponse2002
    """
    return 'do some magic!'
