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
from gfe_db.graph import GfeDB
import os

neo4jpass = ''
if os.getenv("NEO4JPASS"):
    neo4jpass = os.getenv("NEO4JPASS")

neo4juser = ''
if os.getenv("NEO4JUSER"):
    neo4juser = os.getenv("NEO4JUSER")

neo4jurl = "http://neo4j.b12x.org:80"
if os.getenv("NEO4JURL"):
    neo4jurl = os.getenv("NEO4JURL")

gfeurl = "http://gfe.b12x.org"
if os.getenv("GFEURL"):
    gfeurl = os.getenv("GFEURL")


def act_get(locus, sequence, neo4j_url=neo4jurl, user=neo4juser, password=neo4jpass, gfe_url=gfeurl, verbose=None):
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
    typer = GfeDB(graph, hostname=gfe_url)
    gfe_output = typer.type_hla(locus, sequence.upper())
    [gfe_new, hla, gfe] = [gfe_output[0], gfe_output[1], gfe_output[2]]
    allele_call = AlleleCall(gfe=gfe, hla=hla, version='0.0.1')
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
    typer = GfeDB(graph, hostname=gfe_url)
    ars_output = typer.ars_redux(group, allele)
    ars_call = ArsCall(allele=allele, group_type=group, group=ars_output, version='0.0.1')
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
    typer = GfeDB(graph, hostname=gfe_url)
    gfe_output = typer.get_gfe(hla)
    allele_call = AlleleCall(gfe=hla, hla=gfe_output, version='0.0.1')
    return allele_call


def hla_get(gfe, neo4j_url="http://neo4j.b12x.org:80", user=neo4juser, password=neo4jpass, gfe_url=gfeurl, verbose=None):
    """
    hla_get
    Get HLA associated with GFE notation
    :param gfe: GFE Notation
    :type gfe: str
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
    typer = GfeDB(graph, hostname=gfe_url)
    hla_output = typer.get_hla(gfe)
    hla_call = AlleleCall(gfe=gfe, hla=hla_output, version='0.0.1')
    return hla_call


def sequence_get(allele, allele_type, neo4j_url=neo4jurl, user=neo4juser, password=neo4jpass, gfe_url=gfeurl, verbose=None):
    """
    sequence_get
    Get sequence associated with an HLA allele or GFE notation
    :param allele: HLA allele or GFE notation
    :type allele: str
    :param allele_type: Specify whether it&#39;s IMGT or GFE
    :type allele_type: str
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

    :rtype: Sequence
    """
    graph = Graph(neo4j_url, user=user, password=password, bolt=False)
    typer = GfeDB(graph, hostname=gfe_url)
    seq_output = typer.sequence(allele_type, allele)
    seq = Sequence(sequence=seq_output, version='0.0.1')
    return seq
