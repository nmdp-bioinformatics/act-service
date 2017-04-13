import connexion
from swagger_server.models.allele_call import AlleleCall
from swagger_server.models.error import Error
from datetime import date, datetime
from typing import List, Dict
from six import iteritems
from ..util import deserialize_date, deserialize_datetime

from py2neo import Graph
from gfe_typing.gfe import GfeTyping
import os

neo4j_pass = ''
if os.getenv("NEO4JPASS"):
    neo4j_pass = os.getenv("NEO4JPASS")

neo4j_user = ''
if os.getenv("NEO4JUSER"):
    neo4j_user = os.getenv("NEO4JUSER")


def hla_post(locus, sequence, neo4j_url="http://neo4j.b12x.org:80", gfe_url="gfe.b12x.org", user=neo4j_user, password=neo4j_pass, verbose=None):
    """
    hla_post
    Get HLA and GFE from consensus sequence
    :param locus: Valid HLA locus
    :type locus: str
    :param sequence: Consensus sequence
    :type sequence: str
    :param neo4j_url: URL for the neo4j graph
    :type neo4j_url: str
    :param gfe_url: URL for the gfe-service
    :type gfe_url: str
    :param verbose: Flag for running service in verbose
    :type verbose: bool

    :rtype: AlleleCall
    """
    graph = Graph(neo4j_url, user=user, password=password, bolt=False)
    typer = GfeTyping(graph, hostname=gfe_url)
    gfe_output = typer.type_gfe(locus, sequence.upper())
    [gfe_new, hla, gfe] = [gfe_output[0], gfe_output[1], gfe_output[2]]
    act = {"gfe": gfe_new, "hla": hla}
    return act
