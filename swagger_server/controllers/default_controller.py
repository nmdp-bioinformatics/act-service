import connexion
from swagger_server.models.allele_call import AlleleCall
from swagger_server.models.error import Error
from datetime import date, datetime
from typing import List, Dict
from six import iteritems
from ..util import deserialize_date, deserialize_datetime

from py2neo import Graph
import sys
from gfe_typing.gfe import GfeTyping
from Bio import SeqIO
import re
from collections import defaultdict


def hla_post(locus, sequence, neo4j_url=None, gfe_url=None, verbose=None):
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
    user = "dash5"
    password = "chori"
    graph = Graph("http://neo4j.b12x.org:80", user=user, password=password, bolt=False)
    typer = GfeTyping(graph)
    gfe_output = typer.type_gfe(locus, sequence.upper())
    [gfe_new, hla, gfe] = [gfe_output[0], gfe_output[1], gfe_output[2]]
    act = {"gfe": gfe_new, "hla": hla}
    return act
