import connexion
import six
import os

from pygfe.models.error import Error
from neo4j.exceptions import ServiceUnavailable
from swagger_server.models.typing import Typing  # noqa: E501
from swagger_server import util
from py2neo import Graph
from pygfe.pygfe import pyGFE
from seqann.sequence_annotation import BioSeqAnn
import logging
import io
import re

from pandas import DataFrame

seqanns = {}
gfe_feats = None
gfe2hla = None
seq2hla = None


def typealign_get(sequence, locus=None, imgthla_version='3.31.0', neo4j_url="http://neo4j.b12x.org:80", user='neo4j', password='gfedb'):  # noqa: E501
    """typealign_get

    Get HLA and GFE from consensus sequence or GFE notation # noqa: E501

    :param locus: Valid HLA locus
    :type locus: str
    :param sequence: Consensus sequence
    :type sequence: str
    :param imgthla_version: IMGT/HLA DB Version
    :type imgthla_version: str
    :param neo4j_url: URL for the neo4j graph
    :type neo4j_url: str
    :param user: Username for the neo4j graph
    :type user: str
    :param password: Password for the neo4j graph
    :type password: str
    :param verbose: Flag for running service in verbose
    :type verbose: bool

    :rtype: Typing
    """
    global seqanns
    global gfe_feats
    global gfe2hla
    global seq2hla

    log_capture_string = io.StringIO()
    logger = logging.getLogger('')
    logging.basicConfig(datefmt='%m/%d/%Y %I:%M:%S %p',
                        level=logging.INFO)

    # create console handler and set level to debug
    ch = logging.StreamHandler(log_capture_string)
    formatter = logging.Formatter('%(asctime)s - %(name)-35s - %(levelname)-5s - %(funcName)s %(lineno)d: - %(message)s')
    ch.setFormatter(formatter)
    ch.setLevel(logging.INFO)
    logger.addHandler(ch)

    if not re.match(".", imgthla_version):
        imgthla_version = ".".join([list(imgthla_version)[0],
                                    "".join(list(imgthla_version)[1:3]),
                                    list(imgthla_version)[3]])

    db = "".join(imgthla_version.split("."))
    if db in seqanns:
        seqann = seqanns[db]
    else:
        seqann = BioSeqAnn(align=True, verbose=True, dbversion=db, verbosity=3)
        seqanns.update({db: seqann})

    try:
        graph = Graph(neo4j_url,
                      user=user,
                      password=password,
                      bolt=False)
    except ServiceUnavailable as err:
        log_contents = log_capture_string.getvalue()
        log_data = log_contents.split("\n")
        log_data.append(str(err))
        return Error("Failed to connect to graph", log=log_data), 404

    if(not isinstance(gfe_feats, DataFrame)
       or not isinstance(seq2hla, DataFrame)):
        pygfe = pyGFE(graph=graph, seqann=seqann,
                      load_gfe2hla=True, load_seq2hla=True,
                      load_gfe2feat=True, verbose=True)
        gfe_feats = pygfe.gfe_feats
        seq2hla = pygfe.seq2hla
        gfe2hla = pygfe.gfe2hla
    else:
        pygfe = pyGFE(graph=graph, seqann=seqann,
                      gfe2hla=gfe2hla,
                      gfe_feats=gfe_feats,
                      seq2hla=seq2hla,
                      verbose=True)

    try:
        typing = pygfe.type_from_seq(locus, sequence, imgthla_version)
    except:
        log_contents = log_capture_string.getvalue()
        return Error("Type with alignment failed", log=log_contents.split("\n")), 404

    if isinstance(typing, Error):
        log_contents = log_capture_string.getvalue()
        typing.log = log_contents.split("\n")
        return typing, 404

    if not typing:
        log_contents = log_capture_string.getvalue()
        return Error("Type with alignment failed", log=log_contents.split("\n")), 404

    typing.gfedb_version = "2.0.0"
    return typing
