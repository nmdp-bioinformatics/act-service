import connexion
import six

from pygfe.models.error import Error  # noqa: E501
from pygfe.models.seqdiff import Seqdiff  # noqa: E501
from pygfe.models.typing import Typing  # noqa: E501
from pygfe import util

from pygfe.cypher import all_feats

from py2neo import Graph
from pygfe.pygfe import pyGFE
from pygfe.gfe import GFE
import os

from io import StringIO
from Bio import SeqIO

from seqann.sequence_annotation import BioSeqAnn
from BioSQL import BioSeqDatabase
import pymysql
import pickle
import pandas as pd

neo4jpass = 'gfedb'
if os.getenv("NEO4JPASS"):
    neo4jpass = os.getenv("NEO4JPASS")

neo4juser = 'neo4j'
if os.getenv("NEO4JUSER"):
    neo4juser = os.getenv("NEO4JUSER")

neo4jurl = "http://neo4j.b12x.org:80"
if os.getenv("NEO4JURL"):
    neo4jurl = os.getenv("NEO4JURL")

biosqlpass = "my-secret-pw"
if os.getenv("BIOSQLPASS"):
    biosqlpass = os.getenv("BIOSQLPASS")

biosqluser = 'root'
if os.getenv("BIOSQLUSER"):
    biosqluser = os.getenv("BIOSQLUSER")

biosqlhost = "localhost"
if os.getenv("BIOSQLHOST"):
    biosqlhost = os.getenv("BIOSQLHOST")

biosqldb = "bioseqdb"
if os.getenv("BIOSQLDB"):
    biosqldb = os.getenv("BIOSQLDB")

biosqlport = 3306
if os.getenv("BIOSQLPORT"):
    biosqlport = os.getenv("BIOSQLPORT")
    if type(biosqlport) == type(biosqldb):
        biosqlport = int(biosqlport)

loadfeats = True
if os.getenv("LOADFEATS"):
    loadfeats = os.getenv("LOADFEATS")
    if loadfeats == 'True':
        loadfeats = True
    elif loadfeats == 'False':
        loadfeats = False

loadunique = True
if os.getenv("LOADUNIQUE"):
    loadunique = os.getenv("LOADUNIQUE")
    if loadunique == 'True':
        loadunique = True
    elif loadunique == 'False':
        loadunique = False

cached_feats = {}
cached_pygfe = {}
cached_unique = {}


def conn():
    try:
        conn = pymysql.connect(host=biosqlhost,
                               port=biosqlport, user=biosqluser,
                               passwd=biosqlpass, db=biosqldb)
        conn.close()
        return True
    except Exception as e:
        print("Exception while checking MYSQL Connection:" + str(e))
        return False


def typeseq_get(locus, sequence, imgthla_version='3310', neo4j_url=neo4jurl, user=neo4juser, password=neo4jpass, verbose=None):  # noqa: E501
    """typeseq_get

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
    global cached_feats
    global cached_pygfe
    global cached_unique

    graph = Graph(neo4jurl, user=neo4juser, password=neo4jpass,
                  bolt=False)

    if not cached_unique and loadunique:
        tmp_df = pd.DataFrame(graph.data(all_feats()))
        tmp_df['ID'] = tmp_df.apply(lambda row: ":".join(["".join(row['DB'].split(".")),row['LOC'],str(row['RANK']),row['TERM'],row['SEQ']]),axis=1)
        cached_unique = tmp_df['ID'].reset_index().set_index("ID").to_dict()['index']

    if not cached_feats and loadfeats:
        gfe = GFE(load_features=True,
                  loci=[locus],
                  verbose=True)
        cached_feats = gfe.all_feats

    if imgthla_version in cached_pygfe:
        pygfe = cached_pygfe[imgthla_version]
    else:
        pygfe = pyGFE(graph=graph,
                      load_features=False,
                      cached_features=cached_feats,
                      features=cached_unique,
                      load_gfe2hla=True,
                      load_seq2hla=True,
                      load_gfe2feat=True,
                      verbose=True,
                      verbosity=2)
        cached_pygfe.update({imgthla_version: pygfe})

    if conn():
        server = BioSeqDatabase.open_database(driver="pymysql",
                                              user=biosqluser,
                                              passwd=biosqlpass,
                                              host=biosqlhost,
                                              db=biosqldb,
                                              port=biosqlport)
        seqann = BioSeqAnn(server=server,
                           dbversion="".join(imgthla_version.split(".")),
                           align=False, verbose=True)
    else:
        seqann = BioSeqAnn(dbversion="".join(imgthla_version.split(".")),
                           align=False, verbose=True)

    pygfe.seqann = seqann
    typing = pygfe.type_from_seq(locus, sequence, "".join(imgthla_version.split(".")))
    typing.pygfe_version = "0.0.20"
    typing.gfedb_version = "2.0.0"
    return typing


def typealign_get(locus, sequence, imgthla_version='3310', neo4j_url=neo4jurl, user=neo4juser, password=neo4jpass, verbose=None):  # noqa: E501
    """typeseq_get

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
    graph = Graph(neo4jurl, user=neo4juser, password=neo4jpass,
                  bolt=False)

    if imgthla_version in cached_pygfe:
        pygfe = cached_pygfe[imgthla_version]
    else:
        pygfe = pyGFE(graph=graph,
                      load_features=False,
                      load_gfe2hla=True,
                      load_seq2hla=True,
                      load_gfe2feat=True,
                      verbose=True,
                      verbosity=2)
        cached_pygfe.update({imgthla_version: pygfe})

    if conn():
        server = BioSeqDatabase.open_database(driver="pymysql",
                                              user=biosqluser,
                                              passwd=biosqlpass,
                                              host=biosqlhost,
                                              db=biosqldb,
                                              port=biosqlport)
        seqann = BioSeqAnn(server=server,
                           dbversion="".join(imgthla_version.split(".")),
                           align=False, verbose=True)
    else:
        seqann = BioSeqAnn(dbversion="".join(imgthla_version.split(".")),
                           align=False, verbose=True)

    pygfe.seqann = seqann
    typing = pygfe.type_from_seq(locus, sequence, "".join(imgthla_version.split(".")))
    typing.pygfe_version = "0.0.20"
    typing.gfedb_version = "2.0.0"
    return typing

