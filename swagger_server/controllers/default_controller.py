import connexion
import six

from pygfe.models.error import Error  # noqa: E501
from pygfe.models.seqdiff import Seqdiff  # noqa: E501
from pygfe.models.typing import Typing  # noqa: E501
from pygfe import util

from py2neo import Graph
from pygfe.pygfe import pyGFE
import os

from io import StringIO
from Bio import SeqIO

from seqann.sequence_annotation import BioSeqAnn
from BioSQL import BioSeqDatabase
import pymysql

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


def conn():
    try:
        # print(biosqlpass, biosqluser, biosqlhost, biosqldb, biosqlport, sep="\t")
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
    graph = Graph(neo4jurl, user=neo4juser, password=neo4jpass,
                  bolt=False)
    if conn():
        server = BioSeqDatabase.open_database(driver="pymysql",
                                              user=biosqluser,
                                              passwd=biosqlpass,
                                              host=biosqlhost,
                                              db=biosqldb)
        seqann = BioSeqAnn(server=server, align=True, verbose=True)
    else:
        seqann = BioSeqAnn()
    pygfe = pyGFE(graph=graph,
                  seqann=seqann,
                  verbose=True,
                  verbosity=2,
                  loci=[locus])
    typing = pygfe.type_from_seq(locus, sequence)
    typing.pygfe_version = "0.0.14"
    typing.gfedb_version = "2.0.0"
    return typing

