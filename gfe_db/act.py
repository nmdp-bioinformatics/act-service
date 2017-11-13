'''
Created on Feb 8, 2017

@author: mhalagan
'''

from pygfe import pyGFE
from seqann import BioSeqAnn
from pygfe.feature_client.models.feature import Feature

from gfe_db.cypher import sequence_search
from gfe_db.cypher import gfe_search
from gfe_db.cypher import similar_gfe_classI
from gfe_db.cypher import similar_gfe_classII
from gfe_db.cypher import gfe_hla
from gfe_db.cypher import hla_Ggroups
from gfe_db.cypher import hla_gfe
from gfe_db.cypher import hla_ars
from gfe_db.cypher import gfe_ars
from gfe_db.cypher import get_sequence
from gfe_db.cypher import similar_kir
from gfe_db.cypher import get_features
from gfe_db.cypher import ref_query
from gfe_db.cypher import search_hla_features

from gfe_db.cypher import hla_alleleid
from gfe_db.cypher import gfe_alleleid
from gfe_db.cypher import fullseqid
from gfe_db.cypher import seqid
from gfe_db.cypher import search_feature
from gfe_db.cypher import persisted_query

from gfe_db.cypher import groups_classI
from gfe_db.cypher import groups_classII

from swagger_server.models.error import Error
from swagger_server.models.feature import Feature
from swagger_server.models.typing import Typing
from swagger_server.models.gfe_call import GfeCall
from swagger_server.models.gfe_typing import GfeTyping
from swagger_server.models.allele_call import AlleleCall
from swagger_server.models.feature_call import FeatureCall
from swagger_server.models.typing_status import TypingStatus
from swagger_server.models.ars_call import ArsCall
from swagger_server.models.persisted import Persisted
from swagger_server.models.persisted_data import PersistedData

from py2neo import Node, Relationship
import pandas as pa
import os
import glob
import re
import json

import sys
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from BioSQL import BioSeqDatabase
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
from Bio.Alphabet import IUPAC
import pymysql

flatten = lambda l: [item for sublist in l for item in sublist]
is_gfe = lambda x: True if re.search("\d+-\d+-\d+", x) else False
is_kir = lambda x: True if re.search("KIR", x) else False
is_classII = lambda x: True if re.search("HLA-D", x) else False
is_classI = lambda x: True if re.search("HLA-\Dw", x) else False
lc = lambda x: x.lower() if not re.search("UTR", x) else x.lower().replace("utr", "UTR")

biosqlpass = ''
if os.getenv("BIOSQLPASS"):
    biosqlpass = os.getenv("BIOSQLPASS")

biosqluser = ''
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
        print(biosqlpass, biosqluser, biosqlhost, biosqldb, biosqlport, sep="\t")
        conn = pymysql.connect(host=biosqlhost,
                               port=biosqlport, user=biosqluser,
                               passwd=biosqlpass, db=biosqldb)
        conn.close()
        return True
    except Exception as e:
        print("Exception while checking MYSQL Connection:" + str(e))
        return False


class Act(object):
    '''
    classdocs
    '''

    def __init__(self, graph, hostname, user, persist):
        '''
        Constructor
        '''
        seqann = None
        if conn():
            server = BioSeqDatabase.open_database(driver="pymysql",
                                                  user=biosqluser,
                                                  passwd=biosqlpass,
                                                  host=biosqlhost,
                                                  db=biosqldb)
            print("Server found!", file=sys.stderr)
            seqann = BioSeqAnn(server=server)
        else:
            print("No Server found!", file=sys.stderr)
            seqann = BioSeqAnn()

        gfemaker = pyGFE()
        self.gfe = gfemaker
        self.seqann = seqann
        self.user = user
        self.persist = persist
        self.graph = graph
        self.version = '0.0.5'
        api_client = ApiClient(host=hostname)
        self.api = swagger_client.DefaultApi(api_client=api_client)
        structure_dir = os.path.dirname(__file__)
        struture_files = glob.glob(structure_dir + '/data/*.structure')
        self.structures = {}
        for inputfile in struture_files:
            file_path = inputfile.split("/")
            locus = file_path[len(file_path)-1].split(".")[0]
            inputdata = (l.strip() for l in open(inputfile, "r") if l.strip())
            features = []
            for line in inputdata:
                line.rstrip()
                line.strip('\n')
                [feature, rank] = line.split("\t")
                features.append("-".join([feature, rank]))
                if is_kir(locus):
                    self.structures.update({locus: features})
                else:
                    self.structures.update({"HLA-" + locus: features})
        ref_file = structure_dir + '/data/reference_alleles.txt'
        refdata = (l.strip() for l in open(ref_file, "r") if l.strip())

        self.ref_alleles = {}
        for line in refdata:
            line.rstrip()
            line.strip('\n')
            [locus, allele] = line.split("*")
            if not locus in self.ref_alleles:
                self.ref_alleles.update({locus: [line]})
            else:
                self.ref_alleles[locus].append(line)

        max_hlaid = pa.DataFrame(self.graph.data('MATCH(hla:IMGT) RETURN max(hla.alleleId) AS ID'))
        max_gfeid = pa.DataFrame(self.graph.data('MATCH(gfe:GFE) RETURN max(gfe.alleleId) AS ID'))
        max_fullseqid = pa.DataFrame(self.graph.data('MATCH(seq:SEQUENCE) RETURN max(seq.sequenceId) AS ID'))
        max_sequenceid = pa.DataFrame(self.graph.data('MATCH(feat:FEATURE) RETURN max(feat.sequenceId) AS ID'))

        self.max_alleleid = max(max_hlaid['ID'][0], max_gfeid['ID'][0])
        self.max_seqid = max(max_fullseqid['ID'][0], max_sequenceid['ID'][0])

        self.ref_gfe = {}
        for loc in self.ref_alleles:
            ref_q = ref_query(self.ref_alleles[loc])
            self.ref_gfe.update({loc: pa.DataFrame(self.graph.data(ref_q))})

        self.admin = ''
        if os.getenv("NEO4JADMIN"):
            self.admin = os.getenv("NEO4JADMIN")

    def get_seqid(self, sequence, seqtype, rank):

        if seqtype == "SEQUENCE":
            fullseq_q = fullseqid(sequence)
            fullid = pa.DataFrame(self.graph.data(fullseq_q))
            if not fullid.empty:
                return fullid['ID'][0]
            else:
                self.max_seqid += 1
                return self.max_seqid
        else:
            seq_q = seqid(sequence, seqtype, rank)
            seq_id = pa.DataFrame(self.graph.data(seq_q))
            if not seq_id.empty:
                return seq_id['ID'][0]
            else:
                self.max_seqid += 1
                return self.max_seqid

    def get_alleleid(self, typing):

        if is_gfe(typing):
            gfe_q = gfe_alleleid(typing)
            gfeid = pa.DataFrame(self.graph.data(gfe_q))
            if not gfeid.empty:
                return gfeid['ID'][0]
            else:
                self.max_alleleid += 1
                return self.max_alleleid
        else:
            hla_q = hla_alleleid(typing)
            hlaid = pa.DataFrame(self.graph.data(hla_q))
            if not hlaid.empty:
                return hlaid['ID'][0]
            else:
                self.max_alleleid += 1
                return self.max_alleleid

    def persist_typing(self, typing):

        [loc, accesions] = typing.gfe.split("w")

        gfe_alleleid = self.get_alleleid(typing.gfe)
        gfe_node = Node("GFE", name=str(typing.gfe),
                        imgtdb=str(typing.imgtdb_version),
                        locus=str(loc),
                        alleleId=int(gfe_alleleid),
                        status=str("persisted"))

        full_seq = typing.full_gene.sequence
        full_seqid = self.get_seqid(full_seq, "SEQUENCE", 0)
        seq_node = Node("SEQUENCE", name="Sequence", rank=int(0),
                        sequenceId=int(full_seqid),
                        sequence=str(full_seq),
                        length=len(full_seq),
                        nuc=[s for s in list(full_seq)],
                        status=str("persisted"))

        hla_nodes = []
        relationships = []
        for hla_typing in typing.typing:
            hla_alleleid = self.get_alleleid(hla_typing.hla)
            hla_node = Node("IMGT", name=hla_typing.hla,
                            imgtdb=typing.imgtdb_version, locus=loc,
                            alleleId=int(hla_alleleid))
            gfe_rel = Relationship(hla_node, "HAS_GFE", gfe_node,
                                   status="persisted")
            seqgfe_rel = Relationship(hla_node, "HAS_FEATURE", seq_node,
                                      status="persisted")
            seqhla_rel = Relationship(gfe_node, "HAS_FEATURE", seq_node,
                                      status="persisted")
            hla_nodes.append(hla_node)
            relationships.append(gfe_rel)
            relationships.append(seqgfe_rel)
            relationships.append(seqhla_rel)

        for feat in typing.typing_status.novel_features:
            feat_seqid = self.get_seqid(feat.sequence,
                                        feat.term.upper(), feat.rank)
            feature_node = Node("FEATURE", name=str(feat.term.upper()),
                                rank=str(feat.rank),
                                sequenceId=int(feat_seqid),
                                sequence=str(feat.sequence),
                                length=len(feat.sequence),
                                nuc=[s for s in list(feat.sequence)],
                                status=str("persisted"))
            gfefeat_rel = Relationship(gfe_node, "HAS_FEATURE",
                                       feature_node, accession=str(feat.accession),
                                       status="persisted")
            relationships.append(gfefeat_rel)
            for hla_n in hla_nodes:
                hlafeat_rel = Relationship(hla_n, "HAS_FEATURE",
                                           feature_node,
                                           accession=str(feat.accession),
                                           status="persisted")
                relationships.append(hlafeat_rel)

        for feat in typing.features:
            if not feat in typing.typing_status.novel_features:
                feat_seqid = self.get_seqid(feat.sequence,
                                            feat.term.upper(), feat.rank)
                feature_node = Node("FEATURE", name=str(feat.term.upper()),
                                    rank=str(feat.rank),
                                    sequenceId=int(feat_seqid),
                                    sequence=str(feat.sequence),
                                    length=str(len(feat.sequence)),
                                    nuc=[s for s in list(feat.sequence)])
                gfefeat_rel = Relationship(gfe_node, "HAS_FEATURE",
                                           feature_node, accession=str(feat.accession),
                                           status="persisted")
                relationships.append(gfefeat_rel)
                for hla_n in hla_nodes:
                    hlafeat_rel = Relationship(hla_n, "HAS_FEATURE",
                                               feature_node,
                                               accession=str(feat.accession),
                                               status="persisted")
                    relationships.append(hlafeat_rel)

        tx = self.graph.begin()
        for rel in relationships:
            tx.merge(rel)
        tx.commit()

    def unique_features(self, features):

        unique = []
        for feat in features:
            feat_q = search_feature(feat.term, feat.rank, feat.sequence)
            seq_features = pa.DataFrame(self.graph.data(feat_q))
            if seq_features.empty:
                unique.append(feat)
        return unique

    def type_hla(self, locus, sequence, type_gfe):

        ac_object = AlleleCall()
        ac_object.act_version = self.version
        ac_object.gfedb_version = '0.0.2'
        ac_object.typing_status = TypingStatus()

        if type_gfe:
            ac_object.gfe = type_gfe
            seq_features = pa.DataFrame(self.graph.data(get_features(type_gfe)))
            if not seq_features.empty:
                ac_object.typing_status.status = "documented"
                features = list()
                for i in range(0, len(seq_features['term'])):
                    feature = Feature(accession=seq_features['accession'][i], rank=seq_features['rank'][i], sequence=seq_features['sequence'][i], term=lc(seq_features['term'][i]))
                    features.append(feature)
                seq_o = self.gfe_sequence(locus, type_gfe)
                ac_object.full_gene = Feature(rank="1", sequence=seq_o.sequence, term="gene")
                ac_object.features = features
                ac_object.ihiw_ref = self.get_ref_allele(locus, type_gfe, ac_object.features)
                related_gfe = self.gfe_lookup(type_gfe, ac_object.features)
                ac_object.typing = related_gfe
            else:
                seq_o = self.gfe_sequence(locus, type_gfe)
                ac_object.full_gene = Feature(rank="1", sequence=seq_o.sequence, term="gene")
                ac_object.features = [Feature(accession=f.accession, rank=f.rank, sequence=f.sequence, term=lc(f.term)) for f in seq_o.structure]
                ac_object.typing_status.novel_features = self.unique_features(ac_object.features)
                ac_object.ihiw_ref = self.get_ref_allele(locus, type_gfe, ac_object.features)
                ac_object.typing = self.find_similar(ac_object.gfe, ac_object.features)

                if(len(ac_object.typing_status.novel_features) != 0):
                    ac_object.typing_status.status = "novel"
                else:
                    ac_object.typing_status.status = "novel_combination"

                if self.admin == self.user and self.persist:
                    self.persist_typing(ac_object)
            return ac_object
        else:
            sequence = sequence.upper()
            sequence_typing = self.sequence_lookup(locus, sequence)
            if sequence_typing:
                ac_object.typing_status.status = "documented"
                ac_object.typing = [sequence_typing[0]]
                ac_object.gfe = sequence_typing[1]
                ac_object.features = sequence_typing[2]
                ac_object.full_gene = Feature(rank="1", sequence=sequence, term="gene")
                ac_object.ihiw_ref = self.get_ref_allele(locus, ac_object.gfe, ac_object.features)
                return ac_object
            else:

                gfe_o = self.gfe_create(locus, sequence)
                # if hasattr(gfe_o, 'body'):
                #     error = json.loads(gfe_o.body)
                #     error_o = Error()
                #     for k in error:
                #         setattr(error_o, k.lower(), error[k])
                #     return error_o
                ac_object.gfe = gfe_o['gfe']
                #ac_object.full_gene = Feature(accession=gfe_o['fullgene'].accession, rank=gfe_o.fullgene.rank, sequence=gfe_o.fullgene.sequence, term=gfe_o.fullgene.term)
                #ac_object.gfe_version = gfe_o['version']
                #ac_object.ihiw_ref = self.get_ref_allele(locus, ac_object.gfe, gfe_o['structure'])
                ac_object.features = [Feature(accession=f.accession, rank=f.rank, sequence=f.sequence, term=f.term) for f in gfe_o['structure']]
                ac_object.typing_status.novel_features = self.unique_features(ac_object.features)
                related_gfe = self.gfe_lookup(ac_object.gfe, ac_object.features)

                # ac_object.gfe = gfe_o.gfe
                # ac_object.full_gene = Feature(accession=gfe_o.fullgene.accession, rank=gfe_o.fullgene.rank, sequence=gfe_o.fullgene.sequence, term=gfe_o.fullgene.term)
                # ac_object.gfe_version = gfe_o.version
                # ac_object.ihiw_ref = self.get_ref_allele(locus, ac_object.gfe, gfe_o.structure)
                ac_object.features = [Feature(accession=f.accession, rank=f.rank, sequence=f.sequence, term=f.term) for f in gfe_o['structure']]
                ac_object.typing_status.novel_features = self.unique_features(ac_object.features)
                related_gfe = self.gfe_lookup(ac_object.gfe, ac_object.features)

                if(len(ac_object.typing_status.novel_features) != 0):
                    ac_object.typing_status.status = "novel"
                else:
                    ac_object.typing_status.status = "novel_combination"

                if related_gfe:
                    ac_object.typing = related_gfe
                else:
                    ac_object.typing = self.find_similar(ac_object.gfe, ac_object.features)

                if self.admin == self.user and self.persist:
                    self.persist_typing(ac_object)

                return ac_object

    def sequence_lookup(self, locus, sequence):
        """
        Looks up sequence from

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        lookup_query = sequence_search(locus, sequence)
        sequence_data = pa.DataFrame(self.graph.data(lookup_query))
        if not sequence_data.empty:
            features = list()
            gfe = list(set([x for x in sequence_data["GFE"]]))
            hla = list(set([x for x in sequence_data["HLA"]]))
            seq_features = pa.DataFrame(self.graph.data(get_features(gfe[0])))
            for i in range(0, len(seq_features['term'])):
                feature = Feature(accession=seq_features['accession'][i], rank=seq_features['rank'][i], sequence=seq_features['sequence'][i], term=lc(seq_features['term'][i]))
                features.append(feature)

            typing = Typing(hla=hla[0], related_gfe=[GfeTyping(gfe=gfe[0], shares=features, features_shared=len(features))])
            return [typing, gfe[0], features]
        else:
            return

    def gfe_create(self, locus, sequence):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        seq_rec = SeqRecord(Seq(sequence, IUPAC.unambiguous_dna), id="GFE")
        annotation = self.seqann.annotate(seq_rec, locus)
        features, gfe = self.gfe.get_gfe(annotation, locus)
        return {'gfe': gfe, 'structure': features}
        # try:
        #     return self.api.gfe_post(locus=locus, sequence=sequence, verbose=1)
        # except ApiException as e:
        #     print("Exception when calling DefaultApi->gfe_post: %s\n" % e)
        #     return e

    def gfe_sequence(self, locus, gfe):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        try:
            return self.api.sequence_post(locus=locus, gfe=gfe, verbose=1)
        except ApiException as e:
            print("Exception when calling DefaultApi->sequence_post: %s\n" % e)
            return e

    def gfe_lookup(self, gfe, features):

        gfe_data = pa.DataFrame(self.graph.data(gfe_search(gfe)))
        if not gfe_data.empty:
            typing_list = list()
            for hla in gfe_data["HLA"]:
                typing_list.append(Typing(hla=hla, related_gfe=[GfeTyping(gfe=gfe, shares=features, features_shared=len(features))]))
            return typing_list
        else:
            return

    def get_class(self, gfe):

        [loc, accessions] = gfe.split("w")
        if loc in self.loci:
            return self.loci[loc]
        return

    def find_similar(self, gfe, features):

        if is_classI(gfe):
            return self.find_gfe_classI(gfe, features)
        elif is_classII(gfe):
            return self.find_gfe_classII(gfe, features)
        elif is_kir(gfe):
            return self.find_gfe_kir(gfe, features)
        else:
            return

    def find_gfe_classI(self, gfe, features):

        gfe_dict = self.breakup_gfe(gfe)
        [locus, feature_accessions] = gfe.split("w")
        cypher = similar_gfe_classI(gfe, gfe_dict["exon-2"], gfe_dict["exon-3"])
        similar_data = pa.DataFrame(self.graph.data(cypher))
        if not similar_data.empty:
            gfe_similarity = {}
            for i in range(0, len(similar_data['GFE1'])):
                sim1 = self.calcDiff(similar_data['GFE1'][i], gfe)
                sim2 = self.calcDiff(similar_data['GFE2'][i], gfe)
                gfe_similarity.update({similar_data['GFE1'][i]: sim1})
                gfe_similarity.update({similar_data['GFE2'][i]: sim2})

            found_hla = {}
            max_val = max(gfe_similarity.values())
            for gfes in gfe_similarity.keys():
                if gfe_similarity[gfes] == max_val:
                    hlas = self.gfe2hla(gfes)
                    for hla in hlas:
                        if hla not in found_hla:
                            matched_features = self.matching_features(gfe, gfes, self.map_structures(features))
                            hla_typing = Typing(hla=hla, related_gfe=[GfeTyping(gfe=gfes, shares=matched_features, features_shared=len(matched_features))])
                            found_hla.update({hla: hla_typing})
                        else:
                            hla_typing = found_hla[hla]
                            matched_features = self.matching_features(gfe, gfes, self.map_structures(features))
                            hla_typing.related_gfe.append(GfeTyping(gfe=gfes, shares=matched_features, features_shared=len(matched_features)))
                            found_hla.update({hla: hla_typing})

            typing_list = list(found_hla.values())
            return typing_list
        else:
            return list()

    def find_gfe_classII(self, gfe, features):

        gfe_dict = self.breakup_gfe(gfe)
        [locus, feature_accessions] = gfe.split("w")
        cypher = similar_gfe_classII(gfe, gfe_dict["exon-2"])
        similar_data = pa.DataFrame(self.graph.data(cypher))
        if not similar_data.empty:
            gfe_similarity = {}
            for i in range(0, len(similar_data['GFE1'])):
                sim1 = self.calcDiff(similar_data['GFE1'][i], gfe)
                #sim2 = self.calcDiff(similar_data['GFE2'][i], gfe)
                gfe_similarity.update({similar_data['GFE1'][i]: sim1})
                #gfe_similarity.update({similar_data['GFE2'][i]: sim2})

            found_hla = {}
            max_val = max(gfe_similarity.values())
            for gfes in gfe_similarity.keys():
                if gfe_similarity[gfes] == max_val:
                    hlas = self.gfe2hla(gfes)
                    for hla in hlas:
                        if hla not in found_hla:
                            matched_features = self.matching_features(gfe, gfes, self.map_structures(features))
                            hla_typing = Typing(hla=hla, related_gfe=[GfeTyping(gfe=gfes, shares=matched_features, features_shared=len(matched_features))])
                            found_hla.update({hla: hla_typing})
                        else:
                            hla_typing = found_hla[hla]
                            matched_features = self.matching_features(gfe, gfes, self.map_structures(features))
                            hla_typing.related_gfe.append(GfeTyping(gfe=gfes, shares=matched_features, features_shared=len(matched_features)))
                            found_hla.update({hla: hla_typing})

            typing_list = list(found_hla.values())
            return typing_list

        else:
            return list()

    def find_gfe_kir(self, gfe, features):

        [locus, feature_accessions] = gfe.split("w")

        cypher = similar_kir(locus)
        similar_data = pa.DataFrame(self.graph.data(cypher))
        if not similar_data.empty:
            gfe_similarity = {}
            for i in range(0, len(similar_data['GFE'])):
                sim1 = self.calcDiff(similar_data['GFE'][i], gfe)
                gfe_similarity.update({similar_data['GFE'][i]: sim1})

            found_hla = {}
            max_val = max(gfe_similarity.values())
            for gfes in gfe_similarity.keys():
                if gfe_similarity[gfes] == max_val:
                    hlas = self.gfe2hla(gfes)
                    for hla in hlas:
                        if hla not in found_hla:
                            matched_features = self.matching_features(gfe, gfes, self.map_structures(features))
                            hla_typing = Typing(hla=hla, related_gfe=[GfeTyping(gfe=gfes, shares=matched_features, features_shared=len(matched_features))])
                            found_hla.update({hla: hla_typing})
                        else:
                            hla_typing = found_hla[hla]
                            matched_features = self.matching_features(gfe, gfes, self.map_structures(features))
                            hla_typing.related_gfe.append(GfeTyping(gfe=gfes, shares=matched_features, features_shared=len(matched_features)))
                            found_hla.update({hla: hla_typing})

            typing_list = list(found_hla.values())
            return typing_list

        else:
            return list()

    def map_structures(self, gfe_structs):

        feat_map = {}
        for feat in gfe_structs:
            feat_name = "-".join([feat.term, str(feat.rank)])
            feat_map.update({feat_name: feat.sequence})

        return feat_map

    def matching_features(self, gfe1, gfe2, structures):

        gfe_parts1 = self.breakup_gfe(gfe1)
        gfe_parts2 = self.breakup_gfe(gfe2)
        feat_list = list()
        for feat in gfe_parts1:
            if feat in gfe_parts2:
                if gfe_parts1[feat] == gfe_parts2[feat]:
                    feat_term, feat_rank = feat.split('-')
                    shared_feat = Feature(term=feat_term, rank=feat_rank, sequence=structures[feat], accession=gfe_parts1[feat])
                    feat_list.append(shared_feat)

        return(feat_list)

    def breakup_gfe(self, gfe):
        [locus, feature_accessions] = gfe.split("w")
        accessions = feature_accessions.split("-")

        if locus == "HLA-DQB1":
            if len(accessions) < len(self.structures[locus]):
                i = 0
                features = {}
                old_dq = ["intro_1","exon_2","intron_2","exon_3","intro_3"]
                for feature_rank in old_dq:
                    accession = accessions[i]
                    features.update({feature_rank: accession})
                    i += 1

                return(features)
            else:
                i = 0
                features = {}
                for feature_rank in self.structures[locus]:
                    accession = accessions[i]
                    features.update({feature_rank: accession})
                    i += 1

                return(features)
        else:
            i = 0
            features = {}
            for feature_rank in self.structures[locus]:
                accession = accessions[i]
                features.update({feature_rank: accession})
                i += 1

            return(features)

    def gfe2hla(self, gfe):

        cypher = gfe_hla(gfe)
        hlacypher = pa.DataFrame(self.graph.data(cypher))
        hla_alleles = set()

        if not len(hlacypher["HLA"][0]) == 0:
            hla_alleles = flatten(hlacypher["HLA"])
            return hla_alleles
        else:
            return

    def g_group(self, hla):

        g_cypher = hla_Ggroups(hla)
        g_data = pa.DataFrame(self.graph.data(g_cypher))
        g_groups = {}
        if not g_data.empty:
            for g in g_data["G_GROUP"]:
                if g not in g_groups:
                    g_groups.update({g: "G"})
        big_g = [g for g in g_groups]
        if len(big_g) > 1 or len(big_g) == 0:
            return
        else:
            return big_g[0]

    def calcDiff(self, gfe1, gfe2):

        count = 0
        gfe1_parts = gfe1.split("-")
        gfe2_parts = gfe2.split("-")
        for i in range(0, len(gfe1_parts)):
            if gfe1_parts[i] == gfe2_parts[i]:
                count += 1

        return count

    def get_gfe_call(self, hla):
        gfe = self.get_gfe(hla)
        gfe_c = GfeCall(hla=gfe, act_version=self.version,
                        gfedb_version='0.0.2')
        return gfe_c

    def get_gfe(self, hla):
        gfe_query = hla_gfe(hla)
        gfe_data = pa.DataFrame(self.graph.data(gfe_query))
        if not gfe_data.empty:
            gfe = [x for x in gfe_data["GFE"]]
            return gfe
        else:
            return

    def get_hla(self, gfe):
        hla_query = gfe_hla(gfe)
        hla_data = pa.DataFrame(self.graph.data(hla_query))
        if not len(hla_data["HLA"][0]) == 0:
            hla_alleles = flatten(hla_data["HLA"])
            return hla_alleles
        else:
            return

    def sequence(self, seq_type, allele):
        seq_query = get_sequence(seq_type, allele)
        seq_data = pa.DataFrame(self.graph.data(seq_query))
        if not seq_data.empty:
            seq = [x for x in seq_data["SEQ"]]
            return seq
        else:
            return

    def ars_redux(self, group, typing):
        # Have share_ars group

        ars_call = ArsCall()
        ars_call.act_version = self.version
        ars_call.gfedb_version = '0.0.2'
        if is_gfe(typing):
            gfe_query = gfe_ars(group, typing)
            gfe_data = pa.DataFrame(self.graph.data(gfe_query))
            if not gfe_data.empty:
                found_hla = {}
                for i in range(0, len(gfe_data['HLA'])):
                    hla = gfe_data['HLA'][i]
                    hla_typing = Typing(hla=hla)
                    found_hla.update({hla: hla_typing})

                ars_call.allele = typing
                ars_call.group_type = group
                ars_call.group = gfe_data["ARS"][0]
                ars_call.share_allele = list(found_hla.values())
                return ars_call
            else:
                #groups_classI(locus, group, exon2, exon3):
                [loc, typ] = typing.split("w")
                gfe_parts = self.breakup_gfe(typing)
                group_q = groups_classI(loc, group, gfe_parts['exon-2'], gfe_parts['exon-3'])
                ars_data = pa.DataFrame(self.graph.data(group_q))
                if not ars_data.empty:
                    found_hla = {}
                    hla_alleles = flatten(ars_data["HLA"])
                    ars_alleles = flatten(ars_data["ARS"])
                    for i in range(0, len(hla_alleles)):
                        hla = hla_alleles[i]
                        hla_typing = Typing(hla=hla)
                        found_hla.update({hla: hla_typing})

                    ars_call.allele = typing
                    ars_call.group_type = group
                    ars_call.group = ars_alleles[0]
                    ars_call.share_allele = list(found_hla.values())
                    return ars_call
                else:
                    return ars_call
        else:
            hla_query = hla_ars(group, typing)
            hla_data = pa.DataFrame(self.graph.data(hla_query))
            if not hla_data.empty:
                found_hla = {}
                for i in range(0, len(hla_data['HLA'])):
                    hla = hla_data['HLA'][i]
                    gfe = hla_data['GFE'][i]
                    if hla not in found_hla:
                        typ = Typing(hla=hla, related_gfe=[GfeTyping(gfe=gfe)])
                        found_hla.update({hla: typ})
                    else:
                        typ = found_hla[hla]
                        typ.related_gfe.append(GfeTyping(gfe=gfe))
                        found_hla.update({hla: typ})

                ars_call.allele = typing
                ars_call.group_type = group
                ars_call.group = hla_data["ARS"][0]
                ars_call.share_allele = list(found_hla.values())
                return ars_call
            else:
                return ars_call

    def typing_to_bioseq(self, typing):

        # TODO: Add more annotation and qualifiers
        # TODO: Use the full sequence accession as the ID
        sequence = typing.full_gene.sequence
        seqid = 'GFE'
        comments = []
        description = "GFE " + typing.gfe
        comments.append("Typing Status: " + typing.typing_status.status)
        if hasattr(typing.typing_status, 'novel_features'):
            uniq = " ".join(["-".join([feat.term, str(feat.rank)]) for feat in typing.typing_status.novel_features])
            comments.append("Novel features: " + uniq)
        
        allele_typed = "/".join([typ.hla for typ in typing.typing])
        comments.append("Allele Call: " + allele_typed)
        comments.append("IHIW Reference: " + typing.ihiw_ref[0].hla)
        comments.append("")
        comments.append("Typed with ACT Service " + self.version)
        if hasattr(typing, 'gfe_version'):
            comments.append("Annotated with GFE Service " + typing.gfe_version)

        if hasattr(typing.full_gene, 'accession'):
            seqid = 'GFEw' + str(typing.full_gene.accession)
        seqrecord = SeqRecord(Seq(sequence, IUPAC.unambiguous_dna), id=seqid, description=description)
        source_feature = SeqFeature(FeatureLocation(0, len(str(seqrecord.seq))), type="source", strand=1)
        seqrecord.annotations["sequence_version"] = 1
        seqrecord.annotations["molecule_type"] = "DNA"
        seqrecord.annotations["data_file_division"] = "HUM"
        if hasattr(typing.full_gene, 'accession'):
            seqrecord.annotations["accessions"] = [seqid]

        seqrecord.annotations["comment"] = comments
        seqrecord.features.append(source_feature)
        seqrecord.features[0].qualifiers = OrderedDict([('organism', ['Homo sapiens']), ('mol_type', ['genomic DNA']), ('db_xref', ['taxon:9606'])])

        for feat in typing.features:
            start_pos = sequence.find(feat.sequence)
            if len(seqrecord.features) == 1:
                start_pos = 0

            feat_type = feat.term
            if feat_type == "five_prime_UTR" or feat_type == "three_prime_UTR":
                feat_type = "UTR"
            end_pos = start_pos + len(feat.sequence)
            seq_feature = SeqFeature(FeatureLocation(start_pos, end_pos), type=feat_type, strand=1)
            
            if feat.term == 'exon' or feat.term == 'intron':
                seq_feature.qualifiers = OrderedDict([('number', [feat.rank])])
            seqrecord.features.append(seq_feature)

        return seqrecord

    def get_ref_allele(self, locus, gfe, features):

        hla_sim = {}
        gfe_mapped = {}
        for i in range(0, len(self.ref_gfe[locus]['HLA'])):
            sim1 = self.calcDiff(self.ref_gfe[locus]['GFE'][i], gfe)
            hla_sim.update({self.ref_gfe[locus]['HLA'][i]: sim1})
            if not self.ref_gfe[locus]['HLA'][i] in gfe_mapped:
                gfe_mapped.update({self.ref_gfe[locus]['HLA'][i]: [self.ref_gfe[locus]['GFE'][i]]})
            else:
                gfe_mapped[self.ref_gfe[locus]['HLA'][i]].append(self.ref_gfe[locus]['GFE'][i])

        found_hla = {}
        max_val = max(hla_sim.values())
        for hla in hla_sim.keys():
            if hla_sim[hla] == max_val:
                gfes_sim = gfe_mapped[hla]
                for gfes in gfes_sim:
                    if hla not in found_hla:
                        matched_features = self.matching_features(gfe, gfes, self.map_structures(features))
                        hla_typing = Typing(hla=hla, related_gfe=[GfeTyping(gfe=gfes, shares=matched_features, features_shared=len(matched_features))])
                        found_hla.update({hla: hla_typing})
                    else:
                        hla_typing = found_hla[hla]
                        matched_features = self.matching_features(gfe, gfes, self.map_structures(features))
                        hla_typing.related_gfe.append(GfeTyping(gfe=gfes, shares=matched_features, features_shared=len(matched_features)))
                        found_hla.update({hla: hla_typing})

        return list(found_hla.values())

    def get_gfe_features(self, hlas, features):

        gfe_features = {}
        for hla in hlas:
            gfes = self.get_gfe(hla)
            for gfe in gfes:
                gfe_parts = self.breakup_gfe(gfe)
                for feat in features:
                    f = lc(feat)
                    if not f in gfe_features:
                        gfe_features.update({f: [gfe_parts[f]]})
                    else:
                        fa = gfe_features[f]
                        if not gfe_parts[f] in fa:
                            fa.append(gfe_parts[f])
                            gfe_features[f] = fa
        return gfe_features

    def search_features(self, hlas, features):

        [locus, allele] = hlas[0].split("*")
        search_feats = self.get_gfe_features(hlas, features)
        feat_query = search_hla_features("HLA-A", search_feats)
        feat_data = pa.DataFrame(self.graph.data(feat_query))
        fc = FeatureCall(alleles=hlas, features_searched=features,
                         act_version=self.version, gfedb_version='0.0.2')
        found_hla = {}
        if not feat_data.empty:
            for i in range(0, len(feat_data['HLA'])):
                hla = feat_data['HLA'][i]
                gfe = feat_data['GFE'][i]
                if hla not in found_hla:
                    hla_typing = Typing(hla=hla,
                                        related_gfe=[GfeTyping(gfe=gfe)])
                    found_hla.update({hla: hla_typing})
                else:
                    hla_typing = found_hla[hla]
                    hla_typing.related_gfe.append(GfeTyping(gfe=gfe))
                    found_hla.update({hla: hla_typing})
            fc.matched = list(found_hla.values())
            return fc
        else:
            return fc

    def get_persisted(self):

        persisted = Persisted(act_version=self.version, gfedb_version='0.0.2')
        per_data = pa.DataFrame(self.graph.data(persisted_query()))
        persisted_a = []
        if not per_data.empty:
            for i in range(0, len(per_data['HLA'])):
                per = PersistedData(hla=per_data['HLA'][i],gfe=per_data['GFE'][i],
                                    term=per_data['TERM'][i],rank=per_data['RANK'][i],
                                    accession=per_data['ACCESSION'][i],
                                    sequence=per_data['SEQUENCE'][i])
                persisted_a.append(per)

        persisted.persisted_data = persisted_a
        return persisted

