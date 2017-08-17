'''
Created on Feb 8, 2017

@author: mhalagan
'''

from gfe_db.cypher import sequence_search
from gfe_db.cypher import gfe_search
from gfe_db.cypher import similar_gfe
from gfe_db.cypher import gfe_hla
from gfe_db.cypher import groups
from gfe_db.cypher import gfe_Ggroups
from gfe_db.cypher import hla_Ggroups
from gfe_db.cypher import hla_gfe
from gfe_db.cypher import hla_ars
from gfe_db.cypher import gfe_ars
from gfe_db.cypher import get_sequence

import pandas as pa
import swagger_client
from swagger_client.rest import ApiException
from swagger_client.api_client import ApiClient
import os
import glob
import re

flatten = lambda l: [item for sublist in l for item in sublist]
is_gfe = lambda x: True if re.search("\d+-\d+-\d+", x) else False


class Act(object):
    '''
    classdocs
    '''

    def __init__(self, graph, hostname="gfe.b12x.org"):
        '''
        Constructor
        '''
        self.graph = graph
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
                self.structures.update({"HLA-" + locus: features})

    def type_hla(self, locus, sequence):

        sequence_exists = self.sequence_lookup(locus, sequence)
        if sequence_exists:
            return sequence_exists
        else:
            gfe_new = self.gfe_create(locus, sequence)
            gfe_exists = self.gfe_lookup(gfe_new)
            if gfe_exists:
                return gfe_exists
            else:
                gfe_closest = self.find_similar_gfe(gfe_new)
                if gfe_closest:
                    return gfe_closest
                else:
                    return(gfe_new, "", "", "", 0)

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
            gfe_observed = list(set([x for x in sequence_data["GFE"]]))
            hla = list(set([x for x in sequence_data["HLA"]]))
            g_groups = []
            for gfe in gfe_observed:
                g_cypher = gfe_Ggroups(gfe)
                g_data = pa.DataFrame(self.graph.data(g_cypher))
                if not g_data.empty:
                    for g in g_data["G_GROUP"]:
                        g_groups.append(g)
            if len(g_groups) > 1 or len(g_groups) == 0:
                return(gfe, hla, gfe_observed, '', 99)
            else:
                return(gfe_observed, hla, gfe_observed, g_groups[0], 99)
        else:
            return ""

    def gfe_create(self, locus, sequence):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        try:
            r = self.api.gfe_post(locus=locus, sequence=sequence, verbose=1)
            gfe = r.gfe
            return(gfe)
        except ApiException as e:
            print("Exception when calling DefaultApi->gfe_post: %s\n" % e)
            return ""

    def gfe_lookup(self, gfe):

        lookup_query = gfe_search(gfe)
        gfe_data = pa.DataFrame(self.graph.data(lookup_query))
        if not gfe_data.empty:
            gfe_observed = [x for x in gfe_data["GFE"]]
            hla = [x for x in gfe_data["HLA"]]
            g_cypher = gfe_Ggroups(gfe)
            g_data = pa.DataFrame(self.graph.data(g_cypher))
            if not g_data.empty:
                g_groups = [g for g in g_data["G_GROUP"]]
                if len(g_groups) > 1 or len(g_groups) == 0:
                    return(gfe, hla, gfe_observed, '', 99)
                else:
                    return(gfe, hla, gfe_observed, g_groups[0], 99)
            else:
                return(gfe, hla, gfe_observed, '', 99)
        else:
            return ""

    def find_similar_gfe(self, gfe):

        gfe_dict = self.breakup_gfe(gfe)
        [locus, feature_accessions] = gfe.split("w")
        groups_cypher = groups(locus, gfe_dict["exon-2"], gfe_dict["exon-3"])
        groups_data = pa.DataFrame(self.graph.data(groups_cypher))
        if not len(groups_data["HLA"][0]) == 0:
            cypher = similar_gfe(gfe, gfe_dict["exon-2"], gfe_dict["exon-3"])
            similar_data = pa.DataFrame(self.graph.data(cypher))
            if not similar_data.empty:
                gfe_similarity = {}
                for i in range(0, len(similar_data['GFE1'])):
                    sim1 = self.calcDiff(similar_data['GFE1'][i], gfe)
                    sim2 = self.calcDiff(similar_data['GFE2'][i], gfe)
                    gfe_similarity.update({similar_data['GFE1'][i]: sim1})
                    gfe_similarity.update({similar_data['GFE2'][i]: sim2})

                gfe_dict = {}
                hla_dict = {}
                max_val = max(gfe_similarity.values())
                for gfes in gfe_similarity.keys():
                    if gfe_similarity[gfes] == max_val:
                        hlas = self.gfe2hla(gfes)
                        for hla in hlas:
                            if hla not in hla_dict:
                                hla_dict.update({hla: max_val})
                            if gfes not in gfe_dict:
                                gfe_dict.update({gfes: max_val})

                gfe_array = [g for g in gfe_dict]
                hla_array = [h for h in hla_dict]

                g_groups = {}
                for sim_gfe in gfe_array:
                    g_cypher = gfe_Ggroups(sim_gfe)
                    g_data = pa.DataFrame(self.graph.data(g_cypher))
                    if not g_data.empty:
                        for g in g_data["G_GROUP"]:
                            g_groups.update({g: sim_gfe})
                big_g = [g for g in g_groups]
                if len(big_g) > 1 or len(big_g) == 0:
                    return(gfe, hla_array, gfe_array, '', 99)
                else:
                    return(gfe, hla_array, gfe_array, big_g[0], 99)
            else:
                g_groups = {}
                hla_a = flatten(groups_data["HLA"])
                for sim_hla in hla_a:
                    g_cypher = hla_Ggroups(sim_hla)
                    g_data = pa.DataFrame(self.graph.data(g_cypher))
                    if not g_data.empty:
                        for g in g_data["G_GROUP"]:
                            if g not in g_groups:
                                g_groups.update({g: sim_hla})
                big_g = [g for g in g_groups]
                if len(big_g) > 1 or len(big_g) == 0:
                    return(gfe, hla_a, gfe, '', 99)
                else:
                    return(gfe, hla_a, gfe, big_g[0], 99)
        else:
            return ""

    def breakup_gfe(self, gfe):
        [locus, feature_accessions] = gfe.split("w")
        accessions = feature_accessions.split("-")
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
            return ""

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
            return ''
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

    def get_gfe(self, hla):
        gfe_query = hla_gfe(hla)
        gfe_data = pa.DataFrame(self.graph.data(gfe_query))
        if not gfe_data.empty:
            gfe = [x for x in gfe_data["GFE"]]
            return gfe
        else:
            return ''

    def get_hla(self, gfe):
        hla_query = gfe_hla(gfe)
        hla_data = pa.DataFrame(self.graph.data(hla_query))
        if not len(hla_data["HLA"][0]) == 0:
            hla_alleles = flatten(hla_data["HLA"])
            return hla_alleles
        else:
            return ''

    def sequence(self, seq_type, allele):
        seq_query = get_sequence(seq_type, allele)
        seq_data = pa.DataFrame(self.graph.data(seq_query))
        if not seq_data.empty:
            seq = [x for x in seq_data["SEQ"]]
            return seq
        else:
            return ''

    def ars_redux(self, group, typing):
        # Have share_ars group
        if is_gfe(typing):
            gfe_query = gfe_ars(group, typing)
            gfe_data = pa.DataFrame(self.graph.data(gfe_query))
            if not gfe_data.empty:
                gfe = [x for x in gfe_data["ARS"]]
                return gfe
            else:
                return ''
        else:
            hla_query = hla_ars(group, typing)
            hla_data = pa.DataFrame(self.graph.data(hla_query))
            if not hla_data.empty:
                hla = [x for x in hla_data["ARS"]]
                return hla
            else:
                return ''

