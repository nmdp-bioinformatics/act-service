'''
Created on Feb 8, 2017

@author: mhalagan
'''


def search_hla_features(locus, gfe_feats):

    # TODO: Fix error observed with DQB1
    match = "MATCH(hla:IMGT)-[:HAS_GFE]-(gfe:GFE)-[f1:HAS_FEATURE]-(feat1:FEATURE)"
    if(len(gfe_feats)) > 1:
        for i in range(1, len(gfe_feats)):
            j = i + 1
            match = match + ",(hla:IMGT)-[:HAS_GFE]-(gfe:GFE)-[f" + str(j) + ":HAS_FEATURE]-(feat" + str(j) + ":FEATURE)"

    i = 1
    feat_q = "WHERE hla.locus = \"" + locus + "\""

    for feat in gfe_feats:
        [term, rank] = feat.split("-")
        feat_q = feat_q + " AND feat" + str(i) + ".name = \"" + term.upper() + "\""
        feat_q = feat_q + " AND feat" + str(i) + ".rank = \"" + rank + "\""
        acc_q = " AND( f" + str(i) + ".accession = \"" + gfe_feats[feat][0] + "\""
        if(len(gfe_feats[feat]) == 1):
            acc_q = acc_q + ")"
        else:
            for j in range(1, len(gfe_feats[feat])):
                acc_q = acc_q + " OR f" + str(i) + ".accession = \"" + gfe_feats[feat][j] + "\""
            acc_q = acc_q + ")"
        feat_q = feat_q + acc_q
        i += 1

    return_q = " RETURN DISTINCT hla.name AS HLA, gfe.name AS GFE"
    return match + feat_q + return_q


def fullseqid(seq):
    q1 = "MATCH(seq:SEQUENCE)"
    q2 = " WHERE seq.sequence = \"" + seq + "\""
    q3 = " RETURN seq.sequenceId AS ID"
    return q1 + q2 + q3


def seqid(seq, seqtype, rank):
    q1 = "MATCH(feat:FEATURE)"
    q2 = " WHERE feat.sequence = \"" + seq + "\""
    q3 = " AND feat.rank = \"" + str(rank) + "\""
    q4 = " AND feat.name = \"" + seqtype.upper() + "\""
    q5 = "RETURN feat.sequenceId AS ID"
    return q1 + q2 + q3 + q4 + q5


def hla_alleleid(hla):
    q1 = "MATCH(hla:IMGT) "
    q2 = "WHERE hla.name = \"" + hla + "\" "
    q3 = "RETURN hla.alleleId AS ID"
    return q1 + q2 + q3


def gfe_alleleid(gfe):
    q1 = "MATCH(gfe:GFE) "
    q2 = "WHERE gfe.name = \"" + gfe + "\" "
    q3 = "RETURN gfe.alleleId AS ID"
    return q1 + q2 + q3


def ref_query(alleles):
    q1 = "MATCH(hla:IMGT)-[:HAS_GFE]-(gfe:GFE) "
    q2 = "WHERE hla.name = \"" + alleles[0] + "\""
    for i in range(1, len(alleles)):
        q2 = q2 + " OR hla.name = \"" + alleles[i] + "\""
    q3 = " RETURN hla.name AS HLA, gfe.name AS GFE"
    return q1 + q2 + q3


def get_features(gfe):

    q1 = "MATCH(gfe:GFE)-[f1:HAS_FEATURE]-(f:FEATURE)"
    q2 = "WHERE gfe.name = \"" + gfe + "\""
    q3 = " AND NOT(f1.status = \"persisted\")"
    q4 = "RETURN f.name AS term,f.rank AS rank,f1.accession AS accession,f.sequence AS sequence"
    return q1 + q2 + q3 + q4


def sequence_search(locus, sequence):
    seq_query = "MATCH (hla:IMGT)-[f1:HAS_GFE]-(gfe:GFE)-[f2:HAS_FEATURE]-(seq:SEQUENCE)"
    seq_query1 = " WHERE seq.sequence = \"" + sequence + "\""
    seq_query2 = " AND hla.locus = \"" + locus + "\""
    seq_query3 = " AND NOT(f1.status = \"persisted\")"
    seq_query4 = " AND NOT(f2.status = \"persisted\")"
    seq_query5 = " RETURN hla.name AS HLA, gfe.name AS GFE"
    query = seq_query + seq_query1 + seq_query2 + seq_query3 + seq_query4 + seq_query5
    return(query)


def gfe_search(gfe):
    seq_query = "MATCH (hla:IMGT)-[f:HAS_GFE]-(gfe:GFE)"
    seq_query1 = " WHERE gfe.name = \"" + gfe + "\""
    seq_query2 = " AND NOT(f.status = \"persisted\")"
    seq_query3 = " RETURN hla.name AS HLA"
    query = seq_query + seq_query1 + seq_query2 + seq_query3
    return(query)


def similar_gfe_classI(gfe, exon2, exon3):

    [locus, feature_accessions] = gfe.split("w")
    typing_match = "MATCH (gfe1:GFE)-[f1:HAS_FEATURE]->(feat1:FEATURE)," \
        + " (gfe2:GFE)-[f2:HAS_FEATURE]->(feat1:FEATURE)," \
        + " (gfe1:GFE)-[f3:HAS_FEATURE]->(feat2:FEATURE)," \
        + " (gfe2:GFE)-[f4:HAS_FEATURE]->(feat2:FEATURE)," \
        + " (gfe1:GFE)-[f5:HAS_FEATURE]->(feat3:FEATURE)," \
        + " (gfe2:GFE)-[f6:HAS_FEATURE]->(feat3:FEATURE)" \
        + " WHERE gfe1.locus = \"" + locus + "\"" \
        + " AND gfe2.locus = gfe1.locus" \
        + " AND f1.accession = f2.accession" \
        + " AND f3.accession = f4.accession" \
        + " AND f1.accession = \"" + exon2 + "\"" \
        + " AND f3.accession = \"" + exon3 + "\"" \
        + " AND feat1.name = \"EXON\"" \
        + " AND feat2.name = \"EXON\"" \
        + " AND feat1.rank = \"2\"" \
        + " AND feat2.rank = \"3\"" \
        + " AND NOT(f1.status = \"persisted\")" \
        + " AND NOT(f2.status = \"persisted\")" \
        + " AND NOT(f3.status = \"persisted\")" \
        + " AND NOT(f4.status = \"persisted\")" \
        + " AND NOT(f5.status = \"persisted\")" \
        + " AND NOT(f6.status = \"persisted\")" \
        + " AND f5.accession = f6.accession" \
        + " WITH gfe1.name AS GFE1, gfe2.name AS GFE2,collect(DISTINCT feat3.name) AS Names," \
        + " collect(DISTINCT {accesion:f6.accession,rank:feat3.rank, name: feat3.name}) AS Accession " \
        + " return GFE1,GFE2,Names,Accession,size(Accession) AS Count" \
        + " ORDER BY Count DESC"

    return(typing_match)


def similar_kir(locus):
    query = "MATCH(gfe:GFE)" \
        + "WHERE gfe.locus = \"" + locus + "\"" \
        + "RETURN gfe.name AS GFE"
    return query


def similar_gfe_classII(gfe, exon2):

    [locus, feature_accessions] = gfe.split("w")
    typing_match = "MATCH (gfe1:GFE)-[f1:HAS_FEATURE]->(feat1:FEATURE)" \
        + " WHERE gfe1.locus = \"" + locus + "\"" \
        + " AND f1.accession = \"" + exon2 + "\"" \
        + " AND feat1.name = \"EXON\"" \
        + " AND feat1.rank = \"2\"" \
        + " AND NOT(f1.status = \"persisted\")" \
        + " WITH gfe1.name AS GFE1, collect(DISTINCT feat1.name) AS Names," \
        + " collect(DISTINCT {accesion:f1.accession,rank:feat1.rank, name: feat1.name}) AS Accession " \
        + " return GFE1,Names,Accession,size(Accession) AS Count" \
        + " ORDER BY Count DESC"

    return(typing_match)


def hla_Ggroups(hla):
    matchq = "MATCH (G:G_GROUP)-[:IN_GROUP]-(hla:IMGT) WHERE hla.name = \"" + hla + "\""
    returnq = " RETURN DISTINCT G.name as G_GROUP"
    cyper_query = matchq + returnq
    return(cyper_query)


def gfe_Ggroups(gfe):
    matchq = "MATCH (G:G_GROUP)-[:IN_GROUP]-(hla:IMGT)-[:HAS_GFE]-(gfe1:GFE) WHERE gfe1.name = \"" + gfe + "\""
    returnq = " RETURN DISTINCT G.name as G_GROUP"
    cyper_query = matchq + returnq
    return(cyper_query)


def gfe_hla(gfe):
    matchq = "MATCH (hla:IMGT)-[h:HAS_GFE]-(gfe1:GFE) WHERE gfe1.name = \"" + gfe + "\""
    unq = " AND NOT(h.status = \"persisted\") "
    withq = " WITH collect(DISTINCT hla.name) as HLA"
    returnq = " RETURN HLA"
    cyper_query = matchq + unq + withq + returnq
    return(cyper_query)


# def groups_classI(locus, group, exon2, exon3):

#     match = "MATCH (feat2:FEATURE)-[f2:HAS_FEATURE]-(gfe:GFE)-[:HAS_GFE]-(hla:IMGT)-[:IN_GROUP]-(group:" + group + ")-[:IN_GROUP]-(hla:IMGT)-[:HAS_GFE]-(gfe:GFE)-[f1:HAS_FEATURE]-(feat1:FEATURE)"
#     where = " WHERE feat1.rank = '2' AND feat2.name = \"EXON\" AND feat1.name = \"EXON\" AND hla.locus = \"" + locus + "\""
#     where2 = " AND feat2.rank = '3' AND f2.accession = \"" + exon3 + "\" AND f1.accession = \"" + exon2 + "\""
#     returnc = "  RETURN group.name as ARS, hla.name as HLA,gfe.name AS GFE,group.name AS ARS"
#     cypher = match + where + where2 + returnc
#     return(cypher)


# def groups_classII(locus, group, exon3):

#     match = "MATCH (group:" + group + ")-[:IN_GROUP]-(hla:IMGT)-[:HAS_GFE]-(gfe:GFE)-[f1:HAS_FEATURE]->(feat1:FEATURE)"
#     where = " WHERE feat1.rank = '3' AND feat1.name = \"EXON\" AND hla.locus = \"" + locus + "\""
#     where2 = " AND f1.accession = \"" + exon3 + "\""
#     returnc = " RETURN group.name as ARS, hla.name as HLA,gfe.name AS GFE,group.name AS ARS"
#     cypher = match + where + where2 + returnc
#     return(cypher)


def groups_classI(locus, group, exon2, exon3):

    match = "MATCH (group:" + group + ")-[:IN_GROUP]-(hla:IMGT)-[:HAS_GFE]-(gfe1:GFE)-[f1:HAS_FEATURE]->(feat1:FEATURE)"
    match2 = ",(gfe1:GFE)-[f2:HAS_FEATURE]->(feat2:FEATURE) "
    where = " WHERE feat1.rank = '2' AND feat2.name = \"EXON\" AND feat1.name = \"EXON\" AND hla.locus = \"" + locus + "\""
    where2 = " AND feat2.rank = '3' AND f2.accession = \"" + exon3 + "\" AND f1.accession = \"" + exon2 + "\""
    withst = " WITH collect(DISTINCT hla.name) as HLA,collect(DISTINCT gfe1.name) as GFE,collect(DISTINCT group.name) as ARS"
    returnc = " RETURN HLA,GFE,ARS"
    cypher = match + match2 + where + where2 + withst + returnc
    print(cypher)
    return(cypher)


def groups_classII(locus, group, exon3):

    match = "MATCH (group:" + group + ")-[:IN_GROUP]-(hla:IMGT)-[:HAS_GFE]-(gfe1:GFE)-[f1:HAS_FEATURE]->(feat1:FEATURE)"
    match2 = ",(gfe1:GFE)-[f2:HAS_FEATURE]->(feat2:FEATURE) "
    where = " WHERE feat1.rank = '3' AND feat1.name = \"EXON\" AND hla.locus = \"" + locus + "\""
    where2 = " AND f1.accession = \"" + exon3 + "\""
    withst = " WITH collect(DISTINCT hla.name) as HLA,collect(DISTINCT gfe1.name) as GFE"
    returnc = " RETURN HLA,GFE"
    cypher = match + match2 + where + where2 + withst + returnc
    return(cypher)


def search_feature(term, rank, sequence):
    q1 = "MATCH(gfe:GFE)-[h:HAS_FEATURE]-(feat:FEATURE) WHERE feat.name = \"" + term.upper() + "\" "
    q2 = "AND feat.rank = \"" + str(rank) + "\" "
    q3 = "AND feat.sequence = \"" + sequence + "\" "
    unq = " AND NOT(h.status = \"persisted\")"
    q4 = "RETURN feat.name"
    return q1 + q2 + q3 + unq + q4


def hla_gfe(hla):
    seq_query = "MATCH (hla:IMGT)-[h:HAS_GFE]-(gfe:GFE)"
    seq_query1 = " WHERE hla.name = \"" + hla + "\""
    unq = " AND NOT(h.status = \"persisted\") "
    seq_query2 = " RETURN DISTINCT gfe.name AS GFE"
    query = seq_query + seq_query1 + unq + seq_query2
    return(query)


def hla_ars(group, hla):
    matchq = "MATCH (group:" + group + ")-[:IN_GROUP]-(hla:IMGT)-[h:HAS_GFE]-(g:GFE) WHERE hla.name = \"" + hla + "\""
    unq = " AND NOT(h.status = \"persisted\")"
    returnq = " RETURN group.name as ARS, hla.name as HLA,g.name AS GFE"
    cyper_query = matchq + unq + returnq
    return(cyper_query)


def gfe_ars(group, gfe):
    matchq = "MATCH (group:" + group + ")-[:IN_GROUP]-(hla:IMGT)-[:HAS_GFE]-(g:GFE) WHERE g.name = \"" + gfe + "\""
    returnq = " RETURN group.name as ARS, hla.name as HLA,g.name AS GFE"
    cyper_query = matchq + returnq
    return(cyper_query)


def get_sequence(seqtype, allele):
    seq_query = "MATCH (allele:" + seqtype + ")-[h:HAS_FEATURE]-(seq:SEQUENCE)"
    seq_query1 = " WHERE allele.name = \"" + allele + "\""
    unq = " AND NOT(h.status = \"persisted\")"
    seq_query2 = " RETURN DISTINCT seq.sequence AS SEQ"
    query = seq_query + seq_query1 + unq + seq_query2
    return(query)


def persisted_query():
    query = "MATCH(hla:IMGT)-[g:HAS_GFE]-(gfe:GFE)-[f:HAS_FEATURE]-(feat:FEATURE) "
    q2 = "WHERE g.status = \"persisted\" "
    q3 = "AND g.status = \"persisted\" "
    q4 = "RETURN hla.name AS HLA,gfe.name AS GFE,feat.name AS TERM,feat.rank AS RANK,f.accession AS ACCESSION,feat.sequence AS SEQUENCE"
    return query + q2 + q3 + q4

