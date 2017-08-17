# coding: utf-8

from __future__ import absolute_import

from swagger_server.models.allele_call import AlleleCall
from swagger_server.models.ars_call import ArsCall
from swagger_server.models.error import Error
from swagger_server.models.feature_call import FeatureCall
from swagger_server.models.gfe_call import GfeCall
from swagger_server.models.sequence_call import SequenceCall
from . import BaseTestCase
from six import BytesIO
from flask import json


class TestDefaultController(BaseTestCase):
    """ DefaultController integration test stubs """

    def test_act_get(self):
        """
        Test case for act_get

        
        """
        query_string = [('locus', 'locus_example'),
                        ('sequence', 'sequence_example'),
                        ('gfe', 'gfe_example'),
                        ('neo4j_url', 'neo4j_url_example'),
                        ('user', 'user_example'),
                        ('password', 'password_example'),
                        ('gfe_url', 'gfe_url_example'),
                        ('verbose', true),
                        ('persist', true)]
        response = self.client.open('/act',
                                    method='GET',
                                    query_string=query_string)
        self.assert200(response, "Response body is : " + response.data.decode('utf-8'))

    def test_ars_get(self):
        """
        Test case for ars_get

        
        """
        query_string = [('allele', 'allele_example'),
                        ('group', 'group_example'),
                        ('neo4j_url', 'neo4j_url_example'),
                        ('user', 'user_example'),
                        ('password', 'password_example'),
                        ('gfe_url', 'gfe_url_example'),
                        ('verbose', true)]
        response = self.client.open('/ars',
                                    method='GET',
                                    query_string=query_string)
        self.assert200(response, "Response body is : " + response.data.decode('utf-8'))

    def test_feature_get(self):
        """
        Test case for feature_get

        
        """
        query_string = [('hla', 'hla_example'),
                        ('feature', 'feature_example'),
                        ('neo4j_url', 'neo4j_url_example'),
                        ('user', 'user_example'),
                        ('password', 'password_example'),
                        ('gfe_url', 'gfe_url_example'),
                        ('verbose', true)]
        response = self.client.open('/feature_search',
                                    method='GET',
                                    query_string=query_string)
        self.assert200(response, "Response body is : " + response.data.decode('utf-8'))

    def test_gfe_get(self):
        """
        Test case for gfe_get

        
        """
        query_string = [('hla', 'hla_example'),
                        ('neo4j_url', 'neo4j_url_example'),
                        ('user', 'user_example'),
                        ('password', 'password_example'),
                        ('gfe_url', 'gfe_url_example'),
                        ('verbose', true)]
        response = self.client.open('/gfe',
                                    method='GET',
                                    query_string=query_string)
        self.assert200(response, "Response body is : " + response.data.decode('utf-8'))

    def test_seqsrch_get(self):
        """
        Test case for seqsrch_get

        
        """
        query_string = [('locus', 'locus_example'),
                        ('start', 56),
                        ('end', 56),
                        ('hla', 'hla_example'),
                        ('feature', 'feature_example'),
                        ('neo4j_url', 'neo4j_url_example'),
                        ('user', 'user_example'),
                        ('password', 'password_example'),
                        ('gfe_url', 'gfe_url_example'),
                        ('verbose', true)]
        response = self.client.open('/seq_search',
                                    method='GET',
                                    query_string=query_string)
        self.assert200(response, "Response body is : " + response.data.decode('utf-8'))


if __name__ == '__main__':
    import unittest
    unittest.main()
