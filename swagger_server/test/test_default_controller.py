# coding: utf-8

from __future__ import absolute_import

from swagger_server.models.allele_call import AlleleCall
from swagger_server.models.ars_call import ArsCall
from swagger_server.models.error import Error
from swagger_server.models.sequence import Sequence
from . import BaseTestCase
from six import BytesIO
from flask import json


class TestDefaultController(BaseTestCase):
    """ DefaultController integration test stubs """

    def test_act_post(self):
        """
        Test case for act_post

        
        """
        query_string = [('locus', 'locus_example'),
                        ('sequence', 'sequence_example'),
                        ('neo4j_url', 'neo4j_url_example'),
                        ('user', 'user_example'),
                        ('password', 'password_example'),
                        ('gfe_url', 'gfe_url_example'),
                        ('verbose', true)]
        response = self.client.open('/act',
                                    method='POST',
                                    query_string=query_string)
        self.assert200(response, "Response body is : " + response.data.decode('utf-8'))

    def test_ars_post(self):
        """
        Test case for ars_post

        
        """
        query_string = [('allele', 'allele_example'),
                        ('group', 'group_example'),
                        ('neo4j_url', 'neo4j_url_example'),
                        ('user', 'user_example'),
                        ('password', 'password_example'),
                        ('gfe_url', 'gfe_url_example'),
                        ('verbose', true)]
        response = self.client.open('/ars',
                                    method='POST',
                                    query_string=query_string)
        self.assert200(response, "Response body is : " + response.data.decode('utf-8'))

    def test_gfe_post(self):
        """
        Test case for gfe_post

        
        """
        query_string = [('hla', 'hla_example'),
                        ('neo4j_url', 'neo4j_url_example'),
                        ('user', 'user_example'),
                        ('password', 'password_example'),
                        ('gfe_url', 'gfe_url_example'),
                        ('verbose', true)]
        response = self.client.open('/gfe',
                                    method='POST',
                                    query_string=query_string)
        self.assert200(response, "Response body is : " + response.data.decode('utf-8'))

    def test_hla_post(self):
        """
        Test case for hla_post

        
        """
        query_string = [('gfe', 'gfe_example'),
                        ('neo4j_url', 'neo4j_url_example'),
                        ('user', 'user_example'),
                        ('password', 'password_example'),
                        ('gfe_url', 'gfe_url_example'),
                        ('verbose', true)]
        response = self.client.open('/hla',
                                    method='POST',
                                    query_string=query_string)
        self.assert200(response, "Response body is : " + response.data.decode('utf-8'))

    def test_sequence_post(self):
        """
        Test case for sequence_post

        
        """
        query_string = [('allele', 'allele_example'),
                        ('allele_type', 'allele_type_example'),
                        ('neo4j_url', 'neo4j_url_example'),
                        ('user', 'user_example'),
                        ('password', 'password_example'),
                        ('gfe_url', 'gfe_url_example'),
                        ('verbose', true)]
        response = self.client.open('/sequence',
                                    method='POST',
                                    query_string=query_string)
        self.assert200(response, "Response body is : " + response.data.decode('utf-8'))


if __name__ == '__main__':
    import unittest
    unittest.main()
