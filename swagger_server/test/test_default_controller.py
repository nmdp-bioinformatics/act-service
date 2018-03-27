# coding: utf-8

from __future__ import absolute_import

from flask import json
from six import BytesIO

from swagger_server.models.error import Error  # noqa: E501
from swagger_server.models.seqdiff import Seqdiff  # noqa: E501
from swagger_server.models.typing import Typing  # noqa: E501
from swagger_server.test import BaseTestCase


class TestDefaultController(BaseTestCase):
    """DefaultController integration test stubs"""

    def test_actformat_get(self):
        """Test case for actformat_get

        
        """
        query_string = [('locus', 'locus_example'),
                        ('sequence', 'sequence_example'),
                        ('gfe', 'gfe_example'),
                        ('format_type', 'format_type_example'),
                        ('neo4j_url', 'neo4j_url_example'),
                        ('user', 'user_example'),
                        ('password', 'password_example'),
                        ('gfe_url', 'gfe_url_example'),
                        ('verbose', true),
                        ('persist', true)]
        response = self.client.open(
            '/type_format',
            method='GET',
            query_string=query_string)
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_ars_get(self):
        """Test case for ars_get

        
        """
        query_string = [('positions', 'positions_example'),
                        ('alleles', 'alleles_example'),
                        ('neo4j_url', 'neo4j_url_example'),
                        ('user', 'user_example'),
                        ('password', 'password_example'),
                        ('verbose', true)]
        response = self.client.open(
            '/sequence_variations',
            method='GET',
            query_string=query_string)
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_closest_gfe_get(self):
        """Test case for closest_gfe_get

        
        """
        query_string = [('positions', 'positions_example'),
                        ('alleles', 'alleles_example'),
                        ('neo4j_url', 'neo4j_url_example'),
                        ('user', 'user_example'),
                        ('password', 'password_example'),
                        ('verbose', true)]
        response = self.client.open(
            '/closest_gfe',
            method='GET',
            query_string=query_string)
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_gfe_get(self):
        """Test case for gfe_get

        
        """
        query_string = [('hla', 'hla_example'),
                        ('neo4j_url', 'neo4j_url_example'),
                        ('user', 'user_example'),
                        ('password', 'password_example'),
                        ('verbose', true)]
        response = self.client.open(
            '/hla2gfe',
            method='GET',
            query_string=query_string)
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_seqdiff_get(self):
        """Test case for seqdiff_get

        
        """
        query_string = [('allele1', 'allele1_example'),
                        ('allele2', 'allele2_example'),
                        ('neo4j_url', 'neo4j_url_example'),
                        ('user', 'user_example'),
                        ('password', 'password_example'),
                        ('verbose', true)]
        response = self.client.open(
            '/seqdiff/sequences',
            method='GET',
            query_string=query_string)
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_seqdiffalleles_get(self):
        """Test case for seqdiffalleles_get

        
        """
        query_string = [('allele1', 'allele1_example'),
                        ('allele2', 'allele2_example'),
                        ('neo4j_url', 'neo4j_url_example'),
                        ('user', 'user_example'),
                        ('password', 'password_example'),
                        ('verbose', true)]
        response = self.client.open(
            '/seqdiff/alleles',
            method='GET',
            query_string=query_string)
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_sequence_search_get(self):
        """Test case for sequence_search_get

        
        """
        query_string = [('positions', 'positions_example'),
                        ('alleles', 'alleles_example'),
                        ('neo4j_url', 'neo4j_url_example'),
                        ('user', 'user_example'),
                        ('password', 'password_example'),
                        ('verbose', true)]
        response = self.client.open(
            '/sequence_search',
            method='GET',
            query_string=query_string)
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_typegfe_get(self):
        """Test case for typegfe_get

        
        """
        query_string = [('locus', 'locus_example'),
                        ('gfe', 'gfe_example'),
                        ('imgthla_version', 'Latest'),
                        ('neo4j_url', 'neo4j_url_example'),
                        ('user', 'user_example'),
                        ('password', 'password_example'),
                        ('verbose', true)]
        response = self.client.open(
            '/type_gfe',
            method='GET',
            query_string=query_string)
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))

    def test_typeseq_get(self):
        """Test case for typeseq_get

        
        """
        query_string = [('locus', 'locus_example'),
                        ('sequence', 'sequence_example'),
                        ('imgthla_version', 'Latest'),
                        ('neo4j_url', 'neo4j_url_example'),
                        ('user', 'user_example'),
                        ('password', 'password_example'),
                        ('verbose', true)]
        response = self.client.open(
            '/type_seq',
            method='GET',
            query_string=query_string)
        self.assert200(response,
                       'Response body is : ' + response.data.decode('utf-8'))


if __name__ == '__main__':
    import unittest
    unittest.main()
