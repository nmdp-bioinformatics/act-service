# coding: utf-8

from __future__ import absolute_import

from swagger_server.models.allele_call import AlleleCall
from swagger_server.models.error import Error
from . import BaseTestCase
from six import BytesIO
from flask import json


class TestDefaultController(BaseTestCase):
    """ DefaultController integration test stubs """

    def test_hla_post(self):
        """
        Test case for hla_post

        
        """
        query_string = [('locus', 'locus_example'),
                        ('sequence', 'sequence_example'),
                        ('neo4j_url', 'neo4j_url_example'),
                        ('gfe_url', 'gfe_url_example'),
                        ('verbose', true)]
        response = self.client.open('/hla',
                                    method='POST',
                                    query_string=query_string)
        self.assert200(response, "Response body is : " + response.data.decode('utf-8'))


if __name__ == '__main__':
    import unittest
    unittest.main()
