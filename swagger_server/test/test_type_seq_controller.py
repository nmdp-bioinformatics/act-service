# coding: utf-8

from __future__ import absolute_import

from flask import json
from six import BytesIO

from swagger_server.models.error import Error  # noqa: E501
from swagger_server.models.typing import Typing  # noqa: E501
from swagger_server.test import BaseTestCase


class TestTypeSeqController(BaseTestCase):
    """TypeSeqController integration test stubs"""

    def test_typeseq_get(self):
        """Test case for typeseq_get

        
        """
        query_string = [('locus', 'locus_example'),
                        ('sequence', 'sequence_example'),
                        ('imgthla_version', '3310'),
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
