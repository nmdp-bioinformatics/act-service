# coding: utf-8

"""
    Allele Calling Service

    The Allele Calling  service provides an API for converting raw sequence data to GFE and HLA allele calls.

    OpenAPI spec version: 0.0.4
    Contact: mhalagan@nmdp.org
    Generated by: https://github.com/swagger-api/swagger-codegen.git
"""


from pprint import pformat
from six import iteritems
import re


class PersistedData(object):
    """
    NOTE: This class is auto generated by the swagger code generator program.
    Do not edit the class manually.
    """


    """
    Attributes:
      swagger_types (dict): The key is attribute name
                            and the value is attribute type.
      attribute_map (dict): The key is attribute name
                            and the value is json key in definition.
    """
    swagger_types = {
        'hla': 'str',
        'gfe': 'str',
        'term': 'str',
        'rank': 'str',
        'accession': 'str',
        'sequence': 'str'
    }

    attribute_map = {
        'hla': 'hla',
        'gfe': 'gfe',
        'term': 'term',
        'rank': 'rank',
        'accession': 'accession',
        'sequence': 'sequence'
    }

    def __init__(self, hla=None, gfe=None, term=None, rank=None, accession=None, sequence=None):
        """
        PersistedData - a model defined in Swagger
        """

        self._hla = None
        self._gfe = None
        self._term = None
        self._rank = None
        self._accession = None
        self._sequence = None

        if hla is not None:
          self.hla = hla
        if gfe is not None:
          self.gfe = gfe
        if term is not None:
          self.term = term
        if rank is not None:
          self.rank = rank
        if accession is not None:
          self.accession = accession
        if sequence is not None:
          self.sequence = sequence

    @property
    def hla(self):
        """
        Gets the hla of this PersistedData.

        :return: The hla of this PersistedData.
        :rtype: str
        """
        return self._hla

    @hla.setter
    def hla(self, hla):
        """
        Sets the hla of this PersistedData.

        :param hla: The hla of this PersistedData.
        :type: str
        """

        self._hla = hla

    @property
    def gfe(self):
        """
        Gets the gfe of this PersistedData.

        :return: The gfe of this PersistedData.
        :rtype: str
        """
        return self._gfe

    @gfe.setter
    def gfe(self, gfe):
        """
        Sets the gfe of this PersistedData.

        :param gfe: The gfe of this PersistedData.
        :type: str
        """

        self._gfe = gfe

    @property
    def term(self):
        """
        Gets the term of this PersistedData.

        :return: The term of this PersistedData.
        :rtype: str
        """
        return self._term

    @term.setter
    def term(self, term):
        """
        Sets the term of this PersistedData.

        :param term: The term of this PersistedData.
        :type: str
        """

        self._term = term

    @property
    def rank(self):
        """
        Gets the rank of this PersistedData.

        :return: The rank of this PersistedData.
        :rtype: str
        """
        return self._rank

    @rank.setter
    def rank(self, rank):
        """
        Sets the rank of this PersistedData.

        :param rank: The rank of this PersistedData.
        :type: str
        """

        self._rank = rank

    @property
    def accession(self):
        """
        Gets the accession of this PersistedData.

        :return: The accession of this PersistedData.
        :rtype: str
        """
        return self._accession

    @accession.setter
    def accession(self, accession):
        """
        Sets the accession of this PersistedData.

        :param accession: The accession of this PersistedData.
        :type: str
        """

        self._accession = accession

    @property
    def sequence(self):
        """
        Gets the sequence of this PersistedData.

        :return: The sequence of this PersistedData.
        :rtype: str
        """
        return self._sequence

    @sequence.setter
    def sequence(self, sequence):
        """
        Sets the sequence of this PersistedData.

        :param sequence: The sequence of this PersistedData.
        :type: str
        """

        self._sequence = sequence

    def to_dict(self):
        """
        Returns the model properties as a dict
        """
        result = {}

        for attr, _ in iteritems(self.swagger_types):
            value = getattr(self, attr)
            if isinstance(value, list):
                result[attr] = list(map(
                    lambda x: x.to_dict() if hasattr(x, "to_dict") else x,
                    value
                ))
            elif hasattr(value, "to_dict"):
                result[attr] = value.to_dict()
            elif isinstance(value, dict):
                result[attr] = dict(map(
                    lambda item: (item[0], item[1].to_dict())
                    if hasattr(item[1], "to_dict") else item,
                    value.items()
                ))
            else:
                result[attr] = value

        return result

    def to_str(self):
        """
        Returns the string representation of the model
        """
        return pformat(self.to_dict())

    def __repr__(self):
        """
        For `print` and `pprint`
        """
        return self.to_str()

    def __eq__(self, other):
        """
        Returns true if both objects are equal
        """
        if not isinstance(other, PersistedData):
            return False

        return self.__dict__ == other.__dict__

    def __ne__(self, other):
        """
        Returns true if both objects are not equal
        """
        return not self == other
