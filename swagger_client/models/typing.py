# coding: utf-8

"""
    Gene Feature Enumeration Service

    The Gene Feature Enumeration (GFE) Submission service provides an API for converting raw sequence data to GFE. It provides both a RESTful API and a simple user interface for converting raw sequence data to GFE results. Sequences can be submitted one at a time or as a fasta file. This service uses <a href=\"https://github.com/nmdp-bioinformatics/service-feature\">nmdp-bioinformatics/service-feature</a> for encoding the raw sequence data and <a href=\"https://github.com/nmdp-bioinformatics/HSA\">nmdp-bioinformatics/HSA</a> for aligning the raw sequence data. The code is open source, and available on <a href=\"https://github.com/nmdp-bioinformatics/service-gfe-submission\">GitHub</a>.<br><br>Go to <a href=\"http://service-gfe-submission.readthedocs.io\">service-gfe-submission.readthedocs.io</a> for more information

    OpenAPI spec version: 1.0.7
    Contact: mhalagan@nmdp.org
    Generated by: https://github.com/swagger-api/swagger-codegen.git
"""


from pprint import pformat
from six import iteritems
import re


class Typing(object):
    """
    NOTE: This class is auto generated by the swagger code generator program.
    Do not edit the class manually.
    """
    def __init__(self, aligned=None, fullgene=None, gfe=None, imgthla=None, log=None, structure=None):
        """
        Typing - a model defined in Swagger

        :param dict swaggerTypes: The key is attribute name
                                  and the value is attribute type.
        :param dict attributeMap: The key is attribute name
                                  and the value is json key in definition.
        """
        self.swagger_types = {
            'aligned': 'float',
            'fullgene': 'Structure',
            'gfe': 'str',
            'imgthla': 'str',
            'log': 'list[str]',
            'structure': 'list[Structure]'
        }

        self.attribute_map = {
            'aligned': 'aligned',
            'fullgene': 'fullgene',
            'gfe': 'gfe',
            'imgthla': 'imgthla',
            'log': 'log',
            'structure': 'structure'
        }

        self._aligned = aligned
        self._fullgene = fullgene
        self._gfe = gfe
        self._imgthla = imgthla
        self._log = log
        self._structure = structure

    @property
    def aligned(self):
        """
        Gets the aligned of this Typing.

        :return: The aligned of this Typing.
        :rtype: float
        """
        return self._aligned

    @aligned.setter
    def aligned(self, aligned):
        """
        Sets the aligned of this Typing.

        :param aligned: The aligned of this Typing.
        :type: float
        """

        self._aligned = aligned

    @property
    def fullgene(self):
        """
        Gets the fullgene of this Typing.

        :return: The fullgene of this Typing.
        :rtype: Structure
        """
        return self._fullgene

    @fullgene.setter
    def fullgene(self, fullgene):
        """
        Sets the fullgene of this Typing.

        :param fullgene: The fullgene of this Typing.
        :type: Structure
        """

        self._fullgene = fullgene

    @property
    def gfe(self):
        """
        Gets the gfe of this Typing.

        :return: The gfe of this Typing.
        :rtype: str
        """
        return self._gfe

    @gfe.setter
    def gfe(self, gfe):
        """
        Sets the gfe of this Typing.

        :param gfe: The gfe of this Typing.
        :type: str
        """
        if gfe is None:
            raise ValueError("Invalid value for `gfe`, must not be `None`")

        self._gfe = gfe

    @property
    def imgthla(self):
        """
        Gets the imgthla of this Typing.

        :return: The imgthla of this Typing.
        :rtype: str
        """
        return self._imgthla

    @imgthla.setter
    def imgthla(self, imgthla):
        """
        Sets the imgthla of this Typing.

        :param imgthla: The imgthla of this Typing.
        :type: str
        """

        self._imgthla = imgthla

    @property
    def log(self):
        """
        Gets the log of this Typing.

        :return: The log of this Typing.
        :rtype: list[str]
        """
        return self._log

    @log.setter
    def log(self, log):
        """
        Sets the log of this Typing.

        :param log: The log of this Typing.
        :type: list[str]
        """

        self._log = log

    @property
    def structure(self):
        """
        Gets the structure of this Typing.

        :return: The structure of this Typing.
        :rtype: list[Structure]
        """
        return self._structure

    @structure.setter
    def structure(self, structure):
        """
        Sets the structure of this Typing.

        :param structure: The structure of this Typing.
        :type: list[Structure]
        """

        self._structure = structure

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
        return self.__dict__ == other.__dict__

    def __ne__(self, other):
        """
        Returns true if both objects are not equal
        """
        return not self == other
