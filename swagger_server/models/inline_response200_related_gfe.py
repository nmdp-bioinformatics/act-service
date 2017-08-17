# coding: utf-8

from __future__ import absolute_import
from swagger_server.models.inline_response200_shares import InlineResponse200Shares
from .base_model_ import Model
from datetime import date, datetime
from typing import List, Dict
from ..util import deserialize_model


class InlineResponse200RelatedGfe(Model):
    """
    NOTE: This class is auto generated by the swagger code generator program.
    Do not edit the class manually.
    """
    def __init__(self, gfe: str=None, features_shared: int=None, shares: List[InlineResponse200Shares]=None):
        """
        InlineResponse200RelatedGfe - a model defined in Swagger

        :param gfe: The gfe of this InlineResponse200RelatedGfe.
        :type gfe: str
        :param features_shared: The features_shared of this InlineResponse200RelatedGfe.
        :type features_shared: int
        :param shares: The shares of this InlineResponse200RelatedGfe.
        :type shares: List[InlineResponse200Shares]
        """
        self.swagger_types = {
            'gfe': str,
            'features_shared': int,
            'shares': List[InlineResponse200Shares]
        }

        self.attribute_map = {
            'gfe': 'gfe',
            'features_shared': 'features_shared',
            'shares': 'shares'
        }

        self._gfe = gfe
        self._features_shared = features_shared
        self._shares = shares

    @classmethod
    def from_dict(cls, dikt) -> 'InlineResponse200RelatedGfe':
        """
        Returns the dict as a model

        :param dikt: A dict.
        :type: dict
        :return: The inline_response_200_related_gfe of this InlineResponse200RelatedGfe.
        :rtype: InlineResponse200RelatedGfe
        """
        return deserialize_model(dikt, cls)

    @property
    def gfe(self) -> str:
        """
        Gets the gfe of this InlineResponse200RelatedGfe.

        :return: The gfe of this InlineResponse200RelatedGfe.
        :rtype: str
        """
        return self._gfe

    @gfe.setter
    def gfe(self, gfe: str):
        """
        Sets the gfe of this InlineResponse200RelatedGfe.

        :param gfe: The gfe of this InlineResponse200RelatedGfe.
        :type gfe: str
        """

        self._gfe = gfe

    @property
    def features_shared(self) -> int:
        """
        Gets the features_shared of this InlineResponse200RelatedGfe.

        :return: The features_shared of this InlineResponse200RelatedGfe.
        :rtype: int
        """
        return self._features_shared

    @features_shared.setter
    def features_shared(self, features_shared: int):
        """
        Sets the features_shared of this InlineResponse200RelatedGfe.

        :param features_shared: The features_shared of this InlineResponse200RelatedGfe.
        :type features_shared: int
        """

        self._features_shared = features_shared

    @property
    def shares(self) -> List[InlineResponse200Shares]:
        """
        Gets the shares of this InlineResponse200RelatedGfe.

        :return: The shares of this InlineResponse200RelatedGfe.
        :rtype: List[InlineResponse200Shares]
        """
        return self._shares

    @shares.setter
    def shares(self, shares: List[InlineResponse200Shares]):
        """
        Sets the shares of this InlineResponse200RelatedGfe.

        :param shares: The shares of this InlineResponse200RelatedGfe.
        :type shares: List[InlineResponse200Shares]
        """

        self._shares = shares
