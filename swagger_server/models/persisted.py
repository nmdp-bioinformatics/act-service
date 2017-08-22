# coding: utf-8

from __future__ import absolute_import
from swagger_server.models.persisted_data import PersistedData
from .base_model_ import Model
from datetime import date, datetime
from typing import List, Dict
from ..util import deserialize_model


class Persisted(Model):
    """
    NOTE: This class is auto generated by the swagger code generator program.
    Do not edit the class manually.
    """
    def __init__(self, persisted_data: List[PersistedData]=None, act_version: str=None, gfedb_version: str=None):
        """
        Persisted - a model defined in Swagger

        :param persisted_data: The persisted_data of this Persisted.
        :type persisted_data: List[PersistedData]
        :param act_version: The act_version of this Persisted.
        :type act_version: str
        :param gfedb_version: The gfedb_version of this Persisted.
        :type gfedb_version: str
        """
        self.swagger_types = {
            'persisted_data': List[PersistedData],
            'act_version': str,
            'gfedb_version': str
        }

        self.attribute_map = {
            'persisted_data': 'persisted_data',
            'act_version': 'act_version',
            'gfedb_version': 'gfedb_version'
        }

        self._persisted_data = persisted_data
        self._act_version = act_version
        self._gfedb_version = gfedb_version

    @classmethod
    def from_dict(cls, dikt) -> 'Persisted':
        """
        Returns the dict as a model

        :param dikt: A dict.
        :type: dict
        :return: The Persisted of this Persisted.
        :rtype: Persisted
        """
        return deserialize_model(dikt, cls)

    @property
    def persisted_data(self) -> List[PersistedData]:
        """
        Gets the persisted_data of this Persisted.

        :return: The persisted_data of this Persisted.
        :rtype: List[PersistedData]
        """
        return self._persisted_data

    @persisted_data.setter
    def persisted_data(self, persisted_data: List[PersistedData]):
        """
        Sets the persisted_data of this Persisted.

        :param persisted_data: The persisted_data of this Persisted.
        :type persisted_data: List[PersistedData]
        """

        self._persisted_data = persisted_data

    @property
    def act_version(self) -> str:
        """
        Gets the act_version of this Persisted.

        :return: The act_version of this Persisted.
        :rtype: str
        """
        return self._act_version

    @act_version.setter
    def act_version(self, act_version: str):
        """
        Sets the act_version of this Persisted.

        :param act_version: The act_version of this Persisted.
        :type act_version: str
        """

        self._act_version = act_version

    @property
    def gfedb_version(self) -> str:
        """
        Gets the gfedb_version of this Persisted.

        :return: The gfedb_version of this Persisted.
        :rtype: str
        """
        return self._gfedb_version

    @gfedb_version.setter
    def gfedb_version(self, gfedb_version: str):
        """
        Sets the gfedb_version of this Persisted.

        :param gfedb_version: The gfedb_version of this Persisted.
        :type gfedb_version: str
        """

        self._gfedb_version = gfedb_version
