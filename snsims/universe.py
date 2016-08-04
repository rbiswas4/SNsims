"""
Class for describing a collection of SN
"""
from __future__ import absolute_import, print_function
from future.utils import with_metaclass
import abc

__all__ = ['Universe']

class Universe(with_metaclass(abc.ABCMeta, object)):

    @abstractproperty
    def snParams(self):
        """
        """
        pass

    @abstractproperty
    def SN(self):
        """
        """
        pass

