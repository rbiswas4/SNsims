"""
Class for describing a collection of SN
"""
from __future__ import absolute_import, print_function
from future.utils import with_metaclass
import abc

__all__ = ['Universe', 'HomogeneousSNUniverse']

class Universe(with_metaclass(abc.ABCMeta, object)):

    @abc.abstractproperty
    def snParams(self):
        """
        """
        pass

    @abc.abstractproperty
    def SN(self):
        """
        """
        pass


class HomogeneousSNUniverse(with_metaclass(abc.ABCMeta, Universe)):

    @abc.abstractproperty
    def rates(self):
        pass
    
    @abc.abstractproperty
    def paramDistribution(self):
        pass
