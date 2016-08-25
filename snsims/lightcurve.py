"""
Class for holding Light Curve Data
"""
from __future__ import absolute_import, print_function, division
from future.utils import with_metaclass
import abc
import pandas as pd
from astropy.table import Table
from .aliases import aliasDictionary


__all__ = ['BaseLightCurve']


class BaseLightCurve(with_metaclass(abc.ABCMeta, object)):
    """
    Abstract Base Class for Light Curve Data showing methods that need to be
    implemented.
    """
    @abc.abstractmethod
    def __init__(self):
        pass

    @abc.abstractproperty
    def lightCurve(self):
        """
        `pd.DataFrame` holding the lightCurve information. There can be more
        columns, but the following columns are mandatory:
        ['mjd', 'band', 'flux', 'fluxerr', 'zp', 'zpsys']
        """
        pass

    @abc.abstractproperty
    def snCosmoLC(self, coaddTimes=None):
        pass

    @abc.abstractmethod
    def coaddedLC(self, coaddTimes=None, format, *args, **kwargs):
        pass

    def missingColumns(self, lcdf):

        notFound = self.mandatoryColumns - set(lcdf.columns)
        return notFound

    @property
    def mandatoryColumns(self):
        """
        A list of mandatory columns in the light curve dataFrame with
        possible aliases in `self.mandatoryColumnAliases`.

        mjd : time
        band : string
        flux : model flux
        """
        reqd = set(['mjd', 'band', 'flux', 'fluxerr', 'zp', 'zpsys'])
        return reqd

    @property
    def mandatoryColumnAliases(self):
        """
        dictionary that maps standard names as keys to a possible set of
        aliases
        """
        aliases = {}
        aliases['mjd'] = ['time', 'expMJD']
        aliases['band'] = ['filter']
        aliases['fluxerr'] = ['flux_err']
        return aliases
