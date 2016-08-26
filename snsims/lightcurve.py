"""
Class for holding Light Curve Data
"""
from __future__ import absolute_import, print_function, division
from future.utils import with_metaclass
import abc
import pandas as pd
from astropy.table import Table
from .aliases import aliasDictionary


__all__ = ['BaseLightCurve', 'LightCurve']


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

    @abc.abstractmethod
    def snCosmoLC(self, coaddTimes=None):
        pass

    @abc.abstractmethod
    def coaddedLC(self, coaddTimes=None, format=None, *args, **kwargs):
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
    def columnAliases(self):
        """
        dictionary that maps standard names as keys to a possible set of
        aliases
        """
        aliases = {}
        aliases['mjd'] = ['time', 'expMJD']
        aliases['band'] = ['filter', 'filterName', 'bandName', 'bands']
        aliases['fluxerr'] = ['flux_err', 'flux_errs', 'fluxError']
        return aliases

class LightCurve(BaseLightCurve):

    def __init__(self, lcdf):
        """
        Instantiate Light Curve class

        Parameters
        ----------
        lcdf : `pd.DataFrame`, mandatory
            light curve information
        """
        self._lightCurve  = lcdf


    def missingColumns(self, lcdf):

        notFound = self.mandatoryColumns - set(lcdf.columns)
        return notFound

    @property
    def lightCurve(self):
        """
        The lightcurve in native format
        """
        # light curve
        _lc = self._lightCurve

        # Rename columns to standard names if necessary
        aliases = self.columnAliases
        standardNamingDict = aliasDictionary(_lc.columns, aliases)
        if len(standardNamingDict) > 0:
            _lc.rename(columns=standardNamingDict, inplace=True)

        # If all  mandatory columns exist return the light curve, or
        # raise ValueError citing missing columns
        missingColumns = self.missingColumns(_lc)
        if len(missingColumns) > 0:
            raise ValueError('light curve data has missing columns',
                             missingColumns)
        else:
            return _lc

    def snCosmoLC(self, coaddTimes=None):
        lc = self.lightCurve.rename(columns=dict(mjd='time'))
        return Table.from_pandas(lc)

    def coaddedLC(self, coaddTimes=None):
        pass

