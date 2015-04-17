#!/usr/bin/env python

import sncosmo.models
import numpy


class SEDFileSource(sncosmo.models.TimeSeriesSource):
    """A TimeSeriesSource stored in a 3-column ASCII file format, for PHASE,
    LAMBDA, and F_LAMBDA.  The hash symbol # is a comment line.

    The spectral flux density of this model is given by

    .. math::

       F(t, \lambda) = A \\times M(t, \lambda)

    where _M_ is the flux defined on a grid in phase and wavelength and _A_
    (amplitude) is the single free parameter of the model. It should be noted
    that while t and \lambda are in the rest frame of the object, the flux
    density is defined at redshift zero. This means that for objects with the
    same intrinsic luminosity, the amplitude will be smaller for objects at
    larger luminosity distances.

    Parameters
    ----------
    filename : str
        Name of the filename that contains the Time Series
    zero_before : bool, optional
        If True, flux at phases before minimum phase will be zeroed. The
        default is False, in which case the flux at such phases will be equal
        to the flux at the minimum phase (``flux[0, :]`` in the input array).
    version : str, optional
        Version of the model. Default is `None`.

    Returns
    -------
        `~sncosmo.TimeSeriesSource` instance representing the TimeSeriesSource
        in file
    """

    _param_names = ['amplitude']
    param_names_latex = ['A']

    def __init__(self, filename, zero_before=False, version=None):
        phase, wave, flux = numpy.loadtxt(filename, unpack=True)

        # Convert 3 column format to that expected by TimeSeriesSource
        phase_u = numpy.unique(phase)
        wave_u = numpy.unique(wave)

        lenp = len(phase_u)
        lenw = len(wave_u)

        if lenp*lenw != len(flux):
            raise TypeError('File is not a TimeSeriesSource')

        i = numpy.zeros(len(flux), dtype='int')
        j = numpy.zeros(len(flux), dtype='int')
        for index, p in enumerate(phase_u):
            i[phase == p] = index
        for index, w in enumerate(wave_u):
            j[wave == w] = index

        flux = flux[i*lenw+j]
        flux = numpy.reshape(flux, (lenp, lenw))
        super(SEDFileSource, self).__init__(phase_u, wave_u, flux,
                                            zero_before=False,
                                            name=filename, version=None)

# filename = '/Users/akim/project/SNDATA_ROOT/snsed/NON1A/SDSS-019323.SED'

# data = SEDFileSource(filename)

sn = sncosmo.Model(source='snana-2007nc')
print sn.param_names
# wefwe
import matplotlib.pyplot as plt
plt.plot(data._wave, data.flux(0, data._wave))
plt.plot(sn.source._wave, sn.flux(0, sn.source._wave)*0.95)
plt.show()
