import sncosmo
import numpy
import sys


def registry_types():
    """
    get the set of unique supernova types in the SNCosmo registry.

    Parameters
    ----------
    None

    Returns
    -------
    set of strings that correpsond to SN types.


    Examples
    --------
    >>> print sources.registry_types()
    >>> set(['PopIII', 'SN IIP', 'SN IIL/P', 'SN II-pec', 'SN Ic', 'SN Ib', 'SN Ia', 'SN IIn', 'SN Ib/c', 'SN IIL'])

    """
    sources = sncosmo.registry.get_loaders_metadata(sncosmo.Source)
    types = [s['type'] for s in sources]
    return set(types)


def registry_sources(sntype, subclassName=None):
    """
    get a list of sources in the SNCosmo registry that are of the supernova type
    sntype, and also of the same subclass if subclassName is not None

    Parameters
    ----------
    sntype: string, mandatory
        string used to denote a particular supernova type.

    subclassName: string, optional, defaults to `None`
        string used to denote a subclass of the type. 

    Returns
    -------
    list of sources in the SNCosmo registry which have the same type as sntype

    """
    sources = sncosmo.registry.get_loaders_metadata(sncosmo.Source)
    sources = numpy.array(sources)
    types = [s['type'] for s in sources]
    types = numpy.array(types)
    ok = types == sntype
    if subclassName is not None:
        subclasses = [s['subclass'][2:-1] for s in sources]
        subclasses = numpy.array(subclasses)
        ok = numpy.logical_and(ok, subclasses == subclassName)
    sources = sources[ok]

    ans = []
    for s in sources:
        ans.append(sncosmo.registry.retrieve(
            sncosmo.Source, s['name'], version=s['version']))
    return ans


def registry_sources_as_models(sntype, subclassName=None):
    """
    load SNCosmo sources matching sntype as models.

    Parameters
    ----------
    sntype: string, mandatory
        string denoting the type of SN. 

    subclassName: string, optional, defaults to `None`
        sublcass of SN type

    Returns
    -------
    List of SNCosmo models with sources matching sntype
    """
    sources = registry_sources(sntype, subclassName=subclassName)
    ans = []
    for s in sources:
        ans.append(sncosmo.Model(s))
    return ans
