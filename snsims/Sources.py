import sncosmo
import numpy
import sys

def registry_types():
	sources = sncosmo.registry.get_loaders_metadata(sncosmo.Source)
	types = [s['type'] for s in sources]
	return set(types)

def registry_sources(type, subclassName=None):
	sources = sncosmo.registry.get_loaders_metadata(sncosmo.Source)
	sources = numpy.array(sources)
	types = [s['type'] for s in sources]
	types=numpy.array(types)
	ok = types == type
	if subclassName is not None:
		subclasses = [s['subclass'][2:-1] for s in sources]
		subclasses = numpy.array(subclasses)
		ok = numpy.logical_and(ok, subclasses == subclassName)
	sources = sources[ok]

	ans=[]
	for s in sources:
		ans.append(sncosmo.registry.retrieve(sncosmo.Source,s['name'], version=s['version']))
	return ans

def registry_sources_as_models(type, subclassName=None):
	sources = registry_sources(type, subclassName=subclassName)
	ans = []
	for s in sources:
		ans.append(sncosmo.Model(s))
	return ans