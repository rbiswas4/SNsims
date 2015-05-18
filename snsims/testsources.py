import sources

print sources.registry_types()

print sources.registry_sources(sntype='SN IIP')

ms = sources.registry_sources_as_models(sntype='SN IIP')

for m in ms:
    print '================================'
    print m
