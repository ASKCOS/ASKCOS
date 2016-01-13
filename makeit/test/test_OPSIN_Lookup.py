from makeit.utils.OPSIN_Lookup import OPSIN_Lookup
from time import sleep

# Create instance of wrapper class
Looker = OPSIN_Lookup()
print 'initialized OPSIN jar'

# Test simple case
assert(Looker.lookup('methane') == 'C')
print '...tested methane, ok'
assert(Looker.lookup('propane') == 'CCC')
print '...tested propane, ok'

# Close subprocess
del(Looker)
print 'closed OPSIN jar'