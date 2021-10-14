import os
import sys
import create_variant_catalog

EHDN_input = sys.argv[1]
reference = sys.argv[2]
variant_catalog = sys.argv[3]
unmatched = sys.argv[4]
excluded = sys.argv[5]

create_variant_catalog.transform_format_sarah(EHDN_input, reference, variant_catalog, unmatched, excluded)
