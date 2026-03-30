#/usr/bin/env python

import json

fname = '../tables/CMIP7_ocean.json'

with open(fname,'r') as file:
    data = json.load(file)
    unique_dimensions = {item for key in data['variable_entry']
                        for item in data['variable_entry'][key]['dimensions']}
print(unique_dimensions)

