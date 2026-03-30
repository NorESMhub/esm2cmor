#/usr/bin/env python

import json

fname = '../recipes/template/variable_mapping_NorESM3_to_CMIP7.json'
with open(fname,'r') as file:
    data_dict = json.load(file)

sorted_data_dict = dict()
sorted_data_dict['Header'] = data_dict['Header']
sorted_data_dict['variable_entry'] = {k: v for k, v in sorted(data_dict['variable_entry'].items(), key=lambda item: item[0])}

print(sorted_data_dict)
#
with open(fname, 'w') as f:
    json.dump(sorted_data_dict, f, indent=4) #

