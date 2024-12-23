import pandas as pd
import pyarrow.feather as feather

infname = '../data/cellosaurus.txt'

inf = open(infname, 'r')

# find end of header
for line in inf:
    if line[0] == '_':
        break;

# parse records
sep = '|'
objs = []
for line in inf:
    key = line[:2]
    # terminator
    if key == '//':
        objs.append(obj)
        continue
    # initiator
    if key == 'ID':
        obj = {}
    value = line[5:].rstrip()
    if key == 'DR':
        value = value.replace('; ', '=')
    if key in obj:
        obj[key] += sep + value
    else:
        obj[key] = value

inf.close()

df = pd.DataFrame.from_records(objs)

# ID         Identifier (cell line name)     Once; starts an entry
# AC         Accession (CVCL_xxxx)           Once
# AS         Secondary accession number(s)   Optional; once
# SY         Synonyms                        Optional; once
# DR         Cross-references                Optional; once or more
# RX         References identifiers          Optional: once or more
# WW         Web pages                       Optional; once or more
# CC         Comments                        Optional; once or more
# ST         STR profile data                Optional; twice or more
# DI         Diseases                        Optional; once or more
# OX         Species of origin               Once or more
# HI         Hierarchy                       Optional; once or more
# OI         Originate from same individual  Optional; once or more
# SX         Sex of cell                     Optional; once
# AG         Age of donor at sampling        Optional; once
# CA         Category                        Once
# DT         Date (entry history)            Once

df.rename(
    columns = {
        'ID': 'id',
        'AC': 'accession',
        'AS': 'accession2',
        'SY': 'synonyms',
        'DR': 'xref',
        'RX': 'ref',
        'WW': 'web',
        'CC': 'comment',
        'ST': 'str_profile',
        'DI': 'disease',
        'OX': 'species',
        'HI': 'hierarchy',
        'OI': 'origin',
        'SX': 'sex',
        'AG': 'age_donor',
        'CA': 'category',
        'DT': 'date'
    },
    inplace = True
)


feather.write_feather(df, 'cellosaurus.feather')

df.to_csv('cellosaurus.csv', index=False)
