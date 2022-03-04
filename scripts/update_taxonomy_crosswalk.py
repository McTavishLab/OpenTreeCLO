import dendropy
import copy
import csv
from opentree import OT
from chronosynth import chronogram
from opentree import OT, annotations, taxonomy_helpers


infi = open("../taxonomy_info/combined_data_2019.tsv")
header = infi.readline().split('\t')


tax_2019 = {}
no_maps = []
double_maps = {}
ott_dict = {}

for lin in infi:
    lii = lin.strip().split('\t')
    tax_2019[lii[4]]=dict(zip(header, lii))
    if lii[9]=='-' or lii[9]=='':
            no_maps.append(lii)
    else:
        if 'ott'+lii[9] in ott_dict:
            if 'ott'+lii[9] in double_maps:
                double_maps['ott'+lii[9]].append(dict(zip(header, lii)))
            else:
                double_maps['ott'+lii[9]]=[dict(zip(header, lii))]
                double_maps['ott'+lii[9]].append(ott_dict['ott'+lii[9]])
        ott_dict['ott'+lii[9]]=dict(zip(header, lii))

infi.close()


new_in_2021 ={}
infi_2021 = open("../taxonomy_info/eBird_Taxonomy_v2021.csv", encoding='utf-8-sig')
tax_header = infi_2021.readline().strip().split(',')

output_header = header.remove('parent_ottid')
output_header = header.remove('uniqname')

outfi = open("../taxonomy_info/OTT_eBird_combined_taxonomy_2021.tsv", "w")
outfi.write('\t'.join(header))

with open("../taxonomy_info/eBird_Taxonomy_v2021.csv", newline='') as csvfile:
    taxreader = csv.reader(csvfile, delimiter=',')
    next(taxreader)
    for lii in taxreader:
        taxon_info_2021 = dict(zip(tax_header, lii))
        if taxon_info_2021['SCI_NAME'] in tax_2019:
            combo_2021 = tax_2019[taxon_info_2021['SCI_NAME']]
            combo_2021['TAXON_ORDER'] = taxon_info_2021['TAXON_ORDER']
            ##Assert ottid on 3.3
        else:
            new_in_2021[taxon_info_2021['SCI_NAME']]=taxon_info_2021
            combo_2021 = taxon_info_2021
            res = OT.tnrs_match(names = [taxon_info_2021['SCI_NAME']], context_name="Birds")
            if len(res.response_dict['results'][0]['matches']) == 1:
                combo_2021['ott_id'] = str(res.response_dict['results'][0]['matches'][0]['taxon']['ott_id'])
                combo_2021['tax_sources'] = ','.join(res.response_dict['results'][0]['matches'][0]['taxon']['tax_sources'])
                combo_2021['flags'] = ','.join(res.response_dict['results'][0]['matches'][0]['taxon']['flags'])
                combo_2021['rank'] = ','.join(res.response_dict['results'][0]['matches'][0]['taxon']['rank'])
                if res.response_dict['results'][0]['matches'][0]['is_synonym']==False:
                    combo_2021['match_type'] = 'canonical_match'
                elif res.response_dict['results'][0]['matches'][0]['is_synonym']==True:
                    combo_2021['match_type'] = 'synonym_match'
        outfi.write('\t'.join([combo_2021.get(col, '-') for col in header])+'\n')

outfi.close()

doubles = open("double_match_spp.tsv",'w')
doubles.write('\t'.join(header))

for ott_id in double_maps:
    for match in double_maps[ott_id]:
        doubles.write('\t'.join([match[key] for key in header]))

doubles.close()


missing = open("no_match_spp.tsv",'w')
missing.write('\t'.join(header))
for item in no_maps:
    missing.write('\t'.join(item))

missing.close()
