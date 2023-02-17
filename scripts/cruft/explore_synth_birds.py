import dendropy
import copy
from chronosynth import chronogram
from opentree import OT, annotations, taxonomy_helpers


infi = open("../taxonomy_info/OTT_eBird_combined_taxonomy_2021.tsv")
header = infi.readline().split('\t')
ebird_ott_ids = set()



for lin in infi:
    lii = lin.split('\t')
    ebird_ott_ids.add('ott'+lii[9])


infi.close()



## First find out if it is a combined synth tip taxon
aves_all_resp = OT.synth_subtree(node_id = 'ott81461', label_format = 'id')

aves_all_tree = dendropy.Tree.get_from_string(aves_all_resp.response_dict['newick'], schema = 'newick')

all_aves_synth_tip = [tip.taxon.label for tip in aves_all_tree.leaf_node_iter()]


aves_phylo_tree =  chronogram.prune_to_phylo_only(aves_all_tree)

aves_phylo_synth_tip = [tip.taxon.label for tip in aves_phylo_tree.leaf_node_iter()]



curated_tips_no_ebird=set(aves_phylo_synth_tip).difference(ebird_ott_ids)


ott_not_in_ebird_fi = open("ott_curated_tips_not_in_ebird.csv",'w')

for tip in curated_tips_no_ebird:
    ottid=tip.strip('ott')
    info = OT.taxon_info(ottid)
    tip_name = (info.response_dict['name'])
    res = OT.find_trees(ottid, search_property='ot:ottId')
    count_studies = len(res.response_dict['matched_studies'])
#    for study in res.response_dict['matched_studies']:
#        trees = study['matched_trees']
    ott_not_in_ebird_fi.write("{},{},{}\n".format(ottid, tip_name, count_studies))      

ott_not_in_ebird_fi.close()