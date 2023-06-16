#!/usr/bin/env python
import dendropy
from helpers import crosswalk_to_dict

from opentree import OT

taxonomy_crosswalk = "/home/ejmctavish/projects/otapi/OpenTreeCLO/taxonomy_info/OTT_eBird_combined_taxonomy_2021.tsv"
jetz_taxonomy = "/home/ejmctavish/projects/otapi/OpenTreeCLO/taxonomy_info/Jetz-2012-master_taxonomy.csv"


## Generates a dictionary to get Clements names from OTT_ids
clements_name_map = crosswalk_to_dict(taxonomy_crosswalk)
ott_name_map = crosswalk_to_dict(taxonomy_crosswalk, alt_name='name')
taxa_in_clements = [key for key in clements_name_map] 


output = OT.get_tree("ot_809", "tree1", tree_format="newick", label_format="ot:ottId")

jetz_tree_full = dendropy.Tree.get(data=output.tree.decode(),
                                   schema = "newick",
                                   suppress_internal_node_taxa=True,
                                   suppress_leaf_node_taxa=True)

jetz_tree_cleaned = dendropy.Tree.get(path="/home/ejmctavish/projects/otapi/OpenTreeCLO/custom_synth_runs/aves_0.1/cleaned_phylo/tree_ot_809@tree2.tre",
                                      schema = "newick",
                                      suppress_internal_node_taxa=True,
                                      suppress_leaf_node_taxa=True)


jetz_leaves_full = [tip.label for tip in jetz_tree_full.leaf_node_iter()]


jetz_leaves_cleaned = [tip.label.split('_')[-1] for tip in jetz_tree_cleaned.leaf_node_iter()]


mia = [leaf for leaf in jetz_leaves_full if "ott" + str(leaf) not in jetz_leaves_cleaned]

forwards = {}

still_mia = []
for leaf in mia:
    if leaf.isdigit():
        tax_info = OT.taxon_info(leaf).response_dict
        if tax_info['ott_id'] != int(leaf):
            forwards[leaf] = tax_info['ott_id']
            assert "ott"+str(tax_info['ott_id']) in jetz_leaves_cleaned
        else:
            still_mia.append(leaf)
    else:
        print(leaf)




mapped_to_subspp = {}
##Some are type sub spp, where ther is phylo info for the spp.
for tip in ott_name_map:
    if len(ott_name_map[tip].split())==3:
        if tip == "ott6520478":
            mapped_to_subspp[tip] = 'ott515130'
        else:
           parent = OT.taxon_info(tip, include_lineage=True).response_dict['lineage'][0]['ott_id']
           mapped_to_subspp[tip] = "ott" + str(parent)
        

jetz_not_in_ebird = []
for tip in jetz_leaves:
    if tip.startswith("ott*"):
        print(tip.strip('ott'))
        jetz_not_in_ebird.append(tip)
    if
    elif tip not in clements_name_map:
        ott_name = OT.taxon_info(tip).response_dict['name']
        print(ott_name, tip)
        jetz_not_in_ebird.append(tip)

