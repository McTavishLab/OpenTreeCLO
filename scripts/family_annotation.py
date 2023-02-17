import dendropy
import sys
import json
import random
from opentree import OT
from helpers import crosswalk_to_dict

def make_node_url(source, node):
    study, tree = source.split('@')
    return "https://tree.opentreeoflife.org/curator/study/view/{}?tab=home&tree={}&node={}".format(study, tree, node)


input_tree_file = sys.argv[1] 
custom_synth = dendropy.Tree.get(path=input_tree_file, schema="newick")
name_map = crosswalk_to_dict("/home/ejmctavish/projects/otapi/OpenTreeCLO/taxonomy_info/OTT_eBird_combined_taxonomy_2021.tsv")

fam_name_map = crosswalk_to_dict("/home/ejmctavish/projects/otapi/OpenTreeCLO/taxonomy_info/OTT_eBird_combined_taxonomy_2021.tsv", alt_name="FAMILY")
## function to get back walk from ott o clo

for name in fam_name_map:
    fam_name_map[name] = fam_name_map[name].split()[0]

families = set(fam_name_map.values())

families_according_to_CLO = {fam:[] for fam in families}
for ott in fam_name_map:
    families_according_to_CLO[fam_name_map[ott]].append(ott)

tips = set([leaf.taxon.label for leaf in custom_synth.leaf_node_iter()])

non_mono_fams = []

small_fams = []

for family in families:
    fam_tips_in_this_tree = tips.intersection(set(families_according_to_CLO[family]))
    mrca = custom_synth.mrca(taxon_labels=fam_tips_in_this_tree)
    tips_from_mrca = [leaf.taxon.label for leaf in mrca.leaf_iter()]
    intruders = set(tips_from_mrca).difference(set(fam_tips_in_this_tree))
    chars = '0123456789ABCDEF'
    color = '#'+''.join(random.sample(chars,6))
    label = None
    label = mrca.label
    if not label:
        label = mrca.taxon.label
    if label:
        if len(intruders) >= 1:
            non_mono_fams.append(family)
        else:
            if len(fam_tips_in_this_tree) >= 50:
                 print("{n},{n},#ffffff,{c},#000000,solid,0,{f},#000000,4,bold".format(n=label, c=color, f=family))
    else:
        small_fams.append(family)

non_mono_fams.sort()

i = 0
for family in non_mono_fams:
    i += 1
    fam_tips_in_this_tree = tips.intersection(set(families_according_to_CLO[family]))
    for tip in fam_tips_in_this_tree:
        print("{t} {i} #{f}".format(t=tip, i=i*100, f=family))


i = 0
for family in non_mono_fams:
    i += 1
    fam_tips_in_this_tree = tips.intersection(set(families_according_to_CLO[family]))
    if len(fam_tips_in_this_tree) < 20:
        for tip in fam_tips_in_this_tree:
            print("{t} {i} #{f}".format(t=tip, i=i, f=family))


i = 0
for family in small_fams:
    i += 1
    fam_tips_in_this_tree = tips.intersection(set(families_according_to_CLO[family]))
    for tip in fam_tips_in_this_tree:
        print("{t} {i} #{f}".format(t=tip, i=i, f=family))