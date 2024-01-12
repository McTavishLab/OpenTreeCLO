import dendropy
import sys
import os
import csv
import copy
import json
import random
from opentree import OT, annotations
from helpers import crosswalk_to_dict

def make_node_url(source, node):
    study, tree = source.split('@')
    return "https://tree.opentreeoflife.org/curator/study/view/{}?tab=home&tree={}&node={}".format(study, tree, node)





custom_synth_dir = sys.argv[1] 
input_tree_file = sys.argv[2] 
custom_synth = dendropy.Tree.get(path=input_tree_file, schema="newick")
name_map = crosswalk_to_dict("/home/ejmctavish/projects/otapi/OpenTreeCLO/taxonomy_info/OTT_crosswalk_2021.csv")

fam_name_map = crosswalk_to_dict("/home/ejmctavish/projects/otapi/OpenTreeCLO/taxonomy_info/OTT_crosswalk_2021.csv", alt_name="FAMILY")
## function to get back walk from ott o clo

for name in fam_name_map:
    fam_name_map[name] = fam_name_map[name].split()[0]

families = set(fam_name_map.values())

families_according_to_CLO = {fam:[] for fam in families}
for ott in fam_name_map:
    families_according_to_CLO[fam_name_map[ott]].append(ott)

tips = set([leaf.taxon.label for leaf in custom_synth.leaf_node_iter()])


mono_fams = {}
non_mono_fams = []
all_fams = {}

small_fams = []

if not os.path.exists("{}/family_trees".format(custom_synth_dir)):
    os.mkdir("{}/family_trees".format(custom_synth_dir))

annot = json.load(open("{}/annotated_supertree/annotations.json".format(custom_synth_dir)))


# Walk through the nodes, and summaraize annotations
uncontested_taxa = []
node_support_annotation = {}
all_nodes = []

for node in custom_synth:
    label = None
    if not node.is_leaf():
        label = node.label
        all_nodes.append(label)
        assert label
        if label == 'ott81461':
            pass
        else:
            if label in annot['nodes'] and 'supported_by' in annot['nodes'][label].keys():
                    strict_support = annot['nodes'][label]['supported_by']
                    node_support_annotation[label] = strict_support
            else:
                assert label.startswith('ott')
                uncontested_taxa.append(label)


for family in families:
    fam_tips_in_this_tree = tips.intersection(set(families_according_to_CLO[family]))
    fam_tree = copy.deepcopy(custom_synth)
    mrca = custom_synth.mrca(taxon_labels=fam_tips_in_this_tree)
    tips_from_mrca = [leaf.taxon.label for leaf in mrca.leaf_iter()]
    fam_tree.retain_taxa_with_labels(tips_from_mrca)
    fam_tree.write(path ="{}/family_trees/{}.tre".format(custom_synth_dir, family), schema = "newick")
    fam_sources = set()
    study_node_count = {}
    tree_node_count = {}
    for node in fam_tree:
        label = node.label
        if label == 'ott81461':
            pass
        if label in node_support_annotation:
            for source in node_support_annotation[label]:
                study_id = source.split('@')[0]
                if study_id not in study_node_count:
                    study_node_count[study_id] = 1
                else:
                    study_node_count[study_id] += 1
                if source not in tree_node_count:
                    tree_node_count[source] = 1
                else:
                        tree_node_count[source] += 1
    study_cite_file = open("{}/family_trees/{}_citation_node_counts.tsv".format(custom_synth_dir, family), "w")
    for study_id in study_node_count:
        if study_id == 'ot_2019':
            pass
        else:
            cites = OT.get_citations([study_id]).replace('\n','\t')
            study_cite_file.write("{}\t{}\t{}\n".format(study_id, study_node_count[study_id], cites))
    study_cite_file.close()
    intruders = set(tips_from_mrca).difference(set(fam_tips_in_this_tree))
    print("fam {} intruders".format(family))
    print(",".join(intruders))
    chars = '0123456789ABCDEF'
    color = '#'+''.join(random.sample(chars,6))
    label = None
    label = mrca.label
    if not label:
        label = mrca.taxon.label
    if label:
        all_fams[label] = family
        if len(intruders) >= 1:
            non_mono_fams.append(family)
        else:
            if len(fam_tips_in_this_tree) >= 10:
                mono_fams[label]=family
    else:
        small_fams.append(family)


annotations.write_itol_clades(mono_fams, filename="{}/mono_fams_10.txt".format(custom_synth_dir))


annotations.write_itol_clades(all_fams, filename="{}/all_fams.txt".format(custom_synth_dir))

non_mon = open("{}/non_mono_fams.txt".format(custom_synth_dir), 'w')


non_mono_fams.sort()


i = 0
for family in non_mono_fams:
    i += 1
    non_mon.write(family + '\n')
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