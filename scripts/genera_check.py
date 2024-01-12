import dendropy
import sys
import os
import csv
import copy
import json
import random
from opentree import OT, annotations
from helpers import crosswalk_to_dict


custom_synth_dir = sys.argv[1] 
input_tree_file = custom_synth_dir+"/phylo_only.tre"
taxonomy_crosswalk = sys.argv[2] 

custom_synth = dendropy.Tree.get(path=input_tree_file, schema="newick")
name_map = crosswalk_to_dict(taxonomy_crosswalk)
opentree_name_map = crosswalk_to_dict(taxonomy_crosswalk, alt_name="name")


genera_dict = {}
for ott_id in name_map:
    genus = name_map[ott_id].split()[0]
    if genus in genera_dict:
        genera_dict[genus].append(ott_id)
    else:
        genera_dict[genus] = [ott_id]

to_constrain = ['Gracupica','Nesoenas','Heteromyias','Riparia','Ninox','Poecilodryas','Spelaeornis','Nesocharis','Rallicula','Caligavis','Ketupa']
for genus in to_constrain:
     print("('{}');".format("','".join([opentree_name_map[ott_id] for ott_id in genera_dict[genus]])))


tips = set([leaf.taxon.label for leaf in custom_synth.leaf_node_iter()])


mono_genera = {}
non_mono_genera = []
all_genera = {}

annot = json.load(open("{}/annotated_supertree/annotations.json".format(custom_synth_dir)))



for genus in genera_dict:
    if len(genera_dict[genus]) == 1:
        continue
    genus_tips_in_this_tree = tips.intersection(set(genera_dict[genus]))
    genus_tree = copy.deepcopy(custom_synth)
    if len(genus_tips_in_this_tree) >= 2:
        mrca = custom_synth.mrca(taxon_labels=genus_tips_in_this_tree)
        tips_from_mrca = [leaf.taxon.label for leaf in mrca.leaf_iter()]
        intruders = set(tips_from_mrca).difference(set(genus_tips_in_this_tree))
        if len(intruders) == 0:
            print("genus {} is OK\n".format(genus))
        else:
            print("genus {} includes {}\n".format(genus, ", ".join([name_map[ott_id] for ott_id in intruders])))
            non_mono_genera.append(genus)


sources = set()
for node in annot['nodes']:
    for source in annot['nodes'][node].get('supported_by',{}):
        sources.add(source)

genus_info_sheet = open("{}/genus_conflicts.txt".format(custom_synth_dir), "w")

genus_conflicts = {genus:[] for genus in non_mono_genera}
conflict_tree_dir = "{}/genus_conflict_trees".format(custom_synth_dir)

if not os.path.exists(conflict_tree_dir):
    os.mkdir(conflict_tree_dir)


for genus in non_mono_genera:
     for source in sources:
            #print(annot['nodes'][node])
            #print(make_node_url(source, annot['nodes'][node]['supported_by'][source]))
            conf_tree = dendropy.Tree.get(path='{}/cleaned_phylo/tree_{}.tre'.format(custom_synth_dir,source),
                          schema="newick",
                          preserve_underscores=True)
            conf_tips_all = [leaf.taxon.label for leaf in conf_tree.leaf_node_iter()]
            clo_conf_tips = [label for label in conf_tips_all if label.split('_')[-1] in name_map]
            if clo_conf_tips == []:
                print("Source {} has no tips in CLO\n".format(source))
                continue
            conf_tree.retain_taxa_with_labels(clo_conf_tips)
            genus_tips_in_this_tree = [label for label in conf_tips_all if label.split('_')[-1] in genera_dict[genus]]
            if len(genus_tips_in_this_tree) >= 2:
                print(genus)
                print(len(genus_tips_in_this_tree))
                genus_info_sheet.write("\nTips from {} genus in this tree are:\n".format(genus))
                genus_info_sheet.write(",".join(genus_tips_in_this_tree))
                mrca = conf_tree.mrca(taxon_labels=genus_tips_in_this_tree)
                tips_from_mrca = [leaf.taxon.label for leaf in mrca.leaf_iter()]
                intruders = set(tips_from_mrca).difference(set(genus_tips_in_this_tree))
                if len(intruders) >= 1:
                    conf_tree.retain_taxa_with_labels(tips_from_mrca)
                    for tip in conf_tree.leaf_node_iter():
                        ott_id =  tip.taxon.label.split('_')[-1]
                        clo_name = name_map[ott_id]
                        tip.taxon.label =  clo_name +"__ott:" + tip.taxon.label
                    assert len(tips_from_mrca) == len(genus_tips_in_this_tree) + len(intruders)
                    genus_conflicts[genus].append(source)
                    conf_tree.write(path = "{}/{}_{}.tre".format(conflict_tree_dir, genus, source.replace('@', '_')), schema="newick")
                else:
                    genus_info_sheet.write("\nTree {} does not conflict with genus {}\n Contains:".format(source, genus))
                    genus_info_sheet.write(",".join(genus_tips_in_this_tree))
       #     else:
       #         genus_info_sheet.write("\nTree {} contains One or zero tips from family {}\n".format(source, genus))
       #         genus_info_sheet.write(",".join(genus_tips_in_this_tree))



summary_source_page = open("{}/genus_info.txt".format(custom_synth_dir), "w")
for genus in set(non_mono_genera):
    summary_source_page.write("\n\n****************\n"+genus + " is non-monophyletic in \n")
    if len(set(genus_conflicts[genus])) >= 1:
        cites= OT.get_citations(list(set(genus_conflicts[genus])))
        summary_source_page.write(cites)
    else:
        summary_source_page.write("ERROR: No individual study supports non-monophyly")



summary_source_page.close()
"""



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
        genus_tree.retain_taxa_with_labels(tips_from_mrca)
        genus_sources = set()
        study_node_count = {}
        tree_node_count = {}
        for node in genus_tree:
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
"""