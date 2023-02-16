#!/usr/bin/env python
import sys
import os
import json
import dendropy
import copy
import subprocess
from chronosynth import chronogram
import opentree
from opentree import OT, annotations, taxonomy_helpers
from helpers import crosswalk_to_dict

"""
Takes custom synth output and prunes to taxa in Clements taxonony, 
relabels tips to Clements names, 
estimates dates for internal nodes, 
and writes iTOL annotation files.

Requires as inputs the mapping from OTT_ids to Clements names available at: 
https://github.com/McTavishLab/OpenTreeCLO/blob/main/taxonomy_info/OTT_eBird_combined_taxonomy_2021.tsv

And the custom synth directory, untarred.
To perfom custom synth see:
https://opentreeoflife.github.io/CustomSynthesis/

Example:
process_custom_synth_birds.py custom_synth_runs/snacktavish_aves_81461_tmpw3m8cs_b taxonomy_info/OTT_eBird_combined_taxonomy_2021.tsv 
"""


def make_node_url(source, node):
    study, tree = source.split('@')
    return "https://tree.opentreeoflife.org/curator/study/view/{}?tab=home&tree={}&node={}".format(study, tree, node)



def collapse_to_taxa_of_interest(tree, taxa_of_interest):
    """
    Collapses internal nodes if they are present in 'taxa of interest'.
    Modifies tree itself.
    """
    for node in tree:
        if not node.is_leaf():
            if node.label in taxa_of_interest:
                node.clear_child_nodes()
                node.taxon = tree.taxon_namespace.new_taxon(label=node.label)


custom_synth_dir = os.path.abspath(sys.argv[1])
taxonomy_crosswalk = sys.argv[2] 

#-----------------------Prune tree--------------------------------------------

sys.stdout.write("Reading and pruning synth tree\n")

## Generates a dictionary to get Clements names from OTT_ids
clements_name_map = crosswalk_to_dict(taxonomy_crosswalk)
taxa_in_clements = [key for key in clements_name_map] 


## Read in the custom synth tree
custom_synth = dendropy.Tree.get_from_path("{}/labelled_supertree/labelled_supertree.tre".format(custom_synth_dir),
                                            schema = "newick")


## Count the leaves in the starting tree
leaves_start = [tip.taxon.label for tip in custom_synth.leaf_node_iter()]
sys.stdout.write("Total number of tips in synth tree is {}\n".format(len(leaves_start)))


## Custom synth includes subspecies, we collapse those to species.
collapse_to_taxa_of_interest(custom_synth, taxa_in_clements)

leaves_A = [tip.taxon.label for tip in custom_synth.leaf_node_iter()]
assert 'ott3598459' in leaves_A
sys.stdout.write("Total number of tips in synth tree after collapsing subspecies is {}\n".format(len(leaves_A)))


## Prune the tree to only taxa that have Clements names
custom_synth.retain_taxa_with_labels(taxa_in_clements)
leaves_B = [tip.taxon.label for tip in custom_synth.leaf_node_iter()]
assert 'ott3598459' in leaves_B

sys.stdout.write("Total number of tips in synth tree after pruning to taxa in Clements is {}\n".format(len(leaves_B)))
custom_synth.write_to_path(dest="{}/pruned.tre".format(custom_synth_dir), schema = "newick")


#---------------------Prune to phylo only---------------------------------------------------------------------------

annot = json.load(open("{}/annotated_supertree/annotations.json".format(custom_synth_dir)))

## Assess support for each tip. 
no_phylo_info = []
for tip in leaves_B:
    if tip in annot['nodes'].keys():
        if 'terminal' in annot['nodes'][tip].keys():
            if list(annot['nodes'][tip]['terminal'].keys()) == ['ot_2019@tree7']:
                no_phylo_info.append(tip)
        elif 'supported_by' in annot['nodes'][tip].keys():
            if list(annot['nodes'][tip]['supported_by'].keys()) == ['ot_2019@tree7']:
                no_phylo_info.append(tip)
        else:
            #print(tip)
            #print(annot['nodes'][tip])
            pass
    else:
        no_phylo_info.append(tip)

no_phylo_fi = open("{}/tips_without_phylo.txt".format(custom_synth_dir), 'w')
for tip in no_phylo_info:
    no_phylo_fi.write(clements_name_map[tip]+'\n')

no_phylo_fi.close()

sys.stdout.write("{} tips have no phylogenetic information informing their placement.\n".format(len(no_phylo_info)))

phylo_tips_only = copy.deepcopy(custom_synth)

print("pruning to phylo only tree")

phylo_tips_only.prune_taxa_with_labels(no_phylo_info)

phylo_tips = [tip.taxon.label for tip in phylo_tips_only.leaf_node_iter()]
sys.stdout.write("Total number of tips in synth tree after pruning to taxa in phylogenies is {}\n".format(len(phylo_tips)))

custom_synth.write_to_path(dest="{}/phylo_only.tre".format(custom_synth_dir), schema = "newick")

#--------------------------Write Annotations -------------------

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


# Walk through the annotations, 
# and assess how many nodes each study or tree is contributing to
study_node_count = {}
tree_node_count = {}
for node in node_support_annotation:
    for source in node_support_annotation[node]:
        study_id = source.split('@')[0]
        if study_id not in study_node_count:
            study_node_count[study_id] = 1
        else:
            study_node_count[study_id] += 1
        if source not in tree_node_count:
            tree_node_count[source] = 1
        else:
            tree_node_count[source] += 1




study_cite_file = open("{}/citation_node_counts.tsv".format(custom_synth_dir), "w")

for study_id in study_node_count:
    cites = OT.get_citations([study_id]).replace('\n','\t')
    study_cite_file.write("{}\t{}\t{}\n".format(study_id, study_node_count[study_id], cites))

study_cite_file.close()


## Map how branches across the tree support or conlfict with specific input trees.


clem_conf = 0
clem_supp = 0
jetz_conf = 0
jetz_supp = 0
all_conf = 0
phylo_supp = 0
only_clem_supp = 0

tax_conf_nodes = []

node_annotations = annotations.generate_custom_synth_node_annotation(custom_synth, custom_synth_dir)
jetz_annotations = annotations.generate_custom_synth_source_traversal(custom_synth, custom_synth_dir, "ot_809@tree2")
clem_annotations = annotations.generate_custom_synth_source_traversal(custom_synth, custom_synth_dir, "ot_2019@tree7")



for node in node_annotations:
    if jetz_annotations[node].get('support') >= 1:
        jetz_supp += 1
    if jetz_annotations[node].get('conflict') >= 1:
        jetz_conf += 1
    if clem_annotations[node].get('support') >= 1:
        clem_supp += 1
    if clem_annotations[node].get('conflict') >= 1:
        clem_conf += 1
        tax_conf_nodes.append(node)


# Label the ott_id tips to clements labels
annotations.write_itol_relabel(clements_name_map, filename="{}/ottlabel.txt".format(custom_synth_dir))


annotations.write_itol_conflict(node_annotations, filename="{}/conflict.txt".format(custom_synth_dir))
annotations.write_itol_support(node_annotations, filename="{}/support.txt".format(custom_synth_dir), param="support")


annotations.write_itol_conflict(jetz_annotations, filename="{}/jetz_conflict.txt".format(custom_synth_dir), max_conflict=1)
annotations.write_itol_support(jetz_annotations, filename="{}/jetz_support.txt".format(custom_synth_dir), param="support", max_support = 1)

annotations.write_itol_conflict(clem_annotations, filename="{}/clem_conflict.txt".format(custom_synth_dir), max_conflict=1)
annotations.write_itol_support(clem_annotations, filename="{}/clem_support.txt".format(custom_synth_dir), param="support", max_support = 1)



studies_per_tip = open("{}/studies_per_tip.txt".format(custom_synth_dir), 'w')
for tip in custom_synth.leaf_node_iter():
    studies_per_tip.write(clements_name_map[tip.taxon.label] + ":")
    studies_per_tip.write(", ".join(annot['nodes'].get(tip.taxon.label, {}).get('supported_by', ' ')))
    studies_per_tip.write(", ".join(annot['nodes'].get(tip.taxon.label, {}).get('terminal', ' ')))
    studies_per_tip.write('\n')

studies_per_tip.close()


for node in node_support_annotation:
    supp = node_support_annotation.get(node, {})
    if len(set(node_support_annotation[node].keys()).difference(set(['ot_2019@tree7']))) >=1:
        phylo_supp += 1
    if set(node_support_annotation[node].keys()) == (set(['ot_2019@tree7'])):
        only_clem_supp += 1



#----------- Summary information ------------------------------------------
CLO_spp_num = 10824
print("{lt} trees from {ls} published studies contributed information to this tree".format(lt=len(tree_node_count), ls=len(study_node_count)))

summary_statement = """This tree contains {l} leaves and {i} internal nodes.
                        Of those nodes, {asl} are strictly supported 
                        by at least 1 input phylogeny. The rest ({tsl}) 
                        are placed by taxonomy.""".format(l=len(leaves_B),
                                                                i=len(all_nodes),
                                                                asl=phylo_supp,
                                                                tsl=len(all_nodes)-phylo_supp)
print(summary_statement)

print("""This comprises {per:.2} of the {cls} species 
                        in the 2021 clements taxonomy""".format(per = len(leaves_B)/CLO_spp_num,
                                                               cls=CLO_spp_num))

print("""{cc} nodes in the tree disagree with current CLO taxonomy""".format(cc=clem_conf))
                                                                    


#---------------------------- Dating ---------------------------


## Rough estimate on dates

custom_str = chronogram.conflict_tree_str(custom_synth)


dates_dir = "{}/dates".format(custom_synth_dir)
if not os.path.exists(dates_dir):
    os.mkdir(dates_dir)

select_dates_file = "{}/dates/custom_node_ages.json".format(custom_synth_dir)


selected_bird_chrono = ['ot_2018@tree8', 'ot_2013@tree8']
if os.path.exists(select_dates_file):
    select_dates = json.load(open(select_dates_file))
else:
    select_dates = chronogram.combine_ages_from_sources(selected_bird_chrono,#list(bird_chrono)
                                                        json_out = selected_bird_chrono,
                                                        compare_to = custom_str)


dates_cite_file = open("{}/dates/select_dates_citations.txt".format(custom_synth_dir), "w")
cites = OT.get_citations(selected_bird_chrono)
dates_cite_file.write(cites)
dates_cite_file.close()

print("checking mrca")
root_node = OT.synth_mrca(node_ids=leaves_B).response_dict['mrca']['node_id']

max_age_est = 130
#https://academic.oup.com/sysbio/article/63/3/442/1649269

print("dating phylo only tree, custom, mean")

treesfile, sources = chronogram.date_tree(phylo_tips_only,
                                          select_dates,
                                          root_node,
                                          max_age_est,
                                          method='bladj',
                                          output_dir="{}/dates_select_phylo_only".format(custom_synth_dir),
                                          phylo_only=False,
                                          reps=1,
                                          select = "mean")


dated_phylo= dendropy.Tree.get_from_path(treesfile, schema = "newick")

dated_phylo.write(path="{}/phylo_only_select_dates_mean_ott_labels.tre".format(custom_synth_dir), schema="newick")

for tax in dated_phylo.taxon_namespace:
    tax.label = clements_name_map[tax.label]

dated_phylo.write(path="{}/phylo_only_select_dates_mean_clements_labels.tre".format(custom_synth_dir), schema="newick")

print("""date information for {ld} nodes in the tree
         was summarized from {lds} published studies""".format(ld=len(select_dates['node_ages']),
                                                                lds=len(selected_bird_chrono)))