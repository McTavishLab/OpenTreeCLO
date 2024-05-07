#!/usr/bin/env python
import sys
import os
import csv
import json
import dendropy
import subprocess
from opentree import OT
from chronosynth import chronogram
from helpers import crosswalk_to_dict

"""
Pass in tree with ott_id labels
"""




custom_synth_dir = os.path.abspath(sys.argv[1])
input_tree_file = sys.argv[2] 
taxonomy_crosswalk = sys.argv[3] 
label = sys.argv[4]

## Generates a dictionary to get Clements names from OTT_ids
clements_name_map = crosswalk_to_dict(taxonomy_crosswalk)

custom_synth = dendropy.Tree.get(path=input_tree_file, schema = 'newick')
## Estimate dates

custom_str = chronogram.conflict_tree_str(custom_synth)

#now dates



dates_dir = "{}/dates".format(custom_synth_dir)
if not os.path.exists(dates_dir):
    os.mkdir(dates_dir)

all_dates_file = "{}/dates/all_node_ages.json".format(custom_synth_dir)
select_dates_file = "{}/dates/custom_node_ages.json".format(custom_synth_dir)
    
if os.path.exists(all_dates_file):
    all_dates = json.load(open(all_dates_file))
else:
    all_chronograms = chronogram.find_trees(search_property="ot:branchLengthMode", value="ot:time")
    all_birds = chronogram.find_trees(search_property="ot:ottId", value="81461")
    bird_chrono = set(all_chronograms).intersection(set(all_birds))
    problematic_chrono = set(["ot_2008@tree5", "pg_2829@tree6579"])
    bird_chrono = bird_chrono - problematic_chrono
    all_dates = chronogram.combine_ages_from_sources(list(bird_chrono),
                                                 json_out = all_dates_file,
                                                 compare_to = custom_str)


if os.path.exists(select_dates_file):
    select_dates = json.load(open(select_dates_file))
else:
    selected_bird_chrono = ['ot_2018@tree8', 'ot_2018@tree9', 'ot_2013@tree8']
    select_dates = chronogram.combine_ages_from_sources(selected_bird_chrono,#list(bird_chrono)
                                                        json_out = select_dates_file,
                                                        compare_to = custom_str)





if os.path.exists(select_dates_file):
    select_dates = json.load(open(select_dates_file))
else:
    selected_bird_chrono = ['ot_2018@tree9']
    select_dates = chronogram.combine_ages_from_sources(selected_bird_chrono,#list(bird_chrono)
                                                        json_out = select_dates_file,
                                                        compare_to = custom_str)



print("checking mrca")
leaves = [tip.taxon.label for tip in custom_synth.leaf_node_iter()]
root_node = OT.synth_mrca(node_ids=leaves).response_dict['mrca']['node_id']

max_age_est = 130



print("dating full tree, all dates, mean")

for i in range(10):
    treesfile, sources = chronogram.date_tree(custom_synth,
                                          all_dates,
                                          root_node,
                                          max_age_est,
                                          method='bladj',
                                          output_dir="{}/dates/dates_all_rand_sample{}".format(custom_synth_dir, label),
                                          select = "random",
                                          reps = 5
                                          )



treesfile, sources = chronogram.date_tree(custom_synth,
                                          all_dates,
                                          root_node,
                                          max_age_est,
                                          method='bladj',
                                          output_dir="{}/dates/dates_all_{}".format(custom_synth_dir, label),
                                          select = "mean",
                                          reps = 1
                                          )

dated_phylo = dendropy.Tree.get_from_path(treesfile, schema = "newick")
dated_phylo.write(path="{}/dates/dated_mean_all_dates_ott_labels_{}.tre".format(custom_synth_dir, label), schema="newick")

for tax in dated_phylo.taxon_namespace:
    tax.label = clements_name_map[tax.label]

dated_phylo.write(path="{}/dates/dated_mean_all_dates_clements_labels_{}.tre".format(custom_synth_dir, label), schema="newick")


#------------------------------

print("dating full tree, select dates, mean")

treesfile, sources = chronogram.date_tree(custom_synth,
                                          all_dates,
                                          root_node,
                                          max_age_est,
                                          method='bladj',
                                          output_dir="{}/dates/dates_select_{}".format(custom_synth_dir, label),
                                          select = "mean")



dated_phylo = dendropy.Tree.get_from_path(treesfile, schema = "newick")

dated_phylo.write(path="{}/dates/dated_mean_select_dates_ott_labels_{}.tre".format(custom_synth_dir, label), schema="newick")

for tax in dated_phylo.taxon_namespace:
    tax.label = clements_name_map[tax.label]

dated_phylo.write(path="{}/dates/dated_mean_select_dates_clements_labels_{}.tre".format(custom_synth_dir, label), schema="newick")

#-----------------------------------------------------------------------------------------------
#print("dating whole tree, full")

#treesfile, sources = chronogram.date_tree(custom_synth,
#                                           all_dates,
#                                           root_node,
#                                           max_age_est,
#                                           method='bladj',
#                                           output_dir="{}/dates/dates_all_multi".format(custom_synth_dir),
#                                           reps=10,
#                                           select = "random")

# treelist = dendropy.TreeList.get(path=treesfile, schema='newick')

# summaryfilepath = '{}/dates/dated_all_summary_ott_labels.tre'.format(custom_synth_dir)
# subprocess.call(['sumtrees.py', '-f0.95', '-v0.01', '-o{}'.format(summaryfilepath), '-emean-age', treesfile]) 

# dated = dendropy.Tree.get_from_path(summaryfilepath, schema = "nexus")

# for tax in dated.taxon_namespace:
#     tax.label = clements_name_map[tax.label]


# dated.write(path="{}/dates/dated_all_summary_clements_labels.tre".format(custom_synth_dir), schema="newick")

#--------------------Generating citations -------------------------------------

node_date_count = 0
matched_date_studies = set()
for node in all_dates['node_ages']:
    node_date_count += 1
    for source in all_dates['node_ages'][node]:
        matched_date_studies.add(source['source_id'].split('@')[0])


dates_cite_file = open("{}/dates/full_dates_citations.txt".format(custom_synth_dir), "w")
cites = OT.get_citations(matched_date_studies)
dates_cite_file.write(cites)
dates_cite_file.close()


#--------------------------------------------------------------------------------------------------

print("""date information for {ld} nodes in the tree
         was summarized from {lds} published studies""".format(ld=node_date_count,
                                                                lds=len(matched_date_studies)))
