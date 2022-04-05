import sys
import os
import dendropy
import copy
from chronosynth import chronogram
from opentree import OT, annotations, taxonomy_helpers
from helpers import crosswalk_to_dict



custom_synth_dir = os.path.abspath(sys.argv[1])
taxonomy_crosswalk = sys.argv[2] 

# "/home/ejmctavish/projects/otapi/OpenTreeCLO/custom_synth_runs/snacktavish_aves_81461_tmpvu7e5t2x"

#$ curl -X POST https://ot38.opentreeoflife.org/v3/tree_of_life/build_tree -d '{"input_collection":"snacktavish/woodpeckers", "root_id": "ott1020138"}'

#{"opentree_home": "/home/otcetera/custom_synth_repos", "ott_dir": "/home/otcetera/custom_synth_repos/ott3.2-extinct-flagged", "root_ott_id": "1020138", "synth_id": "snacktavish_Woodpeckers_1020138_tmpl_hmudtg", "collections": "snacktavish/Woodpeckers", "cleaning_flags": "major_rank_conflict,major_rank_conflict_inherited,environmental,viral,barren,not_otu,hidden,was_container,inconsistent,hybrid,merged", "addi ~/Desktop/grants/ebirb$ B

#$
#!curl -X GET https://ot38.opentreeoflife.org/v3/tree_of_life/custom_built_tree/snacktavish_woodpeckers_1020138_tmp378jz32z.tar.gz --output woodpecker_synth.tar.gz 

#!tar -xzvf aves_all_plus_tax.tar.gz

name_map = crosswalk_to_dict(taxonomy_crosswalk)
## function to get back walk from ott o clo



current = dendropy.Tree.get_from_path("{}/labelled_supertree/labelled_supertree.tre".format(custom_synth_dir),
                                       schema = "newick")

leaves_start = [tip.taxon.label for tip in current.leaf_node_iter()]

collapsed_tax = []
for node in current:
    if not node.is_leaf():
        if node.label in name_map:
            collapsed_tax.append(node.label)
            node.clear_child_nodes()
            node.taxon = current.taxon_namespace.new_taxon(label=node.label)

leaves_A = [tip.taxon.label for tip in current.leaf_node_iter()]


taxa_retain = [key for key in name_map]

current.retain_taxa_with_labels(taxa_retain)
current.write_to_path(dest="{}/pruned.tre".format(custom_synth_dir), schema = "newick")

leaves_B = [tip.taxon.label for tip in current.leaf_node_iter()]


#now dates
maps = chronogram.map_conflict_nodes_tree(current)


dates = chronogram.build_synth_node_source_ages()


root_node = OT.synth_mrca(node_ids=leaves_B).response_dict['mrca']['node_id']


max_age_est = 130
#https://academic.oup.com/sysbio/article/63/3/442/1649269

chronogram.date_tree(maps['tree'],
                     dates,
                     root_node,
                     max_age_est,
                     method='bladj',
                     output_dir="{}/dates".format(custom_synth_dir),
                     summary='{}/dated.tre'.format(custom_synth_dir),
                     phylo_only=False,
                     reps=5,
                     grid=len(leaves_B))




node_annotations = annotations.generate_custom_synth_node_annotation(current, custom_synth_dir)


annotations.write_itol_relabel(name_map, filename="{}/ottlabel.txt".format(custom_synth_dir))


annotations.write_itol_conflict(node_annotations, filename="{}/conflict.txt".format(custom_synth_dir))
annotations.write_itol_support(node_annotations, filename="{}/support.txt".format(custom_synth_dir), param="support")


jetz_annotations = annotations.generate_custom_synth_source_traversal(current, custom_synth_dir, "ot_809@tree2")
annotations.write_itol_conflict(jetz_annotations, filename="{}/jetz_conflict.txt".format(custom_synth_dir), max_conflict=1)
annotations.write_itol_support(jetz_annotations, filename="{}/jetz_support.txt".format(custom_synth_dir), param="support", max_support = 1)
