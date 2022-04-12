import sys
import os
import json
import dendropy
import copy
from chronosynth import chronogram
from opentree import OT, annotations, taxonomy_helpers
from helpers import crosswalk_to_dict


OT.api_endpoint = "development"
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
leaves_B = [tip.taxon.label for tip in current.leaf_node_iter()]

'''
current.write_to_path(dest="{}/pruned.tre".format(custom_synth_dir), schema = "newick")



#now dates

custom_str = chronogram.conflict_tree_str(current)
custom_dates = chronogram.build_synth_node_source_ages(compare_to = custom_str,
                                                       fresh = True)


root_node = OT.synth_mrca(node_ids=leaves_B).response_dict['mrca']['node_id']


max_age_est = 130
#https://academic.oup.com/sysbio/article/63/3/442/1649269

chronogram.date_tree(current,
                     custom_dates,
                     root_node,
                     max_age_est,
                     method='bladj',
                     output_dir="{}/dates".format(custom_synth_dir),
                     summary='{}/dated.tre'.format(custom_synth_dir),
                     phylo_only=False,
                     reps=5,
                     grid=len(leaves_B))



'''

node_annotations = annotations.generate_custom_synth_node_annotation(current, custom_synth_dir)


annotations.write_itol_relabel(name_map, filename="{}/ottlabel.txt".format(custom_synth_dir))


annotations.write_itol_conflict(node_annotations, filename="{}/conflict.txt".format(custom_synth_dir))
annotations.write_itol_support(node_annotations, filename="{}/support.txt".format(custom_synth_dir), param="support")


jetz_annotations = annotations.generate_custom_synth_source_traversal(current, custom_synth_dir, "ot_809@tree2")
annotations.write_itol_conflict(jetz_annotations, filename="{}/jetz_conflict.txt".format(custom_synth_dir), max_conflict=1)
annotations.write_itol_support(jetz_annotations, filename="{}/jetz_support.txt".format(custom_synth_dir), param="support", max_support = 1)

clem_annotations = annotations.generate_custom_synth_source_traversal(current, custom_synth_dir, "ot_2019@tree7")
annotations.write_itol_conflict(clem_annotations, filename="{}/clem_conflict.txt".format(custom_synth_dir), max_conflict=1)
annotations.write_itol_support(clem_annotations, filename="{}/clem_support.txt".format(custom_synth_dir), param="support", max_support = 1)



## Pruning could happen here 


annot = json.load(open("{}/annotated_supertree/annotations.json".format(custom_synth_dir)))

study_weights = {}

uncontested_taxa = []
node_support_annotation = {}
all_nodes = []
why_no_annot = []
for node in current:
    label = None
    if not node.is_leaf():
        label = node.label
        all_nodes.append(label)
        assert label
        if label == 'ott81461':
            pass
        else:
            if label not in annot['nodes']:
                why_no_annot.append(label)
            elif 'supported_by' in annot['nodes'][label].keys():
                strict_support = annot['nodes'][label]['supported_by']
                node_support_annotation[label] = strict_support
            else:
                assert label.startswith('ott')
                uncontested_taxa.append(label)




node_phylo_count = {}
study_node_count = {}
for node in node_support_annotation:
    for source in node_support_annotation[node]:
        study_id = source.split('@')[0]
        if study_id in study_node_count:
            study_node_count[study_id]+=1
        else:
            study_node_count[study_id] = 1


study_cite_file = open("citation_node_counts.tsv", "w")
for study_id in study_node_count:
    cites = OT.get_citations([study_id]).replace('\n','\t')
    study_cite_file.write("{}\t{}\t{}\n".format(study_id, study_node_count[study_id], cites))

study_cite_file.close()


clem_conf = 0
clem_supp = 0
jetz_conf = 0
jetz_supp = 0
all_conf = 0
phylo_supp = 0
only_clem_supp = 0

for node in node_annotations:
    if jetz_annotations[node].get('support') >= 1:
        jetz_supp += 1
    if jetz_annotations[node].get('conflict') >= 1:
        jetz_conf += 1
    if clem_annotations[node].get('support') >= 1:
        clem_supp += 1
    if clem_annotations[node].get('conflict') >= 1:
        clem_conf += 1


for node in node_support_annotation:
    supp = node_support_annotation.get(node, {})
    if len(set(node_support_annotation[node].keys()).difference(set(['ot_2019@tree7']))) >=1:
        phylo_supp += 1
    if set(node_support_annotation[node].keys()) == (set(['ot_2019@tree7'])):
        only_clem_supp += 1

#    if node_annotations[node].get('conflict') >= 1:
#        all_conf += 1


summary_statement = """This tree contains {l} leaves and {i} internal nodes.\n
                        Of those nodes, {asl} are strictly supported \n
                        by at least 1 input phylogeny. The rest ({tsl}) \n
                        are placed by taxonomy.""".format(l=len(leaves_B),
                                                                i=len(all_nodes),
                                                                asl=phylo_supp,
                                                                tsl=len(all_nodes)-phylo_supp)


print(summary_statement)