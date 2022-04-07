## Get citation weights

import json
import dendropy
from opentree import OT



custom_synth_dir = "custom_synth_runs/snacktavish_aves_81461_tmpvu7e5t2x"
custom_tree = dendropy.Tree.get_from_path("{}/labelled_supertree/labelled_supertree.tre".format(custom_synth_dir), schema="newick")



## Pruning could happen here 


annot = json.load(open("{}/annotated_supertree/annotations.json".format(custom_synth_dir)))

study_weights = {}

nolabel = []
node_support_annotation = {}
for node in custom_tree:
    label = None
    if node.label:
        label = node.label
    elif node.taxon:
        if node.taxon.label:
            label = node.taxon.label
    try:
        strict_support = annot['nodes'][label].get('supported_by', {})
        node_support_annotation[label] = strict_support
    except:
        nolabel.append(label)

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
