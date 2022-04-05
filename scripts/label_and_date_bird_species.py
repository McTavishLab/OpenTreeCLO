import dendropy
import copy
from chronosynth import chronogram
from opentree import OT, annotations, taxonomy_helpers


infi = open("../taxonomy_info/OTT_eBird_combined_taxonomy_2021.tsv")
header = infi.readline().split('\t')


species = {}
no_maps = []
double_maps = {}

for lin in infi:
    lii = lin.split('\t')
    if lii[1]=='species':
        if lii[9]=='-':
            no_maps.append(lii)
        else:
            if 'ott'+lii[9] in species:
                if 'ott'+lii[9] in double_maps:
                    double_maps['ott'+lii[9]].append(dict(zip(header, lii)))
                else:
                    double_maps['ott'+lii[9]]=[dict(zip(header, lii))]
                    double_maps['ott'+lii[9]].append(species['ott'+lii[9]])
            species['ott'+lii[9]]=dict(zip(header, lii))


infi.close()


doubles = open("double_match_spp.tsv",'w')
doubles.write('\t'.join(header))

for ott_id in double_maps:
    for match in double_maps[ott_id]:
        doubles.write('\t'.join([match[key] for key in header]))

doubles.close()

missing = open("no_match_spp.tsv",'w')
missing.write('\t'.join(header))
for item in no_maps:
    missing.write('\t'.join(item))

missing.close()




ott_ids = list(species.keys())

resp = OT.synth_induced_tree(node_ids = ott_ids, label_format = 'id', ignore_unknown_ids=True)


tree = dendropy.Tree.get_from_string(resp.response_dict['newick'], schema= "newick")
tree.write(path='id.tre', schema='newick')

node_annotations = annotations.generate_synth_node_annotation(tree)

annotations.write_itol_conflict(node_annotations)
annotations.write_itol_support(node_annotations)


#tree_dict=taxonomy_helpers.labelled_induced_synth(ott_ids=ott_ids, label_format='name')

#for label in tree_dict['label_map']:
#    ott_id = label.split()[-1]
#    translation_dict[ott_id] = tree_dict['label_map'][label]



translation_dict = {}
for node in node_annotations:
    if node in species:
        translation_dict[node]=species[node]['SCI_NAME']

annotations.write_itol_relabel(translation_dict, filename="ottlabel.txt")


chronogram.date_synth_subtree(node_ids=ott_ids,
                              max_age=165,
                              output_dir='.',
                              method="bladj")


## Summarize facts about tree
## How many taxa, how many branches,
## Crosswalk table for tip labels to clements, ncbi, gbif


tips = [tip.taxon.label for tip in tree.leaf_node_iter()]



missing_tips = set(ott_ids).difference(set(tips))
surprise_tips = set(tips).difference(set(ott_ids))