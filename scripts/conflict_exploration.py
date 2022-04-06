import dendropy
from helpers import crosswalk_to_dict

from opentree import OT, annotations, object_conversion
from chronosynth import chronogram

taxonomy_crosswalk = "../taxonomy_info/OTT_eBird_combined_taxonomy_2021.tsv"

old_custom_synth_dir = "../custom_synth_runs/snacktavish_aves_81461_tmpvu7e5t2x"
new_custom_synth_dir = "../custom_synth_runs/snacktavish_aves_81461_tmplj_gzdsd/"
old_tree = dendropy.Tree.get_from_path("{}/labelled_supertree/labelled_supertree.tre".format(old_custom_synth_dir),
                                       schema = "newick")



new_tree = dendropy.Tree.get_from_path("{}/labelled_supertree/labelled_supertree.tre".format(new_custom_synth_dir),
                                       schema = "newick")

old_tree_str = chronogram.conflict_tree_str(old_tree)

new_tree_str = chronogram.conflict_tree_str(new_tree)

resp = OT.conflict_str(tree_str = old_tree_str,compare_to=new_tree_str)


old_node_annotations = annotations.generate_custom_synth_node_annotation(old_tree, old_custom_synth_dir)

new_node_annotations = annotations.generate_custom_synth_node_annotation(new_tree, new_custom_synth_dir)

jetz_tree = chronogram.as_dendropy("ot_809@tree2")

name_map = crosswalk_to_dict(taxonomy_crosswalk)

