## Get citation weights

import json
import dendropy
from opentree import OT




custom_synth_dir = "custom_synth_runs/snacktavish_aves_81461_tmpvu7e5t2x"
custom_tree = dendropy.Tree.get_from_path("{}/labelled_supertree/labelled_supertree.tre".format(custom_synth_dir), schema="newick")



## Pruning could happen here 


node_annotations = annotations.generate_custom_synth_node_annotation(current, custom_synth_dir)
