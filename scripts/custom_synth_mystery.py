import dendropy
from opentree import OT


#wget http://ot38.opentreeoflife.org/v3/tree_of_life/custom_built_tree/snacktavish_aves_81461_tmplj_gzdsd.tar.gz
# tar -xzvf snacktavish_aves_81461_tmplj_gzdsd.tar.gz

custom_synth_dir = "custom_synth_runs/snacktavish_aves_81461_tmplj_gzdsd"

custom_tree = dendropy.Tree.get_from_path("{}/labelled_supertree/labelled_supertree.tre".format(custom_synth_dir),
                                       schema = "newick")

septophaga = 'ott285198'

node = custom_tree.find_node_with_label(septophaga)
tips = [leaf.taxon.label for leaf in node.leaf_iter()]

custom_tree.retain_taxa_with_labels(tips)


lovette_str = OT.get_tree(study_id='pg_2591',
                          tree_id='tree6024',
                          tree_format='newick',
                          label_format='ot:ottid').response_dict['content']

lovette_str=lovette_str.decode()
lovette_str=lovette_str.replace(",3599881:18.535789",",3599881_B:18.535789") #UGH
lovette_str=lovette_str.replace(",3599828:25.067004",",3599828_B:25.067004") #UGH


lovette_tree = dendropy.Tree.get(data=lovette_str, schema='newick')
for taxon in lovette_tree.taxon_namespace:
    taxon.label = 'ott' + taxon.label

lovette_tree.retain_taxa_with_labels(tips)

custom_tree.print_plot()
lovette_tree.print_plot()