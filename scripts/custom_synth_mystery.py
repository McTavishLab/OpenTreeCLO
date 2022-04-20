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
lovette_str=lovette_str.replace("3599881","3599881_A",1) #UGH
lovette_str=lovette_str.replace("3599828","3599828_A",1) #UGH


lovette_tree = dendropy.Tree.get(data=lovette_str, schema='newick')
for taxon in lovette_tree.taxon_namespace:
    taxon.label = 'ott' + taxon.label

lovette_tree.retain_taxa_with_labels(tips)

custom_tree.print_plot()
lovette_tree.print_plot()

a_str = OT.get_tree(study_id='ot_290',
                          tree_id='tree1',
                          tree_format='newick',
                          label_format='ot:ottid').response_dict['content']
a_str=a_str.decode()

dup_labels = ['31195','21001','136022', '906629', '285198','758213', '5266838', '489145','939936']
for lab in dup_labels:
    a_str = a_str.replace(lab,lab+'A',1)

a_str=a_str.replace("285198","285198B",2) #UGH, '285198'
a_str=a_str.replace("285198B","285198C",1) #UGH, '285198'

a_str=a_str.replace("5266838","5266838B",2) #UGH, '285198'
a_str=a_str.replace("5266838B","5266838C",1)



a_tree = dendropy.Tree.get(data=a_str, schema='newick')


for taxon in a_tree.taxon_namespace:
    taxon.label = 'ott' + taxon.label



a_tree.retain_taxa_with_labels(tips)
a_tree.print_plot()