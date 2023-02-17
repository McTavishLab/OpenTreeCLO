import dendropy
from opentree import OT


#wget http://ot38.opentreeoflife.org/v3/tree_of_life/custom_built_tree/snacktavish_aves_81461_tmplj_gzdsd.tar.gz
# tar -xzvf snacktavish_aves_81461_tmplj_gzdsd.tar.gz



oliveros_str = OT.get_tree(study_id='ot_2079',
                          tree_id='tree1',
                          tree_format='newick',
                          label_format='ot:ottid').response_dict['content']

oliveros_str=oliveros_str.decode()
dup_labels=['3597102', '286089', '1066395', '749724','44866']
for lab in dup_labels:
    oliveros_str = oliveros_str.replace(lab,lab+'A',1)

oliveros_str=oliveros_str.replace("749724","749724B",5) #UGH, '285198'
oliveros_str=oliveros_str.replace("749724B","749724C",1) #UGH, '285198'
oliveros_str=oliveros_str.replace("749724B","749724D",1) #UGH, '285198'
oliveros_str=oliveros_str.replace("749724B","749724E",1) #UGH, '285198'
oliveros_str=oliveros_str.replace("749724B","749724F",1) #UGH, '285198'



oliveros_tree = dendropy.Tree.get(data=oliveros_str, schema='newick')
for taxon in oliveros_tree.taxon_namespace:
    taxon.label = 'ott' + taxon.label


a_str = OT.get_tree(study_id='ot_2013',
                          tree_id='tree8',
                          tree_format='newick',
                          label_format='ot:ottid').response_dict['content']
a_str=a_str.decode()

#dup_labels = ['']
#for lab in dup_labels:
#    a_str = a_str.replace(lab,lab+'A',1)




a_tree = dendropy.Tree.get(data=a_str, schema='newick', taxon_namespace=oliveros_tree.taxon_namespace)

for taxon in a_tree.taxon_namespace:
    taxon.label = 'ott' + taxon.label


tips = set([leaf.taxon.label for leaf in oliveros_tree.leaf_node_iter()]).intersection(set([leaf.taxon.label for leaf in a_tree.leaf_node_iter()]))

oliveros_tree.retain_taxa_with_labels(tip)
a_tree.retain_taxa_with_labels(tips)


