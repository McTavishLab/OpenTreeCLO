import dendropy
import sys
from helpers import crosswalk_to_dict

taxonomy_tre = dendropy.Tree.get_from_path("/home/ejmctavish/projects/otapi/OpenTreeCLO/custom_synth_runs/snacktavish_aves_81461_tmpk3w6j440/cleaned_ott/cleaned_not_updated_ott.tre",schema="newick", preserve_underscores=False)

name_map = crosswalk_to_dict("/home/ejmctavish/projects/otapi/OpenTreeCLO/taxonomy_info/OTT_eBird_combined_taxonomy_2021.tsv")

fam_name_map = crosswalk_to_dict("/home/ejmctavish/projects/otapi/OpenTreeCLO/taxonomy_info/OTT_eBird_combined_taxonomy_2021.tsv", alt_name="FAMILY")
## function to get back walk from ott o clo

for name in fam_name_map:
    fam_name_map[name] = fam_name_map[name].split()[0]


families = set(fam_name_map.values())

families_according_to_CLO = {fam:[] for fam in families}
for ott in fam_name_map:
    families_according_to_CLO[fam_name_map[ott]].append(name_map[ott])



leaves_A = []
for tip in taxonomy_tre.leaf_node_iter():
    if tip.taxon:
        leaves_A.append(tip.taxon.label)
    else:
        leaves_A.append(tip.label)


sys.stdout.write("Total number of tips in taxonomy_tre tree is {}\n".format(len(leaves_A)))

name_to_ott_id = {}
collapsed_tax = []
for node in taxonomy_tre:
    if node.taxon:
        node.taxon.label = node.taxon.label.replace("_"," ")
        name = " ".join(node.taxon.label.split()[:-1])
        ott_id = node.taxon.label.split()[-1]
        node.taxon.label = ott_id
        name_to_ott_id[name] = ott_id
    else:
        node.label = node.label.replace("_"," ")
        name = " ".join(node.label.split()[:-1])
        ott_id = node.label.split()[-1]
        node.label = ott_id
        name_to_ott_id[name] = ott_id
    if not node.is_leaf():
        if node.label in name_map:
            collapsed_tax.append(node.label)
            node.clear_child_nodes()
            node.taxon = taxonomy_tre.taxon_namespace.new_taxon(label=node.label)



leaves_B = []
for tip in taxonomy_tre.leaf_node_iter():
    if tip.taxon:
        leaves_B.append(tip.taxon.label)
    else:
        leaves_B.append(tip.label)

assert 'ott3598459' in leaves_B

sys.stdout.write("Total number of tips in synth tree after collapsing subspecies is {}\n".format(len(leaves_A)))


### NOW find them in the taxonomy
families_according_to_ott = {fam:[] for fam in families}

not_in_tax_tre = []
for fam in families:
    if not name_to_ott_id.get(fam):
        not_in_tax_tre.append(fam)
    else:
        node = taxonomy_tre.find_node_with_label(name_to_ott_id[fam])
        if node:
            families_according_to_ott[fam]=[name_map.get(tip.taxon.label,"None") for tip in node.leaf_iter()]
        else:
            print(name_to_ott_id[fam])

reverse_ott_fam = {}

for fam in families_according_to_ott:
    for name in families_according_to_ott[fam]:
        reverse_ott_fam[name] = fam

reverse_clo_fam = {}

for fam in families_according_to_CLO:
    for name in families_according_to_CLO[fam]:
        reverse_clo_fam[name] = fam

ofi = open("family_names.txt","w")


reverse_name_map = {v: k for k, v in name_map.items()}
ott_id_to_name = {v: k for k, v in name_to_ott_id.items()}

for fam in families:
    ofi.write("\n-----------------------------------\n")
    ofi.write(fam+'\n')
    ofi.write("{} spp in family in clo tax\n".format(len(set(families_according_to_CLO[fam]))))
    if not name_to_ott_id.get(fam):
        ofi.write("!! family missing from OTT taxonomy tree\n")
    ofi.write("\nIn family in CLO but not in fam in OTT:\n")
    names = set(families_according_to_CLO[fam]).difference(set(families_according_to_ott[fam])) 
    for name in names:
        ofi.write("{} is in {} in ott, ({}, {})\n".format(name, reverse_ott_fam.get(name), reverse_name_map[name], ott_id_to_name.get(reverse_name_map[name], 'unk')))
    ofi.write("\nIn family in OTT but not in fam in CLO:\n")
    names = set(families_according_to_ott[fam]).difference(set(families_according_to_CLO[fam])) 
    if "None" in names:
        names.remove("None")
    for name in names:
        ofi.write("{} is in {} in clo, ({}, {})\n".format(name, reverse_clo_fam.get(name), reverse_name_map[name], ott_id_to_name.get(reverse_name_map[name], 'unk')))

ofi.close()


