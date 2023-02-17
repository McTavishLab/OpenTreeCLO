import dendropy
import sys
import json
from opentree import OT
from helpers import crosswalk_to_dict

def make_node_url(source, node):
    study, tree = source.split('@')
    return "https://tree.opentreeoflife.org/curator/study/view/{}?tab=home&tree={}&node={}".format(study, tree, node)


name_map = crosswalk_to_dict("/home/ejmctavish/projects/otapi/OpenTreeCLO/taxonomy_info/OTT_eBird_combined_taxonomy_2021.tsv")

fam_name_map = crosswalk_to_dict("/home/ejmctavish/projects/otapi/OpenTreeCLO/taxonomy_info/OTT_eBird_combined_taxonomy_2021.tsv", alt_name="FAMILY")
## function to get back walk from ott o clo

for name in fam_name_map:
    fam_name_map[name] = fam_name_map[name].split()[0]

families = set(fam_name_map.values())

families_according_to_CLO = {fam:[] for fam in families}
for ott in fam_name_map:
    families_according_to_CLO[fam_name_map[ott]].append(name_map[ott])

annotations = json.load(open("custom_synth_runs/Family_conflict_annotations.json"))


fam_tree = dendropy.Tree.get(path='custom_synth_runs/family_comparisons/cleaned_phylo/tree_ot_2019@tree11.tre',
                             schema="newick",
                             preserve_underscores=True)


families_info_sheet = open("family_conflicts.txt", "w")
for node in annotations['nodes']:
    if "ot_2019@tree11" in annotations['nodes'][node].get('supported_by',[]):
        families_info_sheet.write("\n========================================================================")
        fam_node = annotations['nodes'][node]['supported_by']['ot_2019@tree11']
        tips = [leaf.taxon.label.split('_')[-1] for leaf in fam_tree.find_node_with_label("_{}_".format(fam_node)).leaf_iter()]
        fam_name = fam_name_map[tips[0]]
        assert len(set([fam_name_map.get(tip, fam_name) for tip in tips])) == 1 
        families_info_sheet.write("\n"+fam_name+"\n")
        conflicts = annotations['nodes'][node].get('conflicts_with',{})
        families_info_sheet.write("{} studies conflict with this family\n".format(len(conflicts)))
        for source in conflicts:
            conf_tree = dendropy.Tree.get(path='custom_synth_runs/family_comparisons/cleaned_phylo/tree_{}.tre'.format(source),
                                          schema="newick",
                                          preserve_underscores=True,
                                          taxon_namespace = fam_tree.taxon_namespace)
            families_info_sheet.write("\n\n"+OT.get_citations([source]))
            for conf_node in conflicts[source]:
                families_info_sheet.write(make_node_url(source, conf_node)+'\n')
                conf_tips_all = [leaf.taxon.label for leaf in conf_tree.leaf_node_iter()]
                conf_tips_ingroup = [leaf.taxon.label for leaf in conf_tree.find_node_with_label("_{}_".format(conf_node)).leaf_iter()]
                conf_tips_outgroup = set(conf_tips_all) - set(conf_tips_ingroup)
                fam_tips_in_this_tree = [label for label in conf_tips_all if label.split('_')[-1] in tips]
            if len(fam_tips_in_this_tree) >= 2:
                families_info_sheet.write("\nTips from {} family in this tree are:\n".format(fam_name))
                families_info_sheet.write(",".join(fam_tips_in_this_tree))
                mrca = conf_tree.mrca(taxon_labels=fam_tips_in_this_tree)
                tips_from_mrca = [leaf.taxon.label for leaf in mrca.leaf_iter()]
                intruders = set(tips_from_mrca).difference(set(fam_tips_in_this_tree))
                conf_tree.retain_taxa_with_labels(tips_from_mrca)
                assert len(tips_from_mrca) == len(fam_tips_in_this_tree) + len(intruders)
                for tip in conf_tree.leaf_node_iter():
                    if tip.taxon.label in fam_tips_in_this_tree:
                        tip.taxon.label = fam_name + tip.taxon.label
                    elif tip.taxon.label:
                        tip.taxon.label = "XXXXXXX" +  fam_name_map.get(tip.taxon.label.split('_')[-1],'fam unk') + tip.taxon.label
                    else:
                        pass
                conf_tree.write(path = "family_conflict_trees/{}_{}.tre".format(fam_name, source), schema="newick")
                #families_info_sheet.write("\n\nIntruding tips are:\n")
                #for tax in intruders:
                #    families_info_sheet.write(tax+' ('+fam_name_map.get(tax.split('_')[-1],'fam unk')+"),")
            else:
                families_info_sheet.write("\nERROR: One or zero tips from family {} in tree:\n".format(fam_name))
                families_info_sheet.write(",".join(fam_tips_in_this_tree))


families_info_sheet.close()




"""
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


"""