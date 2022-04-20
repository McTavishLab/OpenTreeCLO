import json
import time
from helpers import crosswalk_to_dict
from opentree import OT

#wget https://api.opentreeoflife.org/v3/collection/snacktavish/aves.json

collection = json.load(open("../taxonomy_info/aves.json"))


name_map = crosswalk_to_dict("../taxonomy_info/OTT_eBird_combined_taxonomy_2021.tsv")

unmapped_tip_labels = open("../taxonomy_info/unmatched_labels.tsv",'w')


for tree in collection['data']['decisions']:
    unmapped_tip_labels.flush()
    time.sleep(10)
    study_id = tree['studyID']
    tree_id = tree['treeID']
    cite = tree.get('ot:studyPublicationReference', tree['name'])
    link = "https://tree.opentreeoflife.org/curator/study/edit/{}/?tab=otu-mapping".format(study_id)
    study_data = OT.get_study(study_id).response_dict['data']['nexml']
    otus_in_tree = set()
    for treesid, treesblock in study_data['treesById'].items():
        if tree_id in treesblock['treeById'].keys():
           otus_in_tree = set([nodeinfo.get('@otu',None) for node, nodeinfo in treesblock['treeById'][tree_id]['nodeById'].items()])
    assert len(otus_in_tree) >= 3
    for otuset in study_data['otusById']:
        for otu, mapping in study_data['otusById'][otuset]['otuById'].items():
            if otu in otus_in_tree:
                if 'ott'+str(mapping.get('^ot:ottId','')) not in name_map:
                    unmapped_tip_labels.write("{}\t{}\t{}\t{}\n".format(mapping['^ot:originalLabel'],
                                                                    study_id,
                                                                    link,
                                                                    cite
                                                                    ))



unmapped_tip_labels.close()