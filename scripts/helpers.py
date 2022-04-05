def crosswalk_to_dict(taxon_crosswalk_file,
                      ott_id_key = "ott_id",
                      alt_name='SCI_NAME',
                      tax_filter = 'species'):
    infi = open(taxon_crosswalk_file)
    header = infi.readline().split('\t')
    ott_id_index = header.index(ott_id_key)
    alt_name_index = header.index(alt_name)
    name_map = {}
    no_maps = []
    double_maps = {}
    for lin in infi:
        lii = lin.split('\t')
        if lii[1] == tax_filter:
            if lii[ott_id_index]=='-':
                no_maps.append(lii)
            else:
                ott_id = 'ott'+lii[ott_id_index]
                if ott_id in name_map:
                    if ott_id in double_maps:
                        double_maps[ott_id].append(dict(zip(header, lii)))
                    else:
                        double_maps[ott_id ]=[dict(zip(header, lii))]
                name_map[ott_id] = lii[alt_name_index]
    infi.close()
    return name_map

