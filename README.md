# OpenTree-CLO


To run notebooks:

    virtualenv -p python3 venv-CLOpenTree
    source venv-CLOpenTree/bin/activate 
    pip install -r requirements.txt 
    pip install ipykernel 
    python -m ipykernel install --user --name=opentree 
    pip install jupyter 
    jupyter notebook


Directory bird_studies_translation contains all the bird trees 
in Phylesystem as of May 4, 2021.


In each subfolder (named for the study_id)
Are at least 3 files - the full citation, a nexus tree, and a .csv mapping across the original label from upload, the label we currently have in opentree, and our best guess of the Clements (2019) label (or "unknown" if we don't even have a guess, and "not_a_bird" for outgroups).
There will be a tree file and a .csv for each tree associated with the study.

If you open the nexus tree in figtree you can choose the tip labels to display, from the original labels that were on the upload, or the OpenTree labels, or the Clements best guess.

Jupyter script which generated these files found in notebooks/OTT_Clements_translator.ipynb
Translation table used is taxonomy_info/combined_data.tsv


