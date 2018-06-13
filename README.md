# Residue_interaction_network

Used to find interaction networks of a set of residues in a protein. The main function is
'find_interactions(pdb,res_to_run,folder_name, save_folder)'

## Input
pdb: 4 letter PDB code.\
res_to_run: array of residue ids to cluster\
folder_name: folder where PDB.mat files are stores\
save_folder: folder to save results to\

## Output
\* _all_core_binned_2.mat
- Column 1: Residue Ids
- Column 2: Which group the amino acid is in
- All pairs that are not ILMTSWYV are removed and put in individual groups. To change this, comment out as indicated

\*_all_paired_data_2.mat: an nxn array with 1s for each interaction pair. All amino acids in res_to_run are included\

## Notes
- Interactions are defined as 2 residues being able to overlap without first overlapping with the backbone
- Set up PDB file using download_preprocess_pdb.py

## Citations
J. C. Gaines, A. Virrueta, D. A. Buch, S. J. Fleishman, C. S. O'Hern, and L. Regan, Collective repacking reveals that the structures of protein cores are uniquely specified by steric repulsive interactions, Protein Eng. Des. Sel. 30 (2017) 387.
