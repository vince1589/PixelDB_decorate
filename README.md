# PixelDB_decorate

run_decorate.pl
To run decorate (some path might have change, double check) -> Check code for arguments!
Ex: perl run_decorate.pl -pdb /full/path/YourPDB.pdb -min 0 -max 200 -pepch B,C,D -wrkdir /home/ironfs/scratch/grigoryanlab/vfrap/decorate/pdb_name/

Random stuff:
It will create folders in the wrkdir: runs (TERMs completion PDB), done (some log of what is done), results (some results)
Many instance of this script can run on parallel, they will check if on TERMs is running or it is done
min and max are based on decorate rank order

analyse_decorate_folder.py
#Script that analyse overlap of TERMs completion to peptide
Ex: pyhton analyse_decorate_folder.py /full/path/YourPDB.pdb B,C,D /absolute/path/pickle_results.pk /wrkdir/pdb/runs/*.pdb

Random stuff:
You might want to go in the source code to see what information is stored


merge_frag_fast.py
Script that merge terms together and store them in a new folder called binding_pose
It is somewhat fast (still take time to load all the PDB)
Ex: python merge_frag_fast.py /wrkdir/

Random stuff: 
To be stitched, fragment need to be within 0.75 Angstrom and RMSD ¸0.5 (these thing can be change in source code)
File name is: PepLen_Score_FragMerge_MinOver_ID


concat_everything.pl
Script that update the run TERMs in the decorate folder (you should run this only once when no job are running)
perl concat_everything.pl /path/where/all/decorate/run/are/*



