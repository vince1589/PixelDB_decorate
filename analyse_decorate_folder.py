import sys
import re
import math
import numpy as np
import glob
import os
import pickle
import copy

#Get the peptide!
def get_peptide(PDB,chain="B"):
    tempPep = []
    for res in PDB:
        if res["cha"] == chain:
            tempPep.append(res)
    return(tempPep)

#Read the PDB
def read_pdb_seq(myF,keep_hydrogen = 1,keep_het = 1,mask=""):
    lines = [line.rstrip('\n') for line in open(myF)]
    Nresnumc = "NA"
    myPDB = []
    myC = 0
    for l in lines:
     # print(l)
        # My hash
        raw = {}
        if (re.search(' H\s*$',l) != None and keep_hydrogen == 0):
            continue
        if (re.search('^HET',l) != None and keep_het == 0):
            continue
        if mask != "":
            pattern = re.compile(mask)
            if not pattern.search(l):
                continue
        m = re.search('^ATOM........(....).(...) (.)\s*(\d+\S*)\s*(-*\d+\.\d{3})\s*(-*\d+.\d{3})\s*(-*\d+.\d{3})', l)
        if (m == None):
            print('Can t find patern in '+l)
            continue
        # Get resnumc
        resnumc = str(m.group(2)) + " " + str(m.group(4)) + " " + str(m.group(3))
     # print(m.group(4,5,6))
        myCord = m.group(5,6,7)
        nCord = [float(i) for i in myCord]
        raw['coord'] = nCord
        raw['res'] = m.group(2)
        raw['num'] = m.group(4)
        raw['cha'] = m.group(3)
        raw['atom'] = m.group(1)
        raw['resnumc'] = resnumc
        raw['ID'] = myC
        myC += 1
        myPDB.append(raw)
     # break
    return(myPDB)

#GET tm between 2 pdb (will slide one versus the other)    
def Get_TM(HoloCord,TermCord,cut=1.0,MinMatch=2,strict=0):
    MaxMatch = 0
    BestAli = "-"*len(HoloCord)
    BestDist = 99999
    BestDistArr = []
    for k in range(-len(HoloCord)+MinMatch,len(HoloCord)-MinMatch-1):
        TotMatch = 0
        Ali1 = ""
        TotDist = []
        for j in range(0,len(HoloCord)):
            if j + k < 0:
                Ali1 += "-"
                continue
            if (j+k) >= len(TermCord):
                Ali1 += "-"
                continue
            Dist = np.sum(np.power(TermCord[j+k]-HoloCord[j],2))
            if Dist < cut*cut:
                TotMatch += 1
                TotDist.append(Dist)
                Ali1 += "!"
            else:
                #TotMatch = -9999
                Ali1 += "."
                if strict == 1:
                    break
            #    break
        #print(cut*cut,TotMatch,np.sqrt(np.mean(TotDist)),TotDist)
        if TotMatch > MaxMatch:
            MaxMatch = TotMatch
            BestAli = Ali1
            BestDist = np.mean(TotDist)
            BestDistArr = copy.copy(TotDist)
        if (TotMatch == MaxMatch) & (MaxMatch != 0):
            if np.mean(TotDist) < BestDist:
                
                BestDist = np.mean(TotDist)
                MaxMatch = TotMatch
                BestAli = Ali1
                BestDistArr = copy.copy(TotDist)
    #if MaxMatch > 1:
    #    print(BestAli,MaxMatch,len(HoloCord),len(TermCord))
                
    return((MaxMatch,BestAli,BestDistArr))
   
   
def get_coord(pdb):
    Coord = []
    for res in pdb:
        Coord.append(np.array(res["coord"]))
    return(Coord)
    
#First argument is the holo PDB
holo = sys.argv[1]

#Peptide chain (comma separated)
holo_ch = sys.argv[2] 

#Where to write the results! (Pickle is your friend!)
picklefile = sys.argv[3]

#List of file to bench
FragTERMListFile = sys.argv[4]




#Pickle need to be pk (This is a weird thing)        
if ".pdb" in picklefile:
    print(picklefile)
    sys.exit()

#Get the peptide        
FullPDB = read_pdb_seq(holo,mask="ATOM.* CA ")
HoloPep = []
for ch in holo_ch.split(","):
    HoloPep.append(get_coord(get_peptide(FullPDB,chain=ch)))


#Other Receptor in the cluster (pdbfile C,H,A,I,N)
AddHolo = []
if len(sys.argv) > 5:
    for TP in [line.rstrip('\n') for line in open(sys.argv[5])]:
        [tpdbf,tch] = TP.split(" ")
        TPdb = read_pdb_seq(tpdbf,mask="ATOM.* CA ")
        for ch in tch.split(","):
            AddHolo.append(get_coord(get_peptide(TPdb,chain=ch)))

#Dict that will have the results data
SimilToHolo = dict()

#Load some data if already ran
if os.path.isfile(picklefile):
    SimilToHolo = pickle.load( open( picklefile, "rb" ))


#Load list of fragment (peptide)
FragTERMList = [line.rstrip('\n') for line in open(FragTERMListFile)]

BestRMSDEVer = 9999999
for p in FragTERMList:
    #If already done next
    #if p in SimilToHolo:
    #    continue
    #Need to be a PDB    
    if ".pdb" not in p:
        continue
    #Need to respect some nomenclature that is part of decorate (hopefullly still the same!)
    if "chids" not in p:
        print(p)
        continue
    #Find chain
    chs = list(re.search("chids(.*)_match",p).group(1))
    #print(p,chs)

    #Get Info
    kp = p #Entry name, could be change for something smaller
    SimilToHolo[kp] = dict()
    for sp in re.split("/",p)[-1].split(".pdb")[0].split("_"):
        m = re.search("^(chids)(.*)$",sp)
        if m != None:
            SimilToHolo[kp][m.group(1)] = str(m.group(2))
            continue
        m = re.search("^([A-Za-z]+)(\d+$)",sp)
        if m != None:
            SimilToHolo[kp][m.group(1)] = int(m.group(2))
            continue
        m = re.search("^([A-Za-z]+)(-*\d+\.\d+$)",sp)
        if m != None:
            SimilToHolo[kp][m.group(1)] = float(m.group(2))                                  
            continue
        m = re.search("^(-*\d+\.\d+$)",sp)
        if m != None:
            SimilToHolo[kp]["Score2"] = float(m.group(1))                                  
            continue
        
    
    #Get TERM pdb seq for comparison
    TermPep = read_pdb_seq(p,mask=" CA .* ")
    
    #Is it overlaping any of the receptor surface peptide
    
    
    #Init the peptide sequence as non-overlaped
    InitAli = ""
    for HoloCord in HoloPep:
        InitAli += "-"*len(HoloCord)
    InitAli = list(InitAli)
    
    #Init Tot overlap position on peptide(s)
    TotAli = 0
    
    #Kept best distance from peptide
    BDarr = []
    
    #For Every Segments chains
    for ch in list(chs):
    
        #Get Coord
        TermCord = get_coord(get_peptide(TermPep,chain=ch))
        
        #Init alignment
        alns = ""
        
        #For every receptor peptide(s)
        BestRMSDarr = []
        BestRMSD = 999999
        for HoloCord in HoloPep:
             #Get Alignment
            (ali,BestAli,peprmsd)= Get_TM(HoloCord,TermCord,cut=1.5,MinMatch=1)
            #Concat Ali
            alns += BestAli
            
            #Store number of node aligned on surface recep
            TotAli += ali
            
            #Get RMSD
            (ali,BestAli,peprmsd)= Get_TM(HoloCord,TermCord,cut=9999)
            if len(peprmsd) == 0:
                print(p)
                die
            if np.mean(peprmsd) < BestRMSD:
                BestRMSD = np.mean(peprmsd)
                BestRMSDarr = copy.copy(peprmsd)
        BDarr += BestRMSDarr
        #Update the alnment string for every segment    
        for i in range(len(alns)):
            if alns[i] == "!":
                InitAli[i] = "!"
        
    InitAli = "".join(InitAli)
    SimilToHolo[kp]["TotAli"] = TotAli
    
    SimilToHolo[kp]["Alignment"] = InitAli
    if "binding_pose" in p:
        SimilToHolo[kp]["Merge"] = 1
    else:
        SimilToHolo[kp]["Merge"] = 0
    if len(BDarr) == 0:
        continue
    SimilToHolo[kp]["PepRMSD"] = np.sqrt(np.mean(BDarr))
    if np.sqrt(np.mean(BDarr)) < BestRMSDEVer:
        BestRMSDEVer =    np.sqrt(np.mean(BDarr))                     
        print(p,SimilToHolo[kp])
    #print(len(SimilToHolo[kp]),SimilToHolo[kp])
    
    
pickle.dump(SimilToHolo, open( picklefile, "wb" ))




