#This script stich fragment together using graph based method!

import glob
import sys
import re
import numpy as np
import os
from collections import defaultdict

#Some function
def get_peptide(PDB,chain="B"):
    tempPep = []
    for res in PDB:
        if res["cha"] == chain:
            tempPep.append(res)
    return(tempPep)
def print_pdb(fname,PDB):
    f = open(fname, 'w')
    for res in PDB:
        #print("ATOM      1  N   ASP A 172      21.902   0.255   5.471  1.00  34.76")
        ToPrint = (int(res["num"]),res["atom"],res["res"],res["cha"],int(res["num"]),res["coord"][0],res["coord"][1],res["coord"][2])
        f.write("ATOM%7d %s %s %s%4d%12.3f%8.3f%8.3f  1.00   1.00\n" % ToPrint)
        #print(res)
        #break
    f.close()

#Function found on the internet that find longest path in graph
def DFS(G,v,seen=None,path=None):
    if seen is None: seen = []
    if path is None: path = [v]

    seen.append(v)

    paths = []
    for t in G[v]:
        if t not in seen:
            t_path = path + [t]
            paths.append(tuple(t_path))
            paths.extend(DFS(G, t, seen[:], t_path))
    return paths

def read_pdb_seq(myF,keep_hydrogen = 1,keep_het = 1,mask=""):
    lines = [line.rstrip('\n') for line in open(myF)]
    Nresnumc = "NA"
    myPDB = []
    myC = 0
    pattern = re.compile(mask)
    for l in lines:
     # print(l)
        # My hash
        raw = {}
        if (re.search(' H\s*$',l) != None and keep_hydrogen == 0):
            continue
        if (re.search('^HET',l) != None and keep_het == 0):
            continue
        if mask != "":
            if not pattern.search(l):
                continue
        m = re.search('^ATOM........(....).(...) (.)\s*(\d+\S*)\s*(-*\d+\.\d{3})\s*(-*\d+.\d{3})\s*(-*\d+.\d{3})', l)
        if (m == None):
            #print('Can t find patern in '+l)
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

#You argument is just the folder where DECORATES was run

#change folder
os.chdir(sys.argv[1])


#Pertubation for grid (+/- 1 grid point around)
##All change
AllChange = []
for x in range(-1,2):
    for y in range(-1,2):
        for z in range(-1,2):
            AllChange.append([x,y,z])
            
#Some grid param that I'm not sure exactly sure what they do
#Don't touch, unless you need to
Step = 1
Prec = 1.0

#Dict that point toward other file (where CA are found)
Sector = dict()

#Find possible overlap
#(Min overlap node need to have to be stich)

Overlap = 5
MinOver = 3 # You don't want to go below 3, cuz you have some real weird structure that get created


#Glob all file
AllFile = glob.glob("./nchains1_nres5_*.pdb") #In the folder where the jobs was run

#Load AllFile
print("Loading %d File(s)" % (len(AllFile)))
for i in range(len(AllFile)):

    #Some output to know what have been run
    if i % int(len(AllFile)/100) == 0:
        print(i,float(i)/float(len(AllFile))*100.0) 
        
    #The file
    f = AllFile[i]
    #Read PDB
    pdb = read_pdb_seq(f,mask=" CA ")
    #Find chain
    ch = re.search("chids(.)",f).group(1)
    pdb = get_peptide(pdb,chain=ch)
    
    #Store CA in grid (sector) and their adjacent sector
    for i in range(len(pdb)):
        if i+1 not in Sector:
            Sector[i+1] = dict()
        if i-len(pdb) not in Sector:
            Sector[i-len(pdb)] = dict()
        r = pdb[i]
        #Build key
        r["name"] = f
        
        #print(i+1,i-len(pdb))
        #print(r)
        
        for c in AllChange:
            xc = np.array(c)*Step+np.array(r['coord'])
            xc = (xc*Prec)
            xc = xc.astype(int).astype(str)
            xc = list(xc)
            key = "_".join(xc)
            if i < Overlap:
                if key not in Sector[i+1]:
                    Sector[i+1][key] = []
                Sector[i+1][key].append(r)
            if i-len(pdb) > -Overlap-1:
                if key not in Sector[i-len(pdb)]:
                    Sector[i-len(pdb)][key] = []
                Sector[i-len(pdb)][key].append(r)


# In[19]:

print("Finding Match")
#PotMatch
AllMatch = dict()
G = defaultdict(list)

#Node distance threshold
dthre = 0.75

#Count comp
count = 0

#Try to find wich node overlap
for NowOver in range(MinOver,Overlap+1):
    print(NowOver)
    
    #Potential match
    PotMatch = dict()
    
    #For all potential over
    for i in range(0,NowOver):
    
        #Nterm and Cterm that overlap
        Nterm = i+1
        Cterm = -NowOver+i
        print("Nterm",Nterm,len(Sector[Nterm]))
        print("Cterm",Cterm,len(Sector[Cterm]))
        
        
        #For each sector in Nterm
        for k in sorted(Sector[Nterm]):
        
            #If nterm sector not in cterm, continue
            if k not in Sector[Cterm]:continue
            
            #For all cterm sector
            for r1 in Sector[Cterm][k]:
            
                #For all nterm sector
                for r2 in Sector[Nterm][k]:
                
                    #If same residue, next
                    if r1["name"] == r2["name"]:continue
                    
                    #Some index
                    ck = r1["name"]+ " " + r2["name"]
                    
                    #Check if already done   
                    
                    #If in match, but not what expected                  
                    if ck in PotMatch:
                        if len(PotMatch[ck]) != (i)*3:
                            continue
                    else:
                        #If not already match and past frst iteration
                        if i != 0:
                            continue
                            
                    #Some dist form RMSD
                    dist = np.array(r1["coord"])-np.array(r2["coord"])
                    
                    #Some threshold
                    if np.sum(np.power(dist,2)) > dthre*dthre:
                        #dist += [9999,9999,9999]
                        if ck in PotMatch:
                            del PotMatch[ck]
                        count += 1
                        continue
                    #print(r1["resnumc"],r2["resnumc"])
                    if ck not in PotMatch:
                        PotMatch[ck] = []
                    PotMatch[ck] = np.append(PotMatch[ck],dist)
                    count += 1
        print("Len PotMatch",len(PotMatch),"Count",count)
    for k in sorted(PotMatch.keys()):
        
        if len(PotMatch[k]) < NowOver*3:
            continue
        rmsd = np.mean(np.power(PotMatch[k],2))
        if rmsd > 0.5:
            continue
        #print(k,rmsd)
        [s,t] = k.split(" ")[0:2]
        G[s].append(t)
        AllMatch[k] = NowOver
        


# In[20]:

#Find all the paths
all_paths = []
for k in list(G.keys()):
    all_paths += DFS(G, k)
max_len   = max(len(p) for p in all_paths)
max_paths = [p for p in all_paths if len(p) == max_len]


# In[25]:

#Write the binding stiched binding pose
AtomToKeep = [" N  "," CA "," C  "," O  "," CB "]

directory = "./binding_pose"
if not os.path.exists(directory):
    os.makedirs(directory)


Id = 0
for p in all_paths:
    if len(p) < 2:
        continue
        
    #Merge stuff
    #print(p)
    MergeCoord = dict()
    PathOver = [0]
    
    #AllOverlap
    MergeOver = []
    
    AllFrag = dict()
    for f in p:
        #Read PDB
        pdb = read_pdb_seq(f,mask=" CA ")
        #Find chain
        ch = re.search("chids(.)",f).group(1)
        pdb = get_peptide(pdb,chain=ch)
        AllFrag[f] = pdb
    
    
    for i in range(1,len(p)):
        k = p[i-1] + " " + p[i]
        pdb1 = AllFrag[p[i-1]]
        #print(i,p[i],AllMatch[k],len(pdb1))
        PathOver.append(PathOver[-1]+len(pdb1)-AllMatch[k])
        MergeOver.append(AllMatch[k])
    #print(PathOver)
    Score = 0.0
    for i in range(len(p)):
        Score += float(re.search("score(-*\d+\.*\d+)_",p[i]).group(1))
        pdb1 = AllFrag[p[i]]
        for j in range(len(pdb1)):
            nci = PathOver[i] + j
            if nci not in MergeCoord:
                MergeCoord[nci] = []
            MergeCoord[nci].append(p[i]+" "+pdb1[j]["resnumc"])
    print(len(p),Score,len(MergeCoord),np.mean(MergeOver),np.min(MergeOver))
    #if len(MergeCoord) < 8:
    #    continue
    MergePDB = []
    for k in MergeCoord:
        AllCord = []
        #print(k,MergeCoord[k])
        for m in range(len(AtomToKeep)):AllCord.append([])
        for m in MergeCoord[k]:
            f = m.split(" ")[0]
            resnumc = " ".join(m.split(" ")[1:])
            #print(k,f,resnumc)
            fpdb = read_pdb_seq(f)
            ch = re.search("chids(.)",f).group(1)
            fpdb = get_peptide(fpdb,chain=ch)
            for r in fpdb:
                if r["resnumc"] == resnumc:
                    if r["atom"] in AtomToKeep:
                        ind = AtomToKeep.index(r["atom"])
                        AllCord[ind].append(r["coord"])

        for m in range(len(AtomToKeep)):
            AvgCoord = np.sum(AllCord[m],axis=0)/float(len(MergeCoord[k]))
            if (len(AllCord[m])) == 0:continue
            r = dict()
            r["ID"] = k
            r["atom"] = AtomToKeep[m]
            r["cha"] = "Z"
            r["coord"] = list(AvgCoord)
            r["num"] = k+1
            r["res"] = "ALA"
            MergePDB.append(r)
    #die
    Id += 1
    
    #nchains1_nres10_score1.922_-0.200_rmsd0.287_term049253_chidsA_match40.pdb
    
    merge_name = ["nchains1","nres"+str(len(MergeCoord)),"score"+"%.3f" % (Score),"0.000","rmsd0.000","term"+str(len(p)),"chidsZ","match"+str(Id)]
    
    #fname = "./binding_pose/" + "_".join(["res"+str(len(MergeCoord)),"score"+"%.3f" %(Score),"frag"+str(len(p)),"MinOver"+str(int(np.min(MergeOver))),str(Id)]) + ".pdb"
    fname = "./binding_pose/" + "_".join(merge_name) + ".pdb"
    print_pdb(fname,MergePDB)


