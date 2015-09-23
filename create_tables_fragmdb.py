#Code creates tables for database

#Python 2.7

import os
import sys
import time
import collections
import fileinput
import base
import nrpdb
import Graph

def organizeFiles(models_dir, mults_dir) :
    #Decompose the BasePairs file to the set of files containing pairs by model
    #e.g. 1.txt contains pairs of 1 model
    #in  : name of future directory for files
    
    #read  BasePairs file
    base = open("BasePairs.txt", 'r')
    #create a directory for models if it doesn't exist
    if not os.path.exists(models_dir):
        os.makedirs(models_dir)
    if not os.path.exists(mults_dir):
        os.makedirs(mults_dir)
    #work in models directory
    os.chdir(models_dir)
    #list of files in the directory
    filelist = [ f for f in os.listdir('.')]
    #delete them
    for f in filelist:
        os.remove(f)
    #start to fill the directory with models files
    cur_model = '1'
    #create a file for the first model
    out = open( str(cur_model) + ".txt", 'w')
    #walk througth the pairs
    for line in base:
        model = line.split('\t')[1]
    #write pairs of the same model
        if model!=cur_model:
            out.close()
            cur_model = model
            out = open( cur_model + ".txt", 'w')
        out.write(line)
    #return to the home directory
    os.chdir("..")
    #close all the files
    out.close()
    base.close()
    return

def nuclStems( ModelPairs,M ):
    #Create dictionary {nucleotide: stems} for every model
    #in  : pairs data of the model
    #out : dictionary {nucleotide: stems}
    dNuclStems = {}
    for line in ModelPairs :
        line_list = line.split('\t')
        stems = []
        fstem = line_list[14]
        rstem = line_list[15]
        nucl1 = line_list[5]
        nucl2 = line_list[7]
        #add different types of stems into stem list
        if fstem == '\\N' and rstem == '\\N' :
            continue 
        if fstem != '\\N' :
            stems.append(int(fstem))
        if rstem != '\\N' :
            stems.append(M+int(rstem))
        #list(set(...)) to delete similar elements
        #build dictionary {nucleotide: stems}
        if nucl1 in dNuclStems: 
            dNuclStems[nucl1] += stems
        else:
            dNuclStems[nucl1] = stems
        if(nucl2 in dNuclStems):
            dNuclStems[nucl2]=list(set(dNuclStems[nucl2]+ stems))
        else:
            dNuclStems[nucl2] = stems
    return dNuclStems


def stemLinks( ModelPairs,dStemLinks, dNuclStems):
    #Create dictionary {stem: links} for every model
    #in  : pairs data of the model
    #out : dictionary {stem: links} , stemGraph - adjacency list
    stemGraph = {}
    for line in ModelPairs :
        line_list = line.split('\t')
        stems = []
        fstem = line_list[14]
        rstem = line_list[15]
        link = line_list[17]
        #if pair doesnt  belong to any stem
        if fstem == '\\N' and rstem == '\\N':
            nucl1 = line_list[5]
            nucl2 = line_list[7]
            #but if these nucleotides are paired with other nucleotides from stems
            if(nucl1 in dNuclStems and nucl2 in dNuclStems) :
                #merge their stems' connections
                 dNuclStems[nucl1] =  list(set(dNuclStems[nucl1]+ dNuclStems[nucl2]))
                 dNuclStems[nucl2] = dNuclStems[nucl1][:]
                 #build dictionary {stem: links}
                 for st in dNuclStems[nucl1]:
                     if st in dStemLinks:
                         dStemLinks[st] += [link]
                     else:
                         dStemLinks[st] = [link]
                     #delete the same elements
                     #dStemLinks[st] = list(set(dStemLinks[st]))
    #build  adjacency list
    for nucl in dNuclStems:
        #get all combinationsof nStem-1 length
        adj = []
        for st in dNuclStems[nucl]:
            vertex_stem = st
            adj = list(set(dNuclStems[nucl]) - set([vertex_stem]))
            if len(adj) > 0:
                if(vertex_stem in stemGraph):
                    stemGraph[ vertex_stem]=list(set(adj+stemGraph[ vertex_stem]))
                else:
                    stemGraph[ vertex_stem] = adj
    return stemGraph

def dfs (v,G,marked):
    #Depth-first search
    #in: vertex to start traverse
    if v in marked:  # if v has been already visited
        return
    marked.add(v) # else mark v as visited now
    for i in G[v]: # check all adjacent vertices
        if not i in marked: # if i hasn't been already visited
            dfs(i,G,marked) # visit it by run dfs

def componentsFind(G, marked):
    #Find all connected components of the graph stemGraph
    #out : lComponents - list of connected components ( list of lists )
    lComponents = []
    component = set() #current component
    for i in G : # for every vertex in graph
        if not i in marked:# check whether it is marked
            marked_prev = set(marked)# save the list of marked vertices before dfs
            dfs(i,G,marked)
            component = marked-marked_prev# compare with the list after dfs
            #and get component
            if(len(component)>1):#if it contains more than one vertex
                lComponents.append(list(component))#add it to the common list
    return lComponents


def writeMults(f,M,dStemLinks, lComponents, ModelPairs, MultsData,FragData, WingsData,EdgeData,ModelNucls, mults_dir, stop_line,statistics):
    #Write obtained multiplets in files, form fragments and edges
    #in  : model file f , lComponents - list of components, basepairs of the model, multiplets directory, stop line, statistics,
    #      multiplets table MultsData,fragments table FragData, wings table WingsData, edges table EdgeData, nucls of the model table
    #out : num of fragments, num of nucleotides involved in more that 1 fragment(errors in BasePairs.txt file)
    #      num of multiplets with nucleotides with a letter in the number
    path = mults_dir + '/' + f
    op = open(path, 'w+')
    nFragTotal = 0
    counter = 0
    counterCharNum = 0
    for comp in lComponents: #for every multiplet in the model
        multid = str(len(MultsData)+1)
        modelid = f.strip('.txt')

        fstem = []
        rstem = []
        links = []
        for st in comp:
            if  st<M:
                fstem.append(str(st))
            else:
                rstem.append(str(st-M))
            if st in dStemLinks:
                links+=dStemLinks[st]
        links = list(set(links))
        op.write(' '.join([str(st) for st in fstem])+'\n')
        op.write(' '.join([str(st) for st in rstem])+'\n')
        op.write(' '.join([str(l) for l in links])+'\n')

        MultPairs = [] #BasePairs of the multiplet
        indexes = []
        for line in ModelPairs:
            line_list = line.split('\t')
            ln_fstem = line_list[14]
            ln_rstem = line_list[15]
            ln_link = line_list[17]
            #construct a list of BaseRairs for the multiplet
            if ln_fstem in fstem or ln_rstem in rstem or ln_link in links:
                MultPairs.append(line)
                indexes.append(ModelPairs.index(line))
        #get fragments of the multiplet
        dStemNucls, Pairs, flagIFcharNum = markWings(MultPairs,M) #get dStemNucls{stem: [left][rigth wing]} vocabulary
        #get fragments
        lFragments = mergeWings(getWings(dStemNucls,len(WingsData)),Pairs,len(FragData), WingsData,ModelNucls,modelid,multid,M)

        #statistics counting
        if(flagIFcharNum):
            counterCharNum +=1
        nFrag = len(lFragments) #num of fragments
        if nFrag not in statistics:
            statistics[nFrag] = 1
        else:
            statistics[nFrag] += 1
        nFragTotal+= nFrag

        #write fragments and form fragments table
        adjList = {} # vocabulary {fragid: [adjacency list]}
        for fragment in lFragments: #for every frgament
                fragid = fragment['ID']
                adjList[fragid] = []
                op.write(' '.join([nucl for nucl in fragment['SEQ']])+'\n')
                fragstr = base.aFragment(fragid, f.strip('.txt'), multid, fragment['SEQ'], fragment['CON'], fragment['NWINGS'], fragment['INDNUCLS'], fragment['BPSNUM'],ModelNucls[fragment['SEQ'][0]],ModelNucls[fragment['SEQ'][-1]])
                #add to the frag table
                FragData.append(fragstr)

        #get edges
        invariant,nEdges,counterMult = getEdges(EdgeData, modelid, multid,MultPairs, lFragments,adjList)

        counter+=counterMult

        #write BasePairs
        for line in MultPairs:
            op.write(line)
        op.write(stop_line)

        #form Multiplets table
        for i in xrange(0, len(indexes)):
            ModelPairs[indexes[i]] = MultPairs[i]

        multstr= base.aMultiplet(multid,modelid,nFrag,nEdges,links,fstem,rstem, MultPairs, invariant)

        MultsData.append(multstr)


    op.close()
    return nFragTotal,counter, counterCharNum

def getEdges(EdgeData, modelid, multid,MultPairs, lFragments,adjList):
    #get edges of the multiplet (form edges table) and invariant of the multiplet
    # in : Edges table, id of the model, id of the multiplet, basepairs of the multiplet, list of fragments,
    #      empty vocabulary {fragid: [adjacency list]} with only keys - ids of all the fragments for the multiplet
    # out: invariant of the multiplet grapth ( nodes are fragments, edges are bonds between (nulceotides related to) them),
    #      num of edges of the multiplet, num of fragments, num of nucleotides involved in more that 1 fragment(errors in BasePairs.txt file)
    counter = 0
    edgeid = len(EdgeData)+1
    dEdges = {} #vocabulary {(frag1id, frag2id): [bond type]}
    for i in xrange(0,len(MultPairs)):
        line_list = MultPairs[i].split('\t')
        nucl1 = line_list[5]
        nucl2 = line_list[7]
        bond = line_list[8]
        if len([f['ID'] for f in lFragments if nucl1 in f['SEQ']])>1:
            counter +=1
        if len([f['ID'] for f in lFragments if nucl2 in f['SEQ']])>1:
            counter +=1
        fragid1 = [f['ID'] for f in lFragments if nucl1 in f['SEQ']][0]
        fragid2 = [f['ID'] for f in lFragments if nucl2 in f['SEQ']][0]
        #suppose the fragment with a bigger id to be on the first place in the edge
        edge  =(fragid1,fragid2) if fragid1>fragid2 else  (fragid2,fragid1)
        if edge not in dEdges:
            dEdges[edge] = [bond]
        else:
            dEdges[edge].append(bond)
        #table of multiplet pairs replenish
        MultPairs[i] +='\t'+'\t'.join([str(el) for el in [multid, edgeid,fragid1, fragid2]])

    #form edges table
    for e in dEdges:
        edge_str = base.aEdge(edgeid, modelid, multid, e[0], e[1], dEdges[e])
        EdgeData.append(edge_str)
        edgeid +=1
    #create adjacency list
    adjList = getAdjList(adjList, dEdges.keys())
    return Graph.invariant(adjList),len(dEdges),counter

def getAdjList(dAdj, edges):
    #create adjacency list, filling in empty vocabulary {fragid: [edges]} with only fragments' ids (nodes labels)
    for e in edges:
        if e[1] not in dAdj[e[0]]:
            dAdj[e[0]].append(e[1])
        if e[0] not in dAdj[e[1]]:
            dAdj[e[1]].append(e[0])
    return dAdj




def markWings(MultPairs,M):
    #build dictionary dStemNucls containing wings : {stem : [[left wing][rigth wing]]} and pairs of nucleotides of ledt-rigth wings
    #in  : Base pairs of a multiplet
    #out : dictionarydStemNucls, list of pairs (nucl1,nucl2), related to left and rigth wings accordingly
    dStemNucls={}
    Pairs = []
    flagIFcharNum = False
    for line in MultPairs:
            line_list = line.split('\t') #split the pair
            fstem = line_list[14]
            rstem = line_list[15]
            nucl1 = line_list[5]
            nucl2 = line_list[7]
            flagIFcharNum = True if (len(nucl1.split('.')[3]) or len(nucl2.split('.')[3])) else False
            #construct dictionary dStemNucls = {stem : [[left wing][rigth wing]]}
            #fullstems:
            if fstem!= '\\N' :
                fstem = int(fstem)
                if(fstem not in dStemNucls):
                    dStemNucls[fstem]=[[nucl1],[nucl2]]
                else:
                    dStemNucls[fstem][0].append(nucl1)
                    dStemNucls[fstem][1] = [nucl2] + dStemNucls[fstem][1]
            #revstems:
            if rstem!= '\\N' :
                rstem = M + int(rstem)
                if(rstem not in dStemNucls):
                    dStemNucls[rstem]=[[nucl1],[nucl2]]
                else:
                    dStemNucls[rstem][0].append(nucl1)
                    dStemNucls[rstem][1].append(nucl2)
            Pairs.append((nucl1,nucl2))
    return dStemNucls, Pairs, flagIFcharNum

def overlap(a, b):
    #get sequence ab from sequences a and b without their intersection
    #in  : list of strings a and b
    #out : maximal amount of common symbols (maxn) and merged line
    maxn = -1
    for n in xrange(1, 1 + min(len(a), len(b))):
        suffix = a[-n:]
        prefix = b[:n]
        if prefix == suffix:
            maxn = n
    return maxn,a + b[maxn:]

def contains(a, b):
    #return bigger list (a or b) that contains another
    #in  : a and b lists
    #out : bigger list if it contains another, one of them if they are equal and -1 otherwise
    sa = str(a).strip('[').strip(']')
    sb = str(b).strip('[').strip(']')
    if(sa==sb): return a
    if len(sa)>len(sb):
        if sa.find(sb)>=0 :
            return a
    else:
        if sb.find(sa)>=0 :
            return b
    return -1

def overlap2side(seq1,seq2):
    #get minimal sequence (ab or ba) from sequences a and b without their intersection
    #in  : list of strings seq1 and seq2
    #out : minimal sequence from subsequences or -1 if there are no intersection
    result = contains(seq1,seq2)
    if result!=-1:
        return result
    nx, rx = overlap(seq1, seq2)
    ny, ry = overlap(seq2,seq1)
    if nx==-1 and ny==-1:
        result = -1 # no intersection
    else:
        result = rx if nx > ny else ry
    return result

def serial(a,b):
    #chech if b nucleotide follows nucleotide a
    #in  : list of strings (nucleotides) a and b
    #out : True / False
    nucl2 = b[0].split('.')
    nucl1 = a[-1].split('.')
    if nucl1[0] == nucl2[0]:
        # A.U.270.; A.G.271.
        if int(nucl2[2])-1==int(nucl1[2]) and not len(nucl1[3]) and not len(nucl2[3]):
            return True
        # A.U.270.Z; A.G.271.
        if int(nucl2[2])-1==int(nucl1[2]) and len(nucl1[3]) and not len(nucl2[3]) and nucl1[3]=='Z':
                return True
        if int(nucl2[2])==int(nucl1[2]) and len(nucl2[3]):
            #A.A.270.A, A.A.270.B
            if len(nucl1[3]):
                if (ord(nucl2[3])-ord(nucl1[3]))==1:
                    return True
            else:
                if(nucl2[3]=='A'):
                    return True
    return False



def serial2side(seq1,seq2, ):
    #chech if seq1 follows seq2 or otherwise
    #in  : list of strings (nucleotides) seq2 and seq2
    #out : whole sequence or -1, if sequences are not consecutive

    #check both sides
    if serial(seq1,seq2):
        return [seq1+seq2, (seq1[-1], seq2[0])]
    else:
        if serial(seq2,seq1):
            return [seq2+seq1, (seq2[-1], seq1[0])] #merge
        else: #are not consecutive
            return [-1, (0, 0)]

def getWings(dStemNucls,nWingsData):
    #return merged wings (fragments)
    #in  : dictionary of wings dStemNucls
    #out : list of Wings, WingsData
    lWings = []
    wingid = nWingsData+1
    for st in dStemNucls:
        # wingstr1 = base.aWing(wingid, modelid, multid, "", 0, st, dStemNucls[st], M)
        # wingstr2 = base.aWing(wingid+1, modelid, multid, "", 1, st, dStemNucls[st], M)
        # WingsData.append(wingstr1)
        # WingsData.append(wingstr2)
        lWings.append({'ID' : wingid, 'SEQ' : dStemNucls[st][0],'FRAGS': 0, 'TYPE': 0, 'ST' : st})
        lWings.append({'ID' : wingid+1, 'SEQ' : dStemNucls[st][1], 'FRAGS' : 0, 'TYPE': 1, 'ST' : st })
        wingid+=2

    return lWings

def sortByFirstNucl(wing):
        return int(wing['SEQ'][0].split('.')[2])

def mergeWings(lWings,Pairs,nFrag,WingsData,NuclsData,modelid,multid,M):
    #return merged wings (fragments)
    #in  : list of wings dStemNucls
    #out : list of dictionaries containing fragments

    lFragments = []
    lWingsTemp = []
    id = nFrag +1
    for wing in lWings:
        lWingsTemp.append({'SEQ': wing['SEQ'], 'WINGS': [wing['ID']], 'INDNUCLS': []})

    while True: #while there are fragments possible to merge
        for wing in lWingsTemp:
            flag = True #True - the wing is not consecutive and does not intersect
            if not len(lFragments):
                lFragments.append({'ID': -1, 'SEQ' : wing['SEQ'],'WINGS': wing['WINGS'], 'NWINGS': len(wing['WINGS']), 'CON' : [], 'INDNUCLS' : wing['INDNUCLS'] , 'BPSNUM' : 0 })
            else:
                for j in range(0, len(lFragments)):
                    frag = lFragments[j]
                    #check if they are consecutive
                    [seq,(nuc1,nuc2)] = serial2side(wing['SEQ'],frag['SEQ'])
                    if seq!=-1: # consecutive
                        flag = False
                        frag['WINGS']+=wing['WINGS']
                        frag['SEQ'] = seq #replace wing1 with merged fragment
                        frag['INDNUCLS']+=wing['INDNUCLS']
                        frag['INDNUCLS'].append((nuc1,nuc2))
                        lFragments[j] = frag
                        break
                    else: #not consecutive
                        #check if they intersect
                        inter = overlap2side(wing['SEQ'],frag['SEQ'])
                        if inter!=-1 : #intersect
                            flag = False
                            frag['WINGS']+=wing['WINGS']
                            frag['SEQ'] = inter
                            frag['INDNUCLS'] += wing['INDNUCLS']
                            if frag['INDNUCLS']:
                                for i in range(0, len(frag['INDNUCLS'])):
                                    nucls = frag['INDNUCLS'][i]
                                    if contains([nucls[0],nucls[1]], wing['SEQ']):
                                         frag['INDNUCLS'][i] = -1
                                frag['INDNUCLS']= [x for x in  frag['INDNUCLS'] if x != -1]
                            lFragments[j] = frag
                if flag:
                    lFragments.append({'ID': -1, 'SEQ' : wing['SEQ'],'WINGS' : wing['WINGS'], 'NWINGS': len(wing['WINGS']), 'CON' : [], 'INDNUCLS' : wing['INDNUCLS'], 'BPSNUM' : 0 })
        #for f in lFragments:
        #    print f
        if(len(lWingsTemp)> len(lFragments)):
            lWingsTemp = lFragments[:]
            lFragments = []
        else:
            break

    for fr in lFragments:
        fr['ID'] = id
        id+=1
        if fr['INDNUCLS']:
            fr['CON'] = len(fr['INDNUCLS'])
            fr['INDNUCLS'] = [fr['SEQ'].index(p[0]) for p in fr['INDNUCLS']]
        else:
            fr['CON'] = 0
        indexes = set()
        for nucl in fr['SEQ']:
            indexes = set( list(indexes) + filter(lambda i: Pairs[i][0] == nucl or Pairs[i][1] == nucl, xrange(0, len(Pairs))))
        fr['BPSNUM'] = len(indexes)

    for w in lWings:
        w['FRAGS'] = [frag['ID'] for frag in lFragments if w['ID'] in frag['WINGS']][0]
        wingstr = base.aWing(w['ID'], modelid, multid, w['FRAGS'], w['TYPE'], w['ST'], w['SEQ'], NuclsData[w['SEQ'][0]],NuclsData[w['SEQ'][-1]], M)
        WingsData.append(wingstr)
    return lFragments

#MAIN


def generateMults():
    start = time.time()
    M = 1000000 #if stems[i]>M then it is revStem and his id =stem[i]-M
    marked = set()
    stemGraph = {}
    dStemLinks = {}
    dNuclStems = {}
    dStemNucls = {}
    MultsData = []
    FragData = []
    WingsData = []
    EdgeData = []
    NuclsData = {}
    #file with models
    ModelsData = open('models.txt').readlines()
    for i in range(0,len(ModelsData)):
        modelstr = ModelsData[i].strip('\n')
        filepdb = modelstr.split('\t')[2]
        if filepdb in nrpdb.files:
            filepdb = 1
        else:
            filepdb = 0
        ModelsData[i] = modelstr + '\t'.join([str(el) for el in [filepdb, 0,0]]) + '\n'
    for line in open('nucls.txt').readlines():
        l = line.split('\t')
        model = int(l[2])
        if model not in NuclsData:
            NuclsData[model] = {}
        NuclsData[model][l[1]] = int(l[10])

    modelstr = "" #string for current model
    #directory to store files by models:
    models_dir = 'models'
    #directory to store result:
    mults_dir = 'mults'
    stop_line = '\n*****\n'

    #-----------FIND MULTIPLETS

    #1. Decompose the BasePairs file to the set of files containing pairs by model
    #e.g. 1.txt contains pairs of 1 model
    organizeFiles(models_dir, mults_dir)

    #work in 'models' directory
    filelist = [ f for f in os.listdir(models_dir)]

    statistics = {} #{num of vertices (containing stems): num of multiplets}
    statistics2 = {} #{num of vertices (containing fragments): num of multiplets}
    nTotalMults = 0
    idPairs = 1
    counter = 0
    counterCharNum = 0
    for f in filelist:
        curmodelnum = int(f.strip('.txt'))
        path = models_dir + '/' + f
        dNuclStems.clear()
        dStemLinks.clear()
        stemGraph.clear()
        inp = file(path, 'r')
        #read file and store content
        ModelPairs = inp.readlines()
        inp.close()

        #2. create dictionary {nucleotide: stems} for every model
        dNuclStems = nuclStems( ModelPairs, M )
        #print {nucleotide: stems}
        #print dNuclStems

        #3. create dictionary {stem: links} and adjacency list for every model
        #example
        #result {1: [100], 2: [3, 101], 3: [2, 101], 100: [1], 101: [2, 3]}
        stemGraph = stemLinks( ModelPairs, dStemLinks, dNuclStems)


        #4. find connected components
        marked.clear() #set of marked vertices
        lComponents = componentsFind(stemGraph, marked) #get the list of components
        #get statistics
        nMults = len(lComponents) #num of multiplets in the model
        nTotalMults += nMults
        for comp in lComponents:
            if len(comp) not in statistics:
                statistics[len(comp)] = 1
            else:
                statistics[len(comp)] += 1

        #5. Create output files
        nFrag = 0
        if len(lComponents)>0:
            ModelNucls = NuclsData[curmodelnum]
            nFrag,counterModel, counterCharNumModel = writeMults(f,M, dStemLinks, lComponents, ModelPairs,MultsData,FragData, WingsData,EdgeData,ModelNucls , mults_dir, stop_line, statistics2)
            counter+=counterModel
            counterCharNum+=counterCharNumModel
        idPairs = base.createBasePairs(ModelPairs, idPairs)
        ModelsData[curmodelnum-1] = ModelsData[curmodelnum-1].replace('\t'.join([str(el) for el in [0,0]]) + '\n', '\t'.join([str(el) for el in [nMults,nFrag]]) + '\n')

    print statistics
    print statistics2
    print counter
    print counterCharNum
    print 'Total = ',nTotalMults

    base.createModels(ModelsData)
    base.createMultiplets(MultsData)
    base.createFragments(FragData)
    base.createWings(WingsData)
    base.createEdges(EdgeData)


    #print taken time
    finish1 = time.time()
    d = int(finish1 - start)
    minutes, seconds = d/60, d % 60
    print minutes,'min ',seconds,'sec\n'

    return

generateMults()
