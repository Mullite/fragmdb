
#Base Module

def createBasePairs( BasePairs, id):
    path = 'BasePairs_fragmdb.txt'
    if id==1:
        op = open(path , 'w+')
    else:
        op = open(path , 'a')
    for j in xrange(0,len(BasePairs)):
        line = aBasePair(BasePairs[j], id)
        op.write('\t'.join([str(line[i]) for i in [ 'ID', 'MODEL', 'MULTIPLET',  'NUM',   'NUCL1',  'BOND', 'NUCL2', 'TYPE', 'CLASS1', 'CLASS2','CLASS3','FULLSTEM', 'REVSTEM', 'LINK','EDGE', 'FRAG1', 'FRAG2','INFO1','INFO2','DIST1','DIST2','DIST3','TOR','HBONDS','HBNUM',  'SHEAR', 'STRETCH','STAGGER','BUCKLE','PROPELLER','OPENING']]))
        op.write('\n')
        id+=1
    return id

def createModels( ModelsData):
    path = 'Models_fragmdb.txt'
    op = open(path, 'w+')
    for line in ModelsData :
        op.write(line)

def createMultiplets( MultsData):
    path = 'Multiplets_fragmdb.txt'
    op = open(path, 'w+')
    for line in MultsData :
        op.write('\t'.join([str(line[i]) for i in ['ID','MODEL', 'NUMVERT','NUMEDGES','NUMPAIRS','FULLSTEMS' ,'REVSTEMS' ,'LINKS' ,'INVARIANT']]))
        op.write('\n')

def createFragments( FragData):
    path = 'Fragments_fragmdb.txt'
    op = open(path, 'w+')
    for line in FragData :
        op.write('\t'.join([str(line[i]) for i in ['ID', 'MODEL', 'CHAIN', 'MULTIPLET','START1','END1','START2','END2','LEN','BPSNUM','CON','SEQ','NWINGS' ] ]))
        op.write('\t')
        op.write(','.join([str(n) for n in line['DETAILS']]))
        op.write('\t')
        if line['INDNUCLS'] :
            op.write(','.join([str(n) for n in line['INDNUCLS']]))
        else:
            op.write('\\N')
        op.write('\n')

def createWings( WingsData):
    path = 'Wings_fragmdb.txt'
    op = open(path, 'w+')
    for line in WingsData :
        op.write('\t'.join([str(line[i]) for i in ['ID', 'MODEL', 'CHAIN', 'MULTIPLET','FRAG', 'TYPE','FS', 'RS', 'LEN' ,'ANOTHER','SEQ', 'START1', 'END1', 'START2', 'END2'] ]))
        op.write('\t')
        op.write(','.join([str(n) for n in line['DETAILS']]))
        op.write('\n')

def createEdges( EdgesData):
    path = 'Edges_fragmdb.txt'
    op = open(path, 'w+')
    for line in EdgesData :
        op.write('\t'.join([str(line[i]) for i in ['ID','MODEL','MULTIPLET','FRAG1','FRAG2','NUMREG','NUMIRREG'] ]))
        op.write('\n')

def aWing(id,model, multid, fragid, type, st, nucls,nucstart, nucend, M):
    return {
        'ID' : id,
        'MODEL': model,
        'CHAIN': nucls[0][0],
        'MULTIPLET': multid,
        'FRAG': fragid,
        'TYPE': 'R' if type else 'L',
        'FS': st if st < M else '\\N',
        'RS': st-M if st> M else '\\N',
        'LEN': len(nucls),
        'ANOTHER': id-1 if type else id+1,
        'SEQ': ','.join([nucl.split('.')[1] for nucl in nucls]),
        'START1': nucls[0],
        'END1':  nucls[-1],
        'START2': nucstart,
        'END2': nucend,
        'DETAILS': nucls
    }

def aFragment(id, model, multid, fragment, con,nwings,indnucls, bpsnum, nucstart, nucend):
    return {
        'ID' : id,
        'MODEL' : model,
        'CHAIN' : fragment[0][0],
        'MULTIPLET' :  multid,
        'START1' : fragment[0],
        'END1' : fragment[-1],
        'START2' : nucstart,
        'END2' : nucend,
        'LEN' : len(fragment),
        'BPSNUM' : bpsnum,
        'CON' : con,
        'SEQ' :  ','.join([nucl.split('.')[1] for nucl in fragment]),
        'DETAILS' : fragment,
        'NWINGS' : nwings,
        'INDNUCLS' : indnucls
    }

def aMultiplet(id, model, nFrag, nEdges,links, fstem, rstem, MultsPairs, invariant):
    return {
        'ID': id,
        'MODEL': model,
        'NUMVERT' : nFrag,
        'NUMEDGES' : nEdges,
        'NUMPAIRS' : len(MultsPairs),
        'FULLSTEMS' : ','.join([str(st) for st in fstem]) if fstem else '\\N',
        'REVSTEMS' : ','.join([str(st) for st in rstem])if rstem else '\\N',
        'LINKS' : ','.join([str(l) for l in links]) if links else '\\N',
        'INVARIANT' : invariant,
        'PAIRS' : MultsPairs
    }

def aEdge(id, model, multid, frag1, frag2, bond):
    nReg = 0
    nIreg = 0
    for b in bond:
        if b == 'WC' or b == 'WB':
            nReg +=1
        else:
            nIreg +=1
    return {
        'ID': id,
        'MODEL': model,
        'MULTIPLET': multid,
        'FRAG1': frag1,
        'FRAG2' : frag2,
        'NUMREG' : nReg,
        'NUMIRREG' : nIreg
    }

def aBasePair(BasePair,id):
    l = BasePair.split('\t')
    flag = False if len(l) == 41 else True
    return {
        'ID': id,
        'MODEL': l[1],
        'MULTIPLET': l[41] if flag else '\\N',
        'NUM': l[2],
        'NUCL1': l[5],
        'BOND': l[6],
        'NUCL2': l[7],
        'TYPE': l[8],
        'CLASS1': l[9],
        'CLASS2': l[10],
        'CLASS3': l[11],
        'FULLSTEM': l[14],
        'REVSTEM': l[15],
        'LINK': l[17],
        'EDGE': l[42] if flag else '\\N',
        'FRAG1': l[43] if flag else '\\N',
        'FRAG2': l[44] if flag else '\\N',
        'INFO1': l[20],
        'INFO2': l[21],
        'DIST1': l[22],
        'DIST2': l[23],
        'DIST3': l[24],
        'TOR': l[25],
        'HBONDS': l[26],
        'HBNUM': l[27],
        'SHEAR': l[28],
        'STRETCH': l[29],
        'STAGGER': l[30],
        'BUCKLE': l[31],
        'PROPELLER': l[32],
        'OPENING': l[33]
    }
