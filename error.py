import multiprocessing as mp,cPickle
from os.path import join,exists
from os import makedirs

topdir = '/cluster/zeng/code/research/lineage/processed_data'
outdir = '/cluster/zeng/code/research/lineage/analysis/error/'
readlen = 150
barcodelen = 5
readlen1 = readlen - barcodelen
readlen2 = readlen
scorecutoff = 30
testnum = 6

def slave(args):
    topdir = args[0]
    i = args[1]
    t_ref1 = args[2]
    t_ref2 = args[3]
    scorecutoff = args[4]
    barcodelen = args[5]

    readlen1 = len(t_ref1)
    readlen2 = len(t_ref2)
    with open(join(topdir,'class'+str(i)+'.tsv')) as f:
        rate1 = [0]*readlen1
        test1 = [0]*readlen1
        rate2 = [0]*readlen2
        test2 = [0]*readlen2
        for x in f:
            line = x.strip().split()[1:]
            seq1 = line[:readlen1]
            score1 = [int(item) for item in line[readlen1:(2*readlen1)]]
            seq2 = line[(2*readlen1):(2*readlen1+readlen2)]
            score2 = [int(item) for item in line[(2*readlen1+readlen2):(2*readlen1+2*readlen2)]]
            for idx in range(readlen1):
                if score1[idx] > scorecutoff:
                    test1[idx] =  test1[idx] + 1
                    if seq1[idx] != t_ref1[idx]:
                        rate1[idx] = rate1[idx] + 1
            for idx in range(readlen2):
                if score2[idx] > scorecutoff:
                    test2[idx] = test2[idx] + 1
                    if seq2[idx] != t_ref2[idx]:
                        rate2[idx] = rate2[idx] + 1
    return ( [[rate1[x]/float(test1[x]),x+1,i+1] for x in range(readlen1)],[[rate2[x]/float(test2[x]),x+1,i+1] for x in range(readlen2)])

with open('refs') as f:
    refs = [list(x.strip()) for x in f]

pool = mp.Pool(processes=testnum)
result = [0]*2
for refidx in range(2):
    ref = refs[refidx]
    t_ref1 = ref[:readlen1]
    t_ref2 = ref[-readlen2:]
    args = []
    for i in range(refidx*testnum,(refidx+1)*testnum):
        args.append([topdir,i,t_ref1,t_ref2,scorecutoff,barcodelen])
    #slave(args[0])
    result[refidx] = pool.map(slave,args)

if not exists(outdir):
	makedirs(outdir)
with open(join(outdir,'error_quality'+str(scorecutoff)+'.pkl'),'wb') as f:
    cPickle.dump(result,f)












