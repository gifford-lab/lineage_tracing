from os import system
from os.path import join
import sys,distance,multiprocessing as mp,pandas as pd,numpy as np,seaborn as sns

topdir = '/cluster/shortreads/reads/160908_lineagetracing_zeng/160902_NS500414_0227_HLLWYBGXY/Data/Intensities/BaseCalls'
outdir = 'index_analysis'

order = sys.argv[1].split('_')

index_mismatch = 0
myindex = ['ATCACG',
        'CGATGT',
        'TTAGGC',
        'TGACCA',
        'ACAGTG',
        'GCCAAT',
        'CAGATC',
        'ACTTGA']

if '0' in order:
    ### re-combine the "undetermined" and "determined" fastq
    cmd = ['cat','cat']
    for i in range(4):
        for j in range(2):
            posfile = join(topdir,'NGEI017/SHE1322A1/lineage-inDrops_HLLWYBGXY_S1_L00'+str(i+1)+'_R'+str(j+1)+'_001.fastq.gz')
            negfile = join(topdir,'Undetermined_HLLWYBGXY_S0_L00'+str(i+1)+'_R'+str(j+1)+'_001.fastq.gz')
            newfile = join(topdir,'all_'+str(i+1)+'_R'+str(j+1)+'_001.fastq.gz')
            system(' '.join(['cat',posfile,negfile,'>',newfile]))
            cmd[j] += ' ' + newfile
    for j in range(2):
        system(' '.join([cmd[j],'>',join(topdir,'all_R'+str(j+1)+'.fastq')]))

if '1' in order:
    ### get only the barcodes for each read
    for i in range(4):
        for j in range(1):
            system(' '.join(['zcat',join(topdir,'all_'+str(i+1)+'_R'+str(j+1)+'_001.fastq.gz'),'|','sed -n \'1~4p\'','|','awk -F\':\' \'{print $10}\'','>',join(outdir,'lane_'+str(i+1)+'_R'+str(j+1))]))

def mydist(args):
    query,mydict,thresh = args[:]
    mydistance = [ (idx,distance.hamming(query,mydict[idx])) for idx in range(len(mydict))]
    mydistance.sort(key = lambda x:x[1])
    return mydistance[0] if mydistance[0][1] != mydistance[1][1] else (-1,mydistance[0][1])

def indexfinder(infile,outfile,unmapfile):
    with open(infile) as f,open(outfile,'w') as fout,open(unmapfile,'w') as funmap:
        counter = [0]*(len(myindex)+1)
        args = [[x.strip(),myindex,index_mismatch] for x in f]
        pool = mp.Pool(processes=10)
        group = pool.map(mydist,args)
        pool.close()
        pool.join()
        for idx in range(len(group)):
            t_idx = group[idx][0]
            counter[t_idx]+=1
            if t_idx == -1:
                funmap.write('%s\n' % args[idx][0])
        for x in counter:
            fout.write('%d\n' % x)

if '2' in order:
    ### for each barcode, find its belonged group. Output to a file if its not uniquely mapped to any
    for i in range(4):
        for j in range(1):
            print i,j
            #indexfinder(join(outdir,'test'),'test','test_mapped')
            indexfinder(join(outdir,'lane_'+str(i+1)+'_R'+str(j+1)),join(outdir,'lane_'+str(i+1)+'_R'+str(j+1)+'_indexsummary'),join(outdir,'lane_'+str(i+1)+'_R'+str(j+1)+'_unmapped'))


if '3' in order:
    ### check the hamming distance between any pair of the desired index
    for i in range(len(myindex)):
        for j in range(i+1,len(myindex)):
            print distance.hamming(myindex[i],myindex[j])

if '4' in order:
    ### analysis the properties of the unmapped barcodes
    from collections import Counter
    def mydist1(query,mydist):
        return np.min([distance.hamming(query,x) for x in myindex])
    for i in range(1):
        for j in range(1):
            with open(join('/cluster/zeng/code/research/lineage/index_analysis',\
                                                           'lane_'+str(i+1)+'_R'+str(j+1)+'_unmapped')) as f:
                data = [x.strip() for x in f]
            t_key =  Counter(data).keys()
            t_value = Counter(data).values()
            mypd = pd.DataFrame([[t_key[idx],t_value[idx],i+1,mydist1(t_key[idx],mydist)]for idx in range(len(t_key))],\
                columns=['kmer','count','lane','minDist2TargetIdx']).sort(['count'],ascending=False)
            mypd.to_csv('unmapped.tsv')
            ax = sns.distplot(mypd['count'])
            ax.get_figure().savefig('unmapped.png')
