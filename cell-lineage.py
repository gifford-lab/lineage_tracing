from os import listdir
from os.path import isfile,join,exists
from collections import Counter
import sys

celltype = sys.argv[1] #mESC or endoderm
analysis = sys.argv[2] #sys or let (for tranditional and barcodelet)
topdir0  = join('/cluster/zeng/research/lineage/',celltype)
topdir = join(topdir0,'post_split/filtered_fastq')
quality_thresh = 20    ### Any bases with quality below this cutoff will be changed to N
clust_mismatch = 2     ### We allow this amount of mistach between any pair of barcode in a cluster
n_thresh = 0.4         ### We ignore all barcode (After cosnodilation) that has more than this portion of N
prefix_suffix_tol = 2  ### The mismatch we tolerant when identifying putative pe2 barcode
nfilter = 13 ### We ingore all barcodes (before consodilation) that have more than this number of N
beforebarcode = 'GAAAGGATGGGAGTACTAAGCT' if analysis == 'sys' else 'AGGCCTTTCGACCTGCATCCA'
barcodelen = 14
afterbarcode = 'TTTCTATAAGT' if analysis == 'sys' else  'AAAAAAAAAAA'
#afterbarcode = 'TGGATGCAGGT'
exptcode = '' if analysis == 'sys' else '.barcodelet'
if not exists(topdir):
    print 'Run preprocess and split first!'
    sys.exit(1)

outputfile = join(topdir0,'post_split','cell-lineage_mapping'+'_quality'+str(quality_thresh)+'_mismatch'+str(clust_mismatch)) + '_prefixsuffixtol'+str(prefix_suffix_tol) + exptcode

### Functions

def indexing(seq,key,equal):
    ### Find the element in the list that matches a key
    if equal:
        return set([i for i in xrange(len(seq)) if seq[i]==key])
    else:
        return set([i for i in xrange(len(seq)) if seq[i]!=key])

def simpledist(seq1,seq2):
    return sum([1 for idx in range(len(seq1)) if seq1[idx] != seq2[idx]])

def dist(seq1,seq2):
    ### Calculate the  Levenshtein distance distance of two string, considering 'N' as wildcard
    n1 = indexing(seq1,'N',False)
    n2 = indexing(seq2,'N',False)
    nindex = list(set.intersection(n1,n2))
    d = sum([1 for x in nindex if seq1[x]!=seq2[x]])
    return d

def Most_Common(lst,cnts):
    t_list = []
    for idx in range(len(lst)):
        if lst[idx]!='N':
            t_list+= [lst[idx]]*cnts[idx]
    if len(t_list)>0:
        data = Counter(t_list)
        return data.most_common(1)[0][0]
    else:
        return 'N'

def consolidate(seqs,cnts):
    ### Consolidate a set of aligned seqs into one
    out = list(seqs[0])
    for x in range(len(out)):
        out[x] = Most_Common([t[x] for t in seqs],cnts)
    return (''.join(out),sum(cnts))

def checkforN(seqs,thresh):
    out = [idx for idx in range(len(seqs)) if list(seqs[idx]).count('N')< thresh*len(seqs[idx])]
    return out

def align(seqs,seq_cnts):
    ### Align the candidate barcode to find the unique ones based on quality score
    out = []
    out_cnt = []
    for seqidx in range(len(seqs)):
        seq = seqs[seqidx]
        seq_cnt = seq_cnts[seqidx]
        if seq.count('N')>nfilter:
            continue;
        if len(out)==0:
            out.append([seq])
            out_cnt.append([seq_cnt])
        else:
            found = -1
            closest = 10000
            for clst_idx in xrange(len(out)):
                clst = out[clst_idx]
                flag = True
                idx = 0
                while flag and idx < len(clst):
                    t_dist = dist(seq,clst[idx])
                    flag = t_dist<=clust_mismatch
                    idx += 1
                if flag and t_dist < closest:
                    found = clst_idx
                    closest = t_dist
                    break
            if found == -1:
                out.append([seq])
                out_cnt.append([seq_cnt])
            else:
                out[found].append(seq)
                out_cnt[found].append(seq_cnt)

    print 'Before consolidate:'
    print out
    print out_cnt
    for idx in xrange(len(out)):
        out[idx],out_cnt[idx] = consolidate(out[idx],out_cnt[idx])

    goodones = checkforN(out,n_thresh)
    out_cnt = [out_cnt[x] for x in goodones]
    out = [out[x] for x in goodones]
    if len(out) == 0:
        out = ['/']
        out_cnt = [0]
    print 'After consolidate:'
    print out
    print out_cnt
    return (out,out_cnt)

#############

expts = [f for f in listdir(topdir) if isfile(join(topdir, f))]

lineage = {}
output = []
for expt in expts:
    cell = expt.split('.')[0]
    if cell == 'unassigned':
        continue

    print 'checking',cell
    barcode_cnt = {}
    umis = []
    with open(join(topdir,expt)) as f:
        cnt = 0
        flag = False
        umicheck = True
        readcnt = 0
        pcr_read_cnt = 0
        for x in f:
            x = x.strip()
            cnt = (cnt + 1)%4
            if cnt == 0 and flag and umicheck:
                quality = [ord(m)-33 for m in x[22:39]]
                for m in xrange(len(barcode)):
                    if quality[m] < quality_thresh:
                        barcode[m] = 'N'
                barcode = "".join(barcode)

                if barcode in barcode_cnt.keys():
                    barcode_cnt[barcode] += 1
                else:
                    barcode_cnt[barcode] = 1
                    if cell in lineage.keys():
                        lineage[cell].append(barcode)
                    else:
                        lineage[cell] = [barcode]
                flag = False
            if cnt == 2 and umicheck:
                readcnt += 1
                d1= simpledist(x[:len(beforebarcode)],beforebarcode)
                if d1 <= prefix_suffix_tol:
                    #d2 = simpledist(x[-len(afterbarcode):] ,afterbarcode)
                    d2 = 0
                    d2 = simpledist(x[(len(beforebarcode)+barcodelen):(len(beforebarcode)+barcodelen+len(afterbarcode))],afterbarcode)
                    if d1 + d2<=prefix_suffix_tol:
                        barcode = list(x[len(beforebarcode):(len(beforebarcode)+barcodelen)])
                        flag = True
            if cnt == 1:
                pcr_read_cnt += 1
                t_umi = x.split(':')[-2]
                if t_umi in umis:
                    umicheck = False
                else:
                    umis.append(t_umi)
                    umicheck = True

    print 'start aligning'
    alignment,align_cnt = align(lineage[cell],[barcode_cnt[x] for x in lineage[cell]]) if cell in lineage.keys() else (['/'],[0])
    toadd =[cell,readcnt,sum(barcode_cnt.values()),','.join(alignment),','.join(map(str,align_cnt)),len(indexing(alignment,'/',False)),pcr_read_cnt]
    print 'Final:'
    print toadd
    output.append(toadd)

with open(outputfile,'w') as f:
    for x in output:
        f.write('%s\n' % ('\t'.join(map(str,x))))

