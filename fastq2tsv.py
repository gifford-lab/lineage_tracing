from itertools import izip
from os.path import join
topdir = '/cluster/shortreads/reads/grna_deepseq_151019'
pair1 = join(topdir,'151019GifA_D15-8859_1_sequence.fastq')
pair2 = join(topdir,'151019GifA_D15-8859_2_sequence.fastq')

outdir = '/cluster/zeng/code/research/lineage/processed_data'
outfile = join(outdir,'concat.tsv')


barcode = ['aactc','ctgga','ggact','tctgc','aaccg','ctctg','ggtaa','aagct','tcgtc','ccaat','gcgta','tgagc']
barcodelen = len(barcode[0])
barcode = [x.upper() for x in barcode]
#revfirst = ['GGTGA','GGTGA','GGTGA','GGTGA','GGTGA','GGTGA','GCTCA','GCTCA','GCTCA','GCTCA','GCTCA','GCTCA']
banbarcode = ['CCCCT','CCCCT','ACATC','GCCTA','GCAGA','GCAGA']
revmap = {'A':'T','T':'A','C':'G','G':'C'}

#ooo = 0
#with open(pair1) as p1,open(pair2) as p2:
with open(pair1) as p1,open(pair2) as p2,open(outfile,'w') as out:
    cnt = 0
    for x ,y in izip(p1,p2):
        x = x.strip()
        y = y.strip()
        cnt = (cnt + 1) % 4
        if cnt == 2:
            flag = -1
            #if x[:5] == barcode[6] and y[:5] == revfirst[6]:
            #    ooo = ooo + 1
            for idx in range(len(barcode)):
                if x[:barcodelen] == barcode[idx]:
                    if idx<6 and y[:barcodelen] != banbarcode[idx] or idx >= 6:
                        seq_x = list(x)[barcodelen:]
                        seq_y = [revmap[t] for t in list(y)][::-1]
                        flag = idx
                        break

        if cnt == 0 and flag != -1:
            score_x = [str(ord(t)-33) for t in x][barcodelen:]
            score_y = [str(ord(t)-33) for t in y][::-1]
            out.write('%s\n' % '\t'.join([str(flag)] + seq_x + score_x + seq_y + score_y))

#print ooo


