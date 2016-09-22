from os import system,makedirs
from os.path import join,exists
topdir = '/cluster/zeng/research/lineage'
for time in range(2):
    for rep in range(2):
	outdir = 'fth1_time{0}_rep{1}'.format(str(time+1),str(rep+1))
	if not exists(outdir):
	    makedirs(outdir)
	system(' '.join(['cp -r',join(topdir,'mESC_fth1_time{0}_rep{1}'.format(str(time+1),str(rep+1)),'post_split','cell*.csv'),outdir]))
	system(' '.join(['cp -r',join(topdir,'mESC_fth1_time{0}_rep{1}'.format(str(time+1),str(rep+1)),'post_split','cell*barcode*'),outdir]))
