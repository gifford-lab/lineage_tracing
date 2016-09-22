from os import system,makedirs,listdir
from os.path import join,exists,isfile
topdir = '/cluster/zeng/research/lineage'
for time in range(2):
    for rep in range(2):
	outdir = 'fth1_time{0}_rep{1}'.format(str(time+1),str(rep+1))
	if exists(outdir):
            system('rm -r '+outdir)
	makedirs(outdir)
	system(' '.join(['cp -r',join(topdir,'mESC_fth1_time{0}_rep{1}'.format(str(time+1),str(rep+1)),'post_split','cell*.csv'),outdir]))
	system(' '.join(['cp -r',join(topdir,'mESC_fth1_time{0}_rep{1}'.format(str(time+1),str(rep+1)),'post_split','cell*barcode*'),outdir]))
        tfiles = [f for f in listdir(outdir) if isfile(join(outdir, f))]
        for x in tfiles:
            system(' '.join(['mv',join(outdir,x),join(outdir,'.'.join(x.split('.')[1:]))]))
