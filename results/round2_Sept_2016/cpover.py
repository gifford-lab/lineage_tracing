from os import system,makedirs,listdir
from os.path import join,exists,isfile
topdir = '/cluster/zeng/research/lineage'
times = ['1','2','both']
reps = [['1','2'],['1','2'],['2']]
for idx,time in enumerate(times):
    for rep in reps[idx]:
	outdir = 'fth1_time{0}_rep{1}'.format(time,rep)
	if exists(outdir):
            system('rm -r '+outdir)
	makedirs(join(outdir,'fth1'))
        makedirs(join(outdir,'barcodelet'))
        system(' '.join(['cp -r',join(topdir,'mESC_fth1_time{0}_rep{1}'.format(time,rep),'post_split','cell*.csv'),outdir]))
	system(' '.join(['cp -r',join(topdir,'mESC_fth1_time{0}_rep{1}'.format(time,rep),'post_split','cell*.raw.tsv'),outdir]))
	system(' '.join(['cp -r',join(topdir,'mESC_fth1_time{0}_rep{1}'.format(time,rep),'post_split','cell*bcCoverage*'),outdir]))
        tfiles = [f for f in listdir(outdir) if isfile(join(outdir, f))]
        for x in tfiles:
            if 'barcodelet' in x:
                system(' '.join(['mv',join(outdir,x),join(outdir,'barcodelet', '.'.join(x.split('.')[2:]))]))
            else:
                system(' '.join(['mv',join(outdir,x),join(outdir,'fth1', '.'.join(x.split('.')[1:]))]))

outdir = 'stats'
if exists(outdir):
    system('rm -r '+outdir)
makedirs(outdir)
system(' '.join(['cp -r',join(topdir,'*stats*'),outdir]))
