import matplotlib as mpl,numpy as np
mpl.use('Agg')
import matplotlib.pyplot as plt,cPickle
import seaborn as sns,pandas
from os.path import join,exists
from os import makedirs,system
scorecutoff = 30
topdir = '/cluster/zeng/code/research/lineage/analysis/error'
outdir = join(topdir,'figures_quality'+str(scorecutoff))
perbaseout = join(outdir,'perbase')
distrout = join(outdir,'distr')

if exists(outdir):
    system('rm -r ' + outdir)
makedirs(perbaseout)
makedirs(distrout)

plotrange = 30


with open(join(topdir,'error_quality'+str(scorecutoff)+'.pkl'),'rb') as f:
    result = cPickle.load(f)


def plotf(all_rate,outfile):
    a = np.asarray(all_rate).shape
    all_rate = pandas.DataFrame(all_rate,columns=('MutRate','BaseIdx','Condition'))
    plt.figure()
    ax = sns.barplot(x="BaseIdx", y="MutRate", hue="Condition", data=all_rate)
    fig = ax.get_figure()
    fig.savefig(outfile)
    plt.clf()

def boxplot(data,outfile):
    plt.figure()
    ax = sns.boxplot(data=data,x='MutRate',y='Condition',orient="h")
    fig = ax.get_figure()
    fig.tight_layout()
    fig.savefig(outfile)
    plt.clf()


def violinplot(data,outfile):
    plt.figure()
    ax = sns.violinplot(data=data,x='MutRate',y='Condition',orient="h")
    fig = ax.get_figure()
    fig.tight_layout()
    fig.savefig(outfile)
    plt.clf()

for refidx in range(2):
    t_r = result[refidx]
    boxdata = []
    for pair in range(2):

        ####
        l = len(t_r[0][pair])
        plotnum = int(np.ceil(float(l)/plotrange))
        all_rate = [[0]]*plotnum
        for conditiondata in t_r:
            for plotidx in range(plotnum):
                s = plotidx * plotrange
                e = np.min([l,s+plotrange])
                all_rate[plotidx] = all_rate[plotidx] + conditiondata[pair][s:e]

        for idx in range(plotnum):
            all_rate[idx] = all_rate[idx][1:]

        for plotidx in range(plotnum):
            plotf(all_rate[plotidx],join(perbaseout,'ref'+str(refidx)+'_'+'pair'+str(pair)+'_part'+str(plotidx)+'_barplot.png'))

        ###
        for condition in range(1,len(t_r)):
            toadd = t_r[condition][pair]
            for pos in range(len(toadd)):
                toadd[pos][0] = toadd[pos][0] - t_r[0][pair][pos][0]
            boxdata = boxdata + toadd
    boxdata = pandas.DataFrame(boxdata,columns=('MutRate','BaseIdx','Condition'))
    boxplot(boxdata,join(distrout,'ref'+str(refidx)+'_ratechange_boxplot.png'))
    violinplot(boxdata,join(distrout,'ref'+str(refidx)+'_ratechage_violinplot.png'))
