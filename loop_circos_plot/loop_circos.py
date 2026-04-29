import fire
import os
import pandas as pd
import pycircos
import seaborn as sns
import matplotlib.pyplot as plt
import patchworklib as pw
from Bio import SeqIO
import bioframe as bf
import matplotlib
from matplotlib import cm
from  multiprocessing import Pool
import multiprocessing
import matplotlib.colors as mcolors
Garc    = pycircos.Garc
Gcircle = pycircos.Gcircle
plt.rcParams['font.size'] = 30



def int_to_roman(num):
    roman_numerals = {
        1: 'I', 4: 'IV', 5: 'V', 9: 'IX', 10: 'X', 40: 'XL', 50: 'L',
        90: 'XC', 100: 'C', 400: 'CD', 500: 'D', 900: 'CM', 1000: 'M'
    }
    roman = ''
    num = int(num[3:])
    for value, numeral in sorted(roman_numerals.items(), reverse=True):
        while num >= value:
            roman += numeral
            num -= value
    return roman


def main(loop_bedpe,chromsize,gene_cov,gc_cov,res,outprefix):
    colors = ["gray", "orange"] 
    cmap = mcolors.LinearSegmentedColormap.from_list("gray_to_orange", colors)
    norm = mcolors.Normalize(vmin=5, vmax=100)

    link = pd.read_csv(loop_bedpe,sep='\t')#,usecols=[0,1,2,3,4,5,6])
    link['color'] = link['score'].apply(lambda x: cmap(norm(x)))
    Chromsize = pd.read_csv(chromsize,sep='\t',header=None)
    Chromsize[0] = Chromsize[0].apply(int_to_roman)
    Chromsize = Chromsize.set_index(0)

    chrlist = Chromsize.index
    gene_cov = pd.read_csv(gene_cov,sep='\t')
    gc_cov = pd.read_csv(gc_cov,sep='\t') 
    gc_cov['chrom'] = gc_cov['chrom'].apply(int_to_roman)
    gene_cov['chrom'] = gene_cov['chrom'].apply(int_to_roman)
    link['chrom1'] = link['chrom1'].apply(int_to_roman)
    link['chrom2'] = link['chrom2'].apply(int_to_roman)

    run_all(Chromsize,gene_cov,gc_cov,link,outprefix)



def run_all(Chromsize,gene_cov,gc_cov,link,outprefix):
    circls = Gcircle(figsize=(7,7))
    add_chromsomes(circls,Chromsize)
    add_link(circls,link,935)
    circls.figure.savefig(f'{outprefix}_genome.png',bbox_inches='tight',dpi=300)
    circls.figure.savefig(f'{outprefix}_genome.pdf',bbox_inches='tight')
    plt.clf()

def run_chr(Chromsize,gene_cov,gc_cov,link,outprefix,chrom):
    circls = Gcircle(figsize=(10, 10))
    add_chromsomes(circls, Chromsize,chrom=chrom)
    add_barplot(circls, gene_cov, (845,925),chrom)
    add_link(circls, link, 835,chrom)
    circls.figure.savefig(f'{outprefix}_{chrom}.png',bbox_inches='tight',dpi=300)
    circls.figure.savefig(f'{outprefix}_{chrom}.pdf',bbox_inches='tight')
    plt.clf()
    plt.cla()


def add_chromsomes(circle,chromsize,chrom='all',label_visible=True,labelposition=80,raxis_range=(955,985)):

    if chrom =='all':
        for chr,length in zip(chromsize.index,chromsize.loc[:,1]):
            arc = Garc(arc_id = chr,size=length,facecolor='white',interspace=2, raxis_range=raxis_range,labelposition=labelposition, label_visible=label_visible,labelsize=30)
            circle.add_garc(arc)
    else:
        chrDict = chromsize[1].to_dict()
        length = chrDict[chrom]
        arc = Garc(arc_id=chrom, size=length, interspace=2, raxis_range=raxis_range,labelposition=labelposition, label_visible=label_visible )
        circle.add_garc(arc)

    circle.set_garcs(start=5,end=355)
    for arc_id in circle.garc_dict:
        circle.tickplot(arc_id,raxis_range=(985,1000), tickinterval=100000, ticklabels=None)

def add_link(circle,link,raxis_range,chrom='all'):
    if chrom == 'all':
        tmplink = link.copy()
    else:
        tmplink = link[link['chrom1']==chrom]
    if tmplink.shape[0] >0:
        for i in tmplink.to_numpy():
            c1,s1,e1,c2,s2,e2,fdr,color = i
            source = (c1,s1-1,e1,raxis_range)
            destination = (c2,s2-1,e2,raxis_range)
            circle.chord_plot(source,destination,facecolor=color)


def add_barplot(circle,df,raxis_range,chrom='all',vmin=None,vmax=None):
    if vmin is  None:
        vmin = df['value'].min()
    if vmax is None:
        vmax = df['value'].max()
    df['width'] = df['end'] - df['start']
    if chrom =='all':
        chrlist = set(df['chrom'].to_list())
        for chr in chrlist:
            tdf = df[df['chrom'] ==chr]
            circle.barplot(chr,data=tdf['value'],positions=tdf['start'],width=tdf['width'],raxis_range=raxis_range,base_value=0,rlim=[vmin-0.05*abs(vmin), vmax+0.05*abs(vmax)],facecolor='r',spine=True)
    else:
        tdf = df[df['chrom'] ==chrom]
        circle.barplot(chrom, tdf['value'], positions=tdf['start'], width=tdf['width'], raxis_range=raxis_range,
                               base_value=0, rlim=[vmin , vmax ], facecolor='r',
                               spine=True)

def add_heatmap(circle,df,raxis_range,chrom='all',vmin=None,vmax=None,color=plt.cm.Reds):
    # df chrom start end value
    if vmin is  None:
        vmin = df['value'].min()
    if vmax is None:
        vmax = df['value'].max()
    df['width'] = df['end'] - df['start']
    if chrom =='all':
        chrlist = set(df['chrom'].to_list())
        for chr in chrlist:
            tdf = df[df['chrom'] ==chr]
            print(tdf)
            circle.heatmap(chr,data=tdf['value'],positions=tdf['start'],width=tdf['width'],raxis_range=raxis_range,vmin=vmin,vmax=vmax,cmap=color)
    else:
        tdf = df[df['chrom'] ==chrom]
        circle.heatmap(chrom, tdf['value'], positions=tdf['start'], width=tdf['width'],raxis_range=raxis_range,vmin=vmin,vmax=vmax,cmap=color)



def get_bins(chromsize,res):
    df = pd.read_csv(chromsize,sep='\t',header=None)
    df[2] = df[1].apply(lambda x: list(range(0, x,res))+[x])
    bins = pd.DataFrame([[x]+[y1]+[y2] for x,y in zip(df[0],df[2]) for y1,y2 in zip(y[:-1],y[1:]) ])
    bins.columns=['chrom','start','end']
    return bins

if __name__ == '__main__':
    fire.Fire(main)
