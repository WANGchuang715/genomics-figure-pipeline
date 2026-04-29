import fire
import cooltools
import numpy as np
import seaborn as sns
import pandas as pd
from os.path import  *
import  matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import cooler
import os
import sys
from matplotlib.ticker import EngFormatter
import cooltools.lib.plotting
from multiprocessing import     Pool
BIN = abspath(dirname(sys.argv[0]))
bp_formatter = EngFormatter('b')
plt.rcParams['font.size'] = 16




def main(zscore1,zscore2,chromsize,res,outprefix):
    out,prefix = abspath(dirname(outprefix)),basename(outprefix)
    sample1,sample2 = prefix.split('_vs_')
    plot_new_all(zscore1,zscore2,f'{outprefix}.diff_genome.{int(int(res)/1000)}kb',chromsize,res,f'{sample1} - {sample2} z-Score')


def int_to_roman(num):
    roman_numerals = {
        1: 'I', 4: 'IV', 5: 'V', 9: 'IX', 10: 'X', 40: 'XL', 50: 'L',
        90: 'XC', 100: 'C', 400: 'CD', 500: 'D', 900: 'CM', 1000: 'M'
    }
    roman = ''
    for value, numeral in sorted(roman_numerals.items(), reverse=True):
        while num >= value:
            roman += numeral
            num -= value
    return roman


def plot_new_all(zscore1,zscore2,outprefix,chromsize,res,label):
    cdf = pd.read_csv(chromsize,sep='\t',header=None,names=['chrom','length'])
    cdf['cumsum'] = cdf['length'].cumsum()
    cdf['bin'] = cdf['cumsum'].apply(lambda x: int(x//res)+1)
    cdf['new_chr'] = cdf['chrom'].apply(lambda x: int_to_roman(int(x[3:])))
    clist = [0]  + cdf['cumsum'].tolist()
    linexticks = cdf['cumsum'].tolist()[:-1]
    xticks = [ (clist[i]+clist[i+1])/2       for i in range(len(clist)-1)]
    chrlist = cdf['new_chr'].tolist()
    

    df1 = pd.read_csv(zscore1,sep='\t',compression='gzip',comment='#',index_col=0)
    df2 = pd.read_csv(zscore2,sep='\t',compression='gzip',comment='#',index_col=0)
    arr1 = df1.to_numpy()
    arr2 = df2.to_numpy()
    mat = arr1 - arr2
    f,ax = plt.subplots(figsize=(7,6))
    ticklabels_fontsize = 'small'
    if len(xticks) > 18:
        ticklabels_fontsize = 6
    ax.xaxis.set_major_formatter(bp_formatter)
    ax.set_yticks(xticks)
    ax.set_yticklabels(chrlist, fontsize=ticklabels_fontsize)
    ax.grid(False)
    for i in ['bottom', 'left', 'top', 'right']:
        ax.spines[i].set_color('black')
        ax.spines[i].set_linewidth(0.5)
    extent = (0,cdf['length'].sum(),cdf['length'].sum(),0)
    im = ax.imshow(mat, interpolation='none', cmap='seismic',vmin=-10,vmax=10,extent=extent)
    cbar = plt.colorbar(im ,fraction=0.046, pad=0.04)
    cbar.set_label(label, fontsize=12) 
    ax.set_title(basename(outprefix),y=1.08)
    plt.savefig(f'{outprefix}.png',dpi=300,bbox_inches='tight')
    plt.savefig(f'{outprefix}.pdf', bbox_inches='tight',dpi=300)
    plt.close()
    plt.cla()
    plt.clf()


def get_losser_zscore_matrix(clr,chr,outprefix):
    print(f'start to deal {outprefix}')
    mat = clr.matrix(balance=True).fetch(chr)
    bin = clr.bins().fetch(chr)
    bin['index']  = bin.index + 1
    header = bin.apply(lambda x: f'HiC__{x[4]}|aaa|{x[0]}:{x[1]}-{x[2]}',axis=1 ).to_list()
    df = pd.DataFrame(mat)
    df.columns = header
    df.index = header
    df.fillna(0,inplace=True)
    in_mat = f'{outprefix}_{chr}.matrix.gz'
    df.to_csv(in_mat,sep='\t',compression='gzip')
    CMD = f'perl -I {BIN}/lib {BIN}/matrix2loess.pl -i {in_mat} --ez --ca 0.005 -o {outprefix}_{chr}'
    os.system(CMD)




def format_ticks(ax, x=True, y=True, rotate=False):
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)


if __name__ == '__main__':
    fire.Fire(main)
