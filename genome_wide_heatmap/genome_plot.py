import fire
import cooltools
import numpy as np
import seaborn as sns
import pandas as pd
from os.path import  *
import  matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import cooltools.lib.plotting
import cooler
from matplotlib.ticker import EngFormatter
bp_formatter = EngFormatter('b')
plt.rcParams['font.size'] = 16


def main(cool,outprefix,chromsize,*res):
    out,prefix = abspath(dirname(outprefix)),basename(outprefix)
    for RES in res:
        clr = cooler.Cooler(f'{cool}::resolutions/{RES}')
        plot(clr,f'{out}/{prefix}-{int(RES/1000)}k_raw',False,chromsize,RES)
        plot(clr, f'{out}/{prefix}-{int(RES / 1000)}k_correct', True,chromsize,RES)

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

def plot(clr,outprefix,balance,chromsize,res):
    cdf = pd.read_csv(chromsize,sep='\t',header=None,names=['chrom','length'])
    cdf['cumsum'] = cdf['length'].cumsum()
    cdf['bin'] = cdf['cumsum'].apply(lambda x: int(x//res)+1)
    cdf['new_chr'] = cdf['chrom'].apply(lambda x: int_to_roman(int(x[3:])))
    
    clist = [0]  + cdf['cumsum'].tolist()
    linexticks = cdf['cumsum'].tolist()[:-1]
    xticks = [ (clist[i]+clist[i+1])/2       for i in range(len(clist)-1)] 

    chrlist = cdf['new_chr'].tolist()
    if balance:
        norm = LogNorm(vmax=0.1,vmin=0.0001)
    else:
        norm = LogNorm(vmin=1,vmax=10_000)

    if balance:
        matrix = clr.matrix(balance=True)[:]
    else:
        matrix = clr.matrix(balance=False)[:]
    fig = plt.figure(figsize=(7,6))
    ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
    ticklabels_fontsize = 'small'
    if len(xticks) > 18:
        ticklabels_fontsize = 6
    
    for i in ['bottom', 'left', 'top', 'right']:
        ax.spines[i].set_color('black')
        ax.spines[i].set_linewidth(0.5)
    extent = (0,cdf['length'].sum(),cdf['length'].sum(),0)

    im = plt.imshow(matrix, interpolation='none', cmap='fall',norm=norm,extent=extent)
    ax.xaxis.set_major_formatter(bp_formatter)
  
    ax.set_yticks(xticks)
    ax.set_yticklabels(chrlist, fontsize=ticklabels_fontsize)
    ax.grid(False)

    cbar = plt.colorbar(im ,fraction=0.046, pad=0.04)

    plt.savefig(f'{outprefix}.png',dpi=300,bbox_inches='tight')
    plt.savefig(f'{outprefix}.pdf', bbox_inches='tight')

if __name__ == '__main__':
    fire.Fire(main)
