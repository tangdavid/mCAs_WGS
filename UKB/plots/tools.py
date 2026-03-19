import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def downsample_gwas(df, log10p_col='LOG10P'):
    ps = 10**(np.ceil(df[log10p_col] * 2)/2) * 1e-2
    df = df[(np.random.uniform(size=len(ps)) < ps) | (df[log10p_col]>4)]
    return df

def manhattan_plot(df, ax, 
    c1='steelblue', c2='lightsteelblue', 
    log10p_col='LOG10P', pos_col='GENPOS', chrom_col='CHROM',
    max_p=300,
    rasterized=True,
):
    pos_adjustment = 0
    pos = []
    xticks = []
    xticklabels = []
    for chrom, group in df.groupby(chrom_col):
        ax.scatter(
            group[group[log10p_col]<max_p][pos_col] + pos_adjustment,  # Adjust positions by chromosome
            group[group[log10p_col]<max_p][log10p_col],
            c=[c1, c2][chrom % 2],
            marker='.',
            s=5,
            rasterized=rasterized
        )

        ax.scatter(
            group[group[log10p_col]>=max_p][pos_col] + pos_adjustment,  # Adjust positions by chromosome
            np.repeat(max_p, len(group[group[log10p_col]>=max_p])),
            c=[c1, c2][chrom % 2],
            marker='^',
            edgecolors='k',
            linewidths=0.3,
            s=5,
            rasterized=rasterized
        )
        pos+= ((pos_adjustment + group[pos_col]).to_list())
        xticks.append(pos_adjustment + group[pos_col].max()/2)
        xticklabels.append(chrom)
        pos_adjustment += group[pos_col].max() + 2e7
    ax.set_xticks(xticks, xticklabels, rotation=90)
    return pos

def plot_gtf(gtf, chrom, start, end, ax, highlight_genes={}, transpose=False, names=False, fontsize=18):
    plus_offset = {i:0 for i in range(1, 10)}
    minus_offset = {i:0 for i in range(0, 9)}
    gene2offset = {}
    for _,row in gtf.query('chr==@chrom and start<@end and end>@start and annot=="transcript"').sort_values(['start', 'end']).iterrows():
        for i in range(1, 10):
            if row['strand'] == '-':
                offset = minus_offset
                i -= 1
            else:
                offset = plus_offset 
            if offset[i] + 0.1e6 < row['start']:
                gene2offset[row['gene']] = 2*i if row['strand'] == '+' else -2*i
                offset[i] = row['end']
                break

    # Plot exons as blue rectangles
    for _, row in gtf\
            .query('chr==@chrom and start<@end and end>@start')\
            .merge(pd.DataFrame(gene2offset.items(), columns=['gene', 'y_offset']), on='gene')\
            .iterrows():
        y_pos = row['y_offset']
        color = 'red' if row['gene'] in highlight_genes else 'grey'

        item_start = max(row['start'], start)/1e6
        item_end = min(row['end'], end)/1e6

        if row['annot'] == 'transcript':
            if not transpose:
                ax.plot([item_start, item_end], [y_pos, y_pos], color="black")
                if names: 
                    text = f'\u25C0 {row["gene"]}' if row['strand'] == '-' else f'{row["gene"]} \u25B6'
                    if len(highlight_genes) != 0:
                        text = text if row['gene'] in highlight_genes else ''
                    ax.text(
                        (row['start']+row['end'])/2e6, y_pos + 0.5 if row['strand'] == '+' else y_pos - 1.5, 
                        text, 
                        fontstyle='italic', 
                        horizontalalignment='center',
                        color='black',
                        fontsize=fontsize
                    )
            else:
                ax.plot([y_pos, y_pos], [item_start, item_end], color="black")
            continue
        elif row['annot'] == 'exon':
            x, y, width, height = (item_start, y_pos - 0.5, item_end - item_start, 1)
        elif row['annot'] == '5UTR':
            x, y, width, height = (item_start, y_pos - 0.1, item_end - item_start, 0.2)
        elif row['annot'] == '3UTR':
            x, y, width, height = (item_start, y_pos - 0.1, item_end - item_start, 0.2)

        if not transpose:
            ax.add_patch(plt.Rectangle((x, y), width, height, color=color))
        else:
            ax.add_patch(plt.Rectangle((y, x), height, width, color=color))
    return ax

def plot_normalized_histogram(data, ax, color, alpha, bins, label, inverted=False, normalize=True):
    hist, bins = np.histogram(data, bins=bins)
    widths = np.diff(bins)
    if normalize: hist = hist/len(data)
    if inverted: hist = -hist
    ax.bar(bins[:-1], hist, widths, color=color, alpha=alpha, label=label, edgecolor='k', align="edge")
    return ax