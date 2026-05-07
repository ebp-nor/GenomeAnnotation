import pysam
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import sys,argparse
import pyranges as pr
import pyfaidx

def get_coverage_bed(bed_file, chromosomes, type, window_size=100000):
    """
    Parsing a BED file
    """
    bed = pr.read_bed(bed_file)
    density_data = []

    for i in chromosomes.keys():
        if 'SUPER' not in i:
            continue
        chr_len = len(chromosomes[i])
        #print (i, chr_len)
        for start in range(0, chr_len, window_size):
            end = min(start + window_size, chr_len)
            interval = bed[i, start:end]
            chr = i.replace('SUPER_', '')
            if interval:
                #print(interval.summary(return_df=True, to_stdout=False)['pyrange']['sum'])
                sum = interval.summary(return_df=True, to_stdout=False)['pyrange']['sum']
                if sum > window_size:
                   density = 1.0
                else:
                   density = sum/window_size
                #remove SUPER_ substring
                density_data.append({'Chromosome': chr, 'Start': start, 'End': end, 'Density': density, 'Type': type})
            else:
                density_data.append({'Chromosome': chr, 'Start': start, 'End': end, 'Density': 0.0, 'Type': type})
            #print(bed[i, start:end])

    return pd.DataFrame(density_data)


def get_coverage_gff(gff_file, chromosomes, type, window_size=100000):
    """
    Parsing a GFF3 file
    """
    gff = pr.read_gff3(gff_file)
    density_data = []

    if type == "exons":
        selector  = (gff.Feature == 'exon')
        exons = gff [selector]
        gff = exons

    for i in chromosomes.keys():
        if 'SUPER' not in i:
            continue
        chr_len = len(chromosomes[i])
        print (i, chr_len)
        for start in range(0, chr_len, window_size):
            end = min(start + window_size, chr_len)
            interval = gff[i, start:end]
            chr = i.replace('SUPER_', '')
            if interval:
                #print(interval.summary(return_df=True, to_stdout=False)['pyrange']['sum'])
                sum = interval.summary(return_df=True, to_stdout=False)['pyrange']['sum']
                if sum > window_size:
                   density = 1.0
                else:
                   density = sum/window_size
                #remove SUPER_ substring
                density_data.append({'Chromosome': chr, 'Start': start, 'End': end, 'Density': density, 'Type': type})
            else:
                density_data.append({'Chromosome': chr, 'Start': start, 'End': end, 'Density': 0.0, 'Type': type})
            #print(bed[i, start:end])

    return pd.DataFrame(density_data)

# Example usage
# coverage_df = get_coverage("your_alignment.bam", "chr1", window_size=10000)

def plot_coverage(cov_df, title="Genome coverage plot"):
    """
    Plots genome coverage using Seaborn.
    """
    plt.figure(figsize=(12, 6))
    sns.lineplot(data=cov_df, x='Start', y='Coverage mean', hue='Sequencing technology', errorbar=None)
    #sns.lineplot(data=cov_pb_df, x='Start', y='Coverage', hue='PacBio', errorbar=None)
    plt.title(title)
    plt.ylim(0,75)
    plt.xlabel("Genomic position (bp)")
    plt.ylabel("Average coverage depth")
    #plt.show()
    return plt


def plot_density(den_df, title="Density"):
    """
    Plots genome density using Seaborn.
    """
    #plt.figure(figsize=(24, 12))
    plt.figure()
    g = sns.FacetGrid(den_df, col="Chromosome", hue='Type', col_wrap=3, height=4, aspect=2)
    g.map(sns.lineplot,'Start', 'Density', linewidth=0.5, alpha=0.8)
    g.set_axis_labels("Genomic position (bp)", "Density")
    g.set_titles(col_template="Chromosome {col_name}")

    g.add_legend()
    #sns.lineplot(data=den_df, x='Start', y='Density', hue='Chromosome', errorbar=None)
    #sns.lineplot(data=cov_pb_df, x='Start', y='Coverage', hue='PacBio', errorbar=None)
    #plt.title(title)
    #plt.ylim(0,75)
    #plt.xlabel("Genomic position (bp)")
    #plt.ylabel("Density")
    #plt.show()

    return plt

# Example usage (assuming coverage_df is already created)
# plot_coverage(coverage_df, title="Coverage on Chromosome 1")

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='')
    #parser.add_argument('-o', '--ont', help='ONT bam file', required=True)
    #parser.add_argument('-o', '--ont', help='ONT depth file (samtools depth -aa)', required=True)
    parser.add_argument('-l', '--longdust', help='longdust output', required=True)
    parser.add_argument('-f', '--fasta', help='assembly fasta file', required=True)
    #parser.add_argument('-p', '--pval', action='store', help='minimum p-value (log2) (not used presently)', required=False, type=int, default=-1)
    #parser.add_argument('-a', '--asm_size', help='Size of the genome assembly (bp)', type=int, default=0, required = False)
    #parser.add_argument('-s', '--score', help='Minimum score', type=float, default=12.0, required = False)
    #parser.add_argument('bed', action='store', help='ULTRA BED file', type=argparse.FileType('r'))
    parser.add_argument('-w', '--window', action='store', help='window size', type=int, default=100000)
    parser.add_argument('-u', '--ultra', help='ULTRA tsv output', required=True)
    parser.add_argument('-t', '--trash', help='TRASH2 GFF output', required=True)
    parser.add_argument('-g', '--genes', help='genes as GFF3 output', required=True)
    parser.add_argument('-e', '--edta', help='transposable elements (EDTA) as GFF3 output', required=True)

    args = parser.parse_args()
    pyf = pyfaidx.Fasta(args.fasta)

    #header = pysam.AlignmentFile(args.ont, "r").header.to_dict()
    #print (header)
    #chrs = header['SQ']
    #print (chrs) 
    print ("Processing longdust track")
    longdust = get_coverage_bed(args.longdust, pyf, "longdust", args.window)
    print ("longdust density mean %g  and std %g" % (longdust['Density'].mean(), longdust['Density'].std()))
    print ("Processing ULTRA track")
    ultra = get_coverage_bed(args.ultra, pyf, "ULTRA", args.window)
    print ("ultra density mean %g  and std %g" % (ultra['Density'].mean(), ultra['Density'].std()))
    print ("Processing TRASH2 track")
    trash = get_coverage_gff(args.trash, pyf, "TRASH", args.window)
    print ("TRASH2 density mean %g  and std %g" % (trash['Density'].mean(), trash['Density'].std()))
    print ("Processing genes")
    genes = get_coverage_gff(args.genes, pyf, "exons", args.window)
    print ("exon density mean %g  and std %g" % (genes['Density'].mean(), genes['Density'].std()))
    print ("Processing transposable elements")
    edta = get_coverage_gff(args.edta, pyf, "transposable elements", args.window)
    print ("transposable element density mean %g  and std %g" % (edta['Density'].mean(), edta['Density'].std()))

    #longdust_20k = get_coverage_bed(args.longdust, pyf, "longdust", 20000)
    #longdust_50k = get_coverage_bed(args.longdust, pyf, "longdust", 50000)
    #longdust_100k = get_coverage_bed(args.longdust, pyf, "longdust", 100000)

    tracks = pd.concat([longdust, ultra, trash, genes, edta], ignore_index=True)

    #longdust.to_csv('longdust.csv', index=False) 

    #need to parse the bed/gff files
    #for i in pyf.keys():
        #print("%s is %d bp long" % (i, len(pyf[i])))
        #if 'SUPER' in i['SN']:
            #print ("sequence %s is %s long" % (i['SN'], i['LN'])) 
            #cov_ont = get_coverage_bam(args.ont, i['SN'], i['LN'], 'ONT', 10000)
            #cov_pb = get_coverage_bam(args.pacbio, i['SN'], i['LN'], 'PacBio', 10000) 
            #cov_both =  pd.concat([cov_ont, cov_pb], ignore_index=True)
            #print(cov_both)
    """
    plt_10k = plot_density(longdust_10k)
    plt.savefig("longdust_10k.png")
    plt_10k.close()

    plt_20k = plot_density(longdust_20k)
    plt.savefig("longdust_20k.png")
    plt_20k.close()

    plt_50k = plot_density(longdust_50k)
    plt.savefig("longdust_50k.png")
    plt_50k.close()

    plt_100k = plot_density(longdust_100k)
    plt.savefig("longdust_100k.png")
    plt_100k.close()
    """
    plt = plot_density(tracks, '')
    plt.savefig("tracks.png")
    plt.savefig("tracks_600.png", dpi=600.0)
    plt.savefig("tracks.svg")
    plt.close()

    plt = plot_density(ultra, '')
    plt.savefig("ultra.png")
    plt.savefig("ultra_600.png", dpi=600.0)
    plt.savefig("ultra.svg")
    plt.close()

    plt = plot_density(longdust, '')
    plt.savefig("longdust.png")
    plt.savefig("longdust_600.png", dpi=600.0)
    plt.savefig("longdust.svg")
    plt.close()

    plt = plot_density(trash, '')
    plt.savefig("trash.png")
    plt.savefig("trash_600.png", dpi=600.0)
    plt.savefig("trash.svg")
    plt.close()

    plt = plot_density(genes, '')
    plt.savefig("genes.png")
    plt.savefig("genes_600.png", dpi=600.0)
    plt.savefig("genes.svg")
    plt.close()

    plt = plot_density(edta, '')
    plt.savefig("tes.png")
    plt.savefig("tes_600.png", dpi=600.0)
    plt.savefig("tes.svg")
    plt.close()
