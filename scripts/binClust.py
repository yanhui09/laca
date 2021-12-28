import sys
import pandas as pd
import argparse
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser()
    # Positional mandatory arguments
    parser.add_argument('paf', help='Alignment file produced by all-vs-all minimap2', type=str)
    # Optional arguments
    parser.add_argument('-p', '--prefix', help='Output file prefix [output]', type=str, default='output')
    parser.add_argument('-s', '--min_score_frac', help='Minimum fraction of maximum score for a cluster [0.8]', type=float, default=0.8)
    parser.add_argument('-n', '--min_reads', help='Minimum reads to create cluster [5]', type=int, default=5)
    parser.add_argument('-R', '--max_recursion', help='Override the max recursion depth in python', type=int, default=-1)
    args = parser.parse_args()
    return args

def run_clustering(df, args, default_sim):
    Z   = sch.linkage(squareform(df, checks=False), 'ward')
    args.max_dist = 3*np.median(Z[:,2])
    clusters = sch.fcluster(Z, args.max_dist, criterion='distance')
    
    df['cluster'] = clusters
    # exclude self comparison
    df.replace(0, np.nan, inplace=True)
    # mean score to others within group and across groups
    scores = df.groupby(['cluster']).mean()
    # extract value by numerical index    
    mean_scores = scores.values[[x-1 for x in clusters], range(scores.shape[1])]
    df_out = pd.DataFrame(data={'qname': df.index, 'clust_read_score': mean_scores, 'cluster': clusters})
    df_out.reset_index(drop=True, inplace=True)
    return df_out

def main(args):
    if args.max_recursion > 0:
        sys.setrecursionlimit(args.max_recursion)

    df = pd.DataFrame()              
    cols = ['qname', 'qlen', 'tname', 'tlen', 'alnlen', 'dv']
    dtypes = {0:"category", 1:"int16", 5:"category", 6:"int16", 10:"int16", 15:"category"}
    chunksize = 10 ** 6
    for chunk in pd.read_csv(args.paf, sep='\t', header=None, usecols=[0,1,5,6,10,15], names=cols, dtype=dtypes, chunksize=chunksize):
        if chunk.shape[0] > 0:
            chunk['sim'] = [1 - float(x.split(':')[-1]) for x in chunk['dv']]
            chunk['score'] = (chunk['alnlen']/10000)**2 * chunk['sim']
            chunk = chunk.drop(columns=['dv', 'alnlen'])
            df = pd.concat([df, chunk])
        else:
            # There are no alignments to parse
            # print('keep\tread\tall_v_all_aligns')
            # print('{}\t{}\t{}'.format('False', 'none', 0))
            pass
    # keep the longest in case of supplementary alignments, minimap sort this by default  
    df.drop_duplicates(['qname', 'tname'], inplace = True)
    #create similarity matrix out of pairwise alignment scores
    nan_value = 0
    # make unique, sorted, common index
    idx = sorted(set(df['qname']).union(df['tname']))
    # reshape
    df_pivot = (df.pivot(index='qname', columns='tname', values='score')
    .reindex(index=idx, columns=idx)
    .fillna(nan_value, downcast='infer')
    .pipe(lambda x: x+x.values.T)
    )
    # clustering     
    cluster_df = run_clustering(df_pivot, args, nan_value)
    
    # rls dict
    rls = dict(zip(df['qname'].append(df['tname']), df['qlen'].append(df['tlen'])))
    cluster_df['qlen'] = cluster_df['qname'].map(rls)
    
    cluster_df['mean_qlen'] = cluster_df.groupby('cluster')['qlen'].transform('mean')
    cluster_df['frac_max_score'] = cluster_df['clust_read_score']/(cluster_df['mean_qlen']/10000)**2
    # add kmer bin info
    bin_id = args.paf.split('/')[-2]
    cluster_df['bin_id'] = bin_id
    # keep reads > args.min_score_frac
    cluster_df = cluster_df[cluster_df['frac_max_score'] > args.min_score_frac]
    # keep clusters > args.min_reads
    cluster_df = cluster_df[cluster_df.groupby('cluster')['qname'].transform('count') >= args.min_reads]
    
    cluster_df = cluster_df[['bin_id', 'cluster', 'qname', 'qlen', 'clust_read_score', 'frac_max_score']]
    cluster_df.rename(columns={'qname': 'read_id', 'qlen': 'read_len'}, inplace=True)
    cluster_df['clust_read_score'] = cluster_df['clust_read_score'].round(3)
    cluster_df['frac_max_score'] = cluster_df['frac_max_score'].round(3)
    # to csv
    cluster_df.to_csv('{}.csv'.format(args.prefix), index=False)

if __name__ == '__main__':
    args = parse_args()
    main(args)
