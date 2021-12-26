#import os,sys
import sys
import pandas as pd
import argparse
#import collections
#import matplotlib; matplotlib.use('Agg')
#import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
#from scipy.spatial.distance import pdist
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
    Z   = sch.linkage(squareform(df, checks=False), 'ward', optimal_ordering=True)
    args.max_dist = 3*np.median(Z[:,2])
    clusters = sch.fcluster(Z, args.max_dist, criterion='distance')
    
    df.loc[:,'cluster'] = clusters
    df                  = df.reset_index()
    cluster_dfs         = df.groupby('cluster')
    return cluster_dfs

def main(args):
    if args.max_recursion > 0:
        sys.setrecursionlimit(args.max_recursion)

    bin_id = args.paf.split('/')[-2]
    
    df = pd.DataFrame()              
    cols = ['qname', 'qlen', 'tname', 'tlen', 'alnlen', 'dv']
    dtypes = {0:"category", 1:"int16", 5:"category", 6:"int16", 10:"int16", 15:"category"}
    chunksize = 10 ** 6
    for chunk in pd.read_csv(args.paf, sep='\t', header=None, usecols=[0,1,5,6,10,15], names=cols, dtype=dtypes, chunksize=chunksize):
        if chunk.shape[0] > 0:
            chunk['sim'] = 1 - pd.Series([x.split(':')[-1] for x in chunk['dv']], dtype=np.float32)
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
    cluster_dfs = run_clustering(df_pivot, args, nan_value)

    # prepare clustering info
    clust_info_out = '{}.info.csv'.format(args.prefix)
    # dict to store sequence length
    rls = dict(zip(df['qname'].append(df['tname']), df['qlen'].append(df['tlen'])))
    with open(clust_info_out, 'w') as f:
        f.write('bin_id,cluster,read_id,read_len,clust_read_score,frac_max_score\n')
        for cluster,df_c in cluster_dfs:
            keep_cols = df_c['qname'].append(pd.Series(['qname', 'cluster']))
            df_c = df_c.loc[:,keep_cols] # only keep the columns of reads from the cluster
            rls_c = []
            read_means = []
            for idx,row in df_c.iterrows():
                nonself   = row.replace(1.00, np.nan).dropna()
                dist_idx   = [i for i in nonself.index if (i!='qname' and i!='cluster')]
                read_means.append(nonself.loc[dist_idx].mean())
                rls_c.append(rls[row['qname']])

            df_c.loc[:,'score_mean'] = read_means
            clust_mean_score = df_c.loc[:,'score_mean'].mean()
            clust_mean_rl = np.mean(rls_c)
            max_score = (clust_mean_rl/10000)**2
            max_score_frac = clust_mean_score / max_score
            
            df_c = df_c.set_index('qname')
            if df_c.shape[0]>=args.min_reads and max_score_frac>=args.min_score_frac:
                for r in df_c.index.values:
                    f.write('{},{},{},{},{},{}\n'.format(bin_id, \
                                                         cluster, \
                                                         r, \
                                                         rls[r], \
                                                         round(df_c.loc[r,'score_mean'],3), \
                                                         round(max_score_frac,3)))

if __name__ == '__main__':
    args = parse_args()
    main(args)
