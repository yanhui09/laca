import sys
import pandas as pd
import argparse
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('paf', help='Alignment file produced by all-vs-all minimap2', type=str)
    parser.add_argument('-p', '--prefix', help='Output file prefix [output]', type=str, default='output')
    parser.add_argument('-s', '--min_score_frac', help='Minimum fraction of maximum score for a cluster [0.8]', type=float, default=0.8)
    parser.add_argument('-n', '--min_reads', help='Minimum reads to create cluster [5]', type=int, default=5)
    parser.add_argument('-R', '--max_recursion', help='Override the max recursion depth in python', type=int, default=-1)
    args = parser.parse_args()
    return args

def process_chunk(chunk):
    chunk['sim'] = 1 - chunk['dv'].str.split(':').str[-1].astype(float)
    chunk['score'] = (chunk['alnlen'] / 10000) ** 2 * chunk['sim']
    return chunk.drop(columns=['dv', 'alnlen'])

def cluster_data(df):
    Z = sch.linkage(squareform(df, checks=False), 'ward')
    max_dist = 3 * np.median(Z[:, 2])
    clusters = sch.fcluster(Z, max_dist, criterion='distance')
    df['cluster'] = clusters
    # exclude self comparison
    df.replace(0, np.nan, inplace=True)
    # mean score to others within group and across groups
    scores = df.groupby(['cluster']).mean()
    mean_scores = scores.values[[x - 1 for x in clusters], range(scores.shape[1])]
    df_out = pd.DataFrame(data={'qname': df.index, 'clust_read_score': mean_scores, 'cluster': clusters})
    df_out.reset_index(drop=True, inplace=True)
    return df_out

def main(args):
    if args.max_recursion > 0:
        sys.setrecursionlimit(args.max_recursion)

    cols = ['qname', 'qlen', 'tname', 'tlen', 'alnlen', 'dv']
    dtypes = {0:"category", 1:"int16", 5:"category", 6:"int16", 10:"int16", 15:"category"}
    chunksize = 10 ** 6
    df_list = []

    for chunk in pd.read_csv(args.paf, sep='\t', header=None, usecols=[0, 1, 5, 6, 10, 15], names=cols, dtype=dtypes, chunksize=chunksize):
        if chunk.shape[0] > 0:
            df_list.append(process_chunk(chunk))

    df = pd.concat(df_list)
    df.drop_duplicates(['qname', 'tname'], inplace=True)

    # create similarity matrix out of pairwise alignment scores
    idx = sorted(set(df['qname']).union(df['tname']))
    df_pivot = (df.pivot(index='qname', columns='tname', values='score')
                .reindex(index=idx, columns=idx)
                .fillna(0, downcast='infer')
                .pipe(lambda x: x+x.values.T))
    # hierarchical clustering
    df_out = cluster_data(df_pivot)
    rls = dict(zip(df['qname'].append(df['tname']), df['qlen'].append(df['tlen'])))
    df_out['qlen'] = df_out['qname'].map(rls)
    df_out['mean_qlen'] = df_out.groupby('cluster')['qlen'].transform('mean')
    df_out['frac_max_score'] = df_out['clust_read_score'] / ((df_out['mean_qlen'] / 10000) ** 2)
    # bin id, basename
    bin_id = args.paf.split('/')[-1].split('.')[0]
    df_out['bin_id'] = bin_id
    # filter
    df_out = df_out[df_out['frac_max_score'] > args.min_score_frac]
    df_out = df_out[df_out.groupby('cluster')['qname'].transform('count') >= args.min_reads]

    df_out = df_out[['bin_id', 'cluster', 'qname', 'qlen', 'clust_read_score', 'frac_max_score']]
    df_out.rename(columns={'qname': 'read_id', 'qlen': 'read_len'}, inplace=True)
    df_out.to_csv(f'{args.prefix}.csv', index=False)

if __name__ == '__main__':
    args = parse_args()
    main(args)
