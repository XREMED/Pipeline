from collections import defaultdict
from itertools import groupby

import pandas as pd
import plotly.graph_objects as go
import pysam
from plotly.subplots import make_subplots

from pipeline.toolkits import utils


class COUNT():
    """
    Features:
    - Count umi for each gene in each barcode.
    - Filter UMI: 
        1. Cannot contain 'N'.
        2. Cannot be a multimer, such as 'AAAAAAAAAA'.
        3. Cannot have base quality lower than 10.
        
    
    Arguments:
    - `bam` Featurecounts output bam file, containing gene info. Required.
    - `gtf` GTF file path. Required.
    
    Outputs:
    - `{sample}_count.tsv` UMI, read count raw file.
    - `{sample}_matrix.txt` Gene expression matrix.
    """
    def __init__(self, step, args):
        
        # init
        self.step = step
        self.sample = args.sample
        self.outdir = args.outdir
        
        # required parameters
        self.bam = args.bam
        self.gtf = args.gtf
        
        # default parameters
        
        # output files
        utils.check_dir(f'{self.outdir}')
        self.outprefix = f'{self.outdir}/{self.sample}'
        self.count_detail_file = f'{self.outprefix}_count.tsv'
        self.count_matrix = f'{self.outprefix}_matrix.txt'
        self.count_summary = f'{self.outprefix}_metadata.txt'
        
    @staticmethod
    def correct_umi(umi_dict, percent=0.1):
        """
        Correct umi_dict in place.
        Args:
            umi_dict: {umi_seq: umi_count}
            percent: if hamming_distance(low_seq, high_seq) == 1 and
                low_count / high_count < percent, merge low to high.
        Returns:
            n_corrected_umi: int
            n_corrected_read: int
        """
        n_corrected_umi = 0
        n_corrected_read = 0

        # sort by value(UMI count) first, then key(UMI sequence)
        umi_arr = sorted(
            umi_dict.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)
        while True:
            # break when only highest in umi_arr
            if len(umi_arr) == 1:
                break
            umi_low = umi_arr.pop()
            low_seq = umi_low[0]
            low_count = umi_low[1]

            for umi_kv in umi_arr:
                high_seq = umi_kv[0]
                high_count = umi_kv[1]
                if float(low_count / high_count) > percent:
                    break
                if utils.hamming_distance(low_seq, high_seq) == 1:
                    n_low = umi_dict[low_seq]
                    n_corrected_umi += 1
                    n_corrected_read += n_low
                    # merge
                    umi_dict[high_seq] += n_low
                    del (umi_dict[low_seq])
                    break
        return n_corrected_umi, n_corrected_read

    @utils.logit
    def bam2table(self):
        """
        bam to detail table
        must be used on name_sorted bam
        """
        samfile = pysam.AlignmentFile(self.bam, "rb")
        with open(self.count_detail_file, 'wt') as fh1:
            fh1.write('\t'.join(['Barcode', 'geneID', 'UMI', 'count']) + '\n')

            def keyfunc(x): 
                return x.query_name.split('_', 1)[0]
            for _, g in groupby(samfile, keyfunc):
                gene_umi_dict = defaultdict(lambda: defaultdict(int))
                for seg in g:
                    (barcode, umi, umi_qual) = seg.query_name.split('_')[:3]
                    umi_qual = [ord(i)-33 for i in umi_qual]
                    if not seg.has_tag('XT'):
                        continue
                    gene_id = seg.get_tag('XT')
                    # filter umi
                    # Must not be a homopolymer, e.g. AAAAAAAAAA
                    # Must not contain N
                    # Must not contain bases with base quality < 10
                    if len(set(umi))==1 or 'N' in umi or min(umi_qual)<10:
                        continue
                    gene_umi_dict[gene_id][umi] += 1
                for gene_id in gene_umi_dict:
                    # gene_umi_dict[gene_id] = COUNT.correct_umi(gene_umi_dict[gene_id]) # test1
                    COUNT.correct_umi(gene_umi_dict[gene_id])

                # output
                for gene_id in gene_umi_dict:
                    for umi in gene_umi_dict[gene_id]:
                        fh1.write('%s\t%s\t%s\t%s\n' % (barcode, gene_id, umi,
                                                        gene_umi_dict[gene_id][umi]))
        samfile.close()
            
            
    @staticmethod
    def get_df_sum(df, col='UMI'):
        def num_gt2(x):
            return pd.Series.sum(x[x > 1])

        df_sum = df.groupby('Barcode', as_index=False).agg({
            'count': ['sum', num_gt2],
            'UMI': 'count',
            'geneID': 'nunique'
        })
        df_sum.columns = ['Barcode', 'readcount', 'UMI2', 'UMI', 'geneID']
        df_sum = df_sum.sort_values(col, ascending=False)
        return df_sum
    
    @utils.logit
    def write_matrix(self, df):
        # output count matrix and count summary
        df_UMI = df.groupby(['Barcode', 'geneID'], as_index=False).agg({'UMI': 'count'})
        mtx = df_UMI.pivot(values='UMI', 
                        columns='Barcode',
                        index='geneID',).fillna(0).astype(int)
        mtx.insert(0, 'gene_id', mtx.index)

        mtx.insert(0, 'gene_name', mtx['gene_id'].apply(lambda x: self.id_name[x]))
        
        mtx.to_csv(self.count_matrix, sep='\t', index=False)


    @utils.logit
    def plot_violin(self, df_sum):
        fig = go.Figure()
        fig = make_subplots(rows=1,  
                        cols=3,  
                        subplot_titles=["Read count", 
                                        "UMI count", 
                                        "trace2的标题", 
                                        "trace3的标题"], 
                    )
        for i in ['readcount', 'UMI', 'geneID']:
            fig.add_trace(go.Violin(y=df_sum[i],
                                    x=i,
                                    name=i,
                                    box_visible=True,
                                    meanline_visible=True))
        fig.update_layout()
        

    @utils.logit
    def run(self):
        self.id_name = utils.get_id_name_dict(self.gtf)
        self.bam2table()
        df = pd.read_csv(self.count_detail_file, sep='\t')
        self.write_matrix(df)
        df_sum = self.get_df_sum(df)
        df_sum.to_csv(self.count_summary, sep='\t', index=False)
        
        
        
def count(args):
    step = 'count'
    count_obj = COUNT(step, args)
    count_obj.run()
    
    
def get_count_para(parser, optional=False):
    parser.add_argument("--bam", help="Sorted featureCounts output bamfile.",
                        required=True)
    parser.add_argument("--gtf", help="GTF file path.",
                        required=True)
    if optional:
        parser = utils.common_args(parser)
    return(parser)


