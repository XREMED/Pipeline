import gzip
from collections import defaultdict
from itertools import combinations, product

import mappy as mpp
import pandas as pd
from pipeline.toolkits.utils import check_dir, check_file, common_args, logit


def findall_mismatch(seq, n_mismatch=1, bases="ACGTN"):
    """
    choose locations where there's going to be a mismatch using combinations
    and then construct all satisfying lists using product
    Return:
    all mismatch <= n_mismatch set. 
    >>> answer = set(["TCG", "AAG", "ACC", "ATG", "ACT", "ACN", "GCG", "ANG", "ACA", "ACG", "CCG", "AGG", "NCG"])
    >>> seq_set = seq.findall_mismatch("ACG")
    >>> seq_set == answer
    True
    """
    seq_set = set()
    seq_len = len(seq)
    if n_mismatch > seq_len:
        n_mismatch = seq_len
    for locs in combinations(range(seq_len), n_mismatch):
        seq_locs = [[base] for base in seq]
        for loc in locs:
            seq_locs[loc] = list(bases)
        for poss in product(*seq_locs):
            seq_set.add("".join(poss))
    return list(seq_set)


def generate_seq_dict(seq_list):
    """_summary_

    Args:
        seq_list (str): text file, containing the seqs list. 
    
    Return:
        seq_dict: {mismatch_seq: raw_seq}      
    """
    
    ### check file path
    check_file([seq_list])
    seq_dict = {}
    with open(seq_list) as fh:
        lines = fh.readlines()
        for line in lines:
            bc = line.strip("\n").split("\t")[1]
            mismatch_seqs = findall_mismatch(bc, n_mismatch=1)
            for i in mismatch_seqs:
                seq_dict[i] = bc
    return seq_dict
    
    
def correct_seq(mis_seq, seq_dict):
    return seq_dict[mis_seq]


class BARCODE:
    """
    Features: 
    - Extract the Barcode and UMI information in R1, and use it as the header of R2 read.
        
    Arguments:
    - fq1: R1 read path, required.
    - fq2: R2 read path, required.
    - barcode_list: 
    - barcode_range: Barcode range in the R1 read. Default: 1,10
    - umi_range: UMI range in the R1 read. Default:11,20
    - qual_filter: Whether to filter barcodes by base quality.
    - qual_min: Minimum base quality in a barcode sequence, only useful if "qual_filter" is present. Default:30
    - gzip: Whether to output fastq files in compressed format.

    Outputs:
    - {sample}.fq(.gz): R2 data with modified read header.
    - {samle}_readcount.tsv: Read count statistics for each barcode.
    """
    def __init__(self, args, step):
        self.step = step
        ### required para
        self.fq1 = args.fq1
        self.fq2 = args.fq2
        self.barcode_list = args.barcode_list
        self.sample = args.sample
        ### para with a default value 
        self.outdir = args.outdir
        self.barcode_range = args.barcode_range
        self.umi_range = args.umi_range
        self.qual_filter = args.qual_filter
        if self.qual_filter:
            self.qual_min = args.qual_min
        else:
            self.qual_min = None
        self.gzip = args.gzip
        
    @logit    
    def run(self):
        ### INPUT
        check_file([self.fq1, self.fq2])
        
        f1, f2 = mpp.fastx_read(self.fq1), mpp.fastx_read(self.fq2)
        barcode_dict = generate_seq_dict(self.barcode_list)
        
        ### barcode statistics
        barcode_range = self.barcode_range.split(",")
        umi_range = self.umi_range.split(",")
        
        total_reads = 0
        valid_barcode_reads = 0
        incorrect_barcode_reads = 0
        low_qual_barcode_reads = 0
        
        count_dict = defaultdict(lambda: defaultdict(int))

        ### make outdir
        check_dir([f'{self.outdir}/{self.step}'])
            
        if self.gzip:
            out_fq = gzip.open(f'{self.outdir}/{self.step}/{self.sample}.fq.gz', "wt")
        else:
            out_fq = open(f'{self.outdir}/{self.step}/{self.sample}.fq', "wt")
            
        for entry1, entry2 in zip(f1, f2):
            total_reads += 1
            if total_reads % 5000000 == 0:
                BARCODE.run.logger.info(f'Processed {total_reads} reads.')
            f1_seq = entry1[1]
            f1_qual = entry1[2]
            barcode = f1_seq[int(barcode_range[0])-1:int(barcode_range[1])]
            umi = f1_seq[int(umi_range[0])-1:int(umi_range[1])]
            
            if barcode in barcode_dict:
                if self.qual_min is not None:
                    ### check barcode quality:
                    bc_qual = f1_qual[int(barcode_range[0])-1:int(barcode_range[1])]
                    bc_qual = [ord(i)-33 for i in bc_qual]
                    if min(bc_qual) < int(self.qual_min):
                        low_qual_barcode_reads += 1
                        continue   
                    else:
                        count_dict[barcode_dict[barcode]][umi] += 1
                        valid_barcode_reads += 1
                else:
                    count_dict[barcode_dict[barcode]][umi] += 1
                    valid_barcode_reads += 1
                new_seq = f'@{barcode_dict[barcode]}_{umi}_{total_reads}\n{entry2[1]}\n+\n{entry2[2]}\n'  
                out_fq.write(f'{new_seq}')     
            else:
                incorrect_barcode_reads += 1
                continue  
        
        out_fq.close()
                
        df = pd.DataFrame.from_dict(count_dict, orient='index')
        df.to_csv(f'{self.outdir}/{self.step}/{self.sample}_readcount.tsv', 
                  sep='\t', header=False)
        
        ### sum barcode:
        barcode_summary = {
            "Total reads": total_reads,
            "Reads with valid barcode": valid_barcode_reads,
            "Reads with incorrect barcode": incorrect_barcode_reads,
            "Reads with low quality barcode": low_qual_barcode_reads
        }
        barcode_summary = pd.DataFrame.from_dict(barcode_summary, orient="index")
        barcode_summary.to_csv(f'{self.outdir}/{self.step}/summary.txt', sep='\t', header=False)
        
def barcode(args):
    step = "barcode"
    barcode_obj = BARCODE(args, step)
    barcode_obj.run()
    
    
def get_barcode_para(parser):
    parser.add_argument("--fq1", help="R1 read.", required=True)
    parser.add_argument("--fq2", help="R2 read.", required=True)
    parser.add_argument("--barcode_list", help="Barcode list file.", required=True)
    parser.add_argument("--qual_filter", help="Filter barcode by base quality.", action='store_true')
    parser.add_argument("--gzip", help="Output gzip fastq file.", action="store_true")
    parser.add_argument("--barcode_range", help="Barcode range in Read 1.", default="1,10")
    parser.add_argument("--umi_range", help="UMI range in Read 1.", default="11,20")
    parser.add_argument("--qual_min", help="Min barcode base quality, only used when qual_filter is specified", default=30)

    
    return common_args(parser)
