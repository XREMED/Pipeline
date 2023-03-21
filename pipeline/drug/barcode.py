import gzip
from itertools import combinations, product

import mappy as mpp
import pandas as pd
from pipeline.toolkits.utils import check_dir, check_file, common_args, logit


def findall_mismatch(seq, n_mismatch, bases="ACGTN"):
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


def generate_seq_dict(seq_list, n_mismatch):
    """_summary_

    Args:
        seq_list (str): text file, containing the seqs list. 
        n_mismatch (int): Maxium allowed mismatch num.
    
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
            mismatch_seqs = findall_mismatch(bc, n_mismatch)
            for i in mismatch_seqs:
                seq_dict[i] = bc
    return seq_dict
    
    
def correct_seq(mis_seq, seq_dict):
    return seq_dict[mis_seq]


# def rescue(seq, qual, barcode_range, umi_range):
#     barcode = seq[barcode_range[0]-1:barcode_range[1]]
#     umi = seq[umi_range[0]-1:umi_range[1]]
#     # check last base in umi
#     if umi[-1]=='T':
#         # check the 5 bases after umi
#         # if the 5 bases are all 'T'
#         # insert 'N' to the last position of barcode, and shift umi
#         if seq[umi_range[1]:umi_range[1]+5]=='TTTTT':
#             new_seq = barcode[:-1] + 'N' + seq[barcode_range[1]-1:]
#         else:
#             new_seq = seq
#     else:
#         new_seq = seq


class BARCODE:
    """
    Features: 
    - Extract the Barcode and UMI information in R1, and use it as the header of R2 read.
    - Filter barcode: Only one base mismatch is allowed at most, and the mismatched base must be a low-quality base.
        
    Arguments:
    - `fq1` R1 read path, required.
    - `fq2` R2 read path, required.
    - `barcode_list` Barcode file path. E.g. "barcode_name\tbarcode_seq". Required.
    - `barcode_range` Barcode range in the R1 read. Default: `1,10`.
    - `umi_range` UMI range in the R1 read. Default: `11,20`.
    - `min_qual` Minimum base quality in a barcode sequence. Default: `20`.
    - `gzip` Output fastq files in compressed format.

    Outputs:
    - `{sample}.fq(.gz)` R2 data with modified read header.
    - `stat.txt` Barcode summary.
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
        self.n_mismatch = 1
        self.min_qual = int(args.min_qual)

        self.gzip = args.gzip
        
    @logit    
    def run(self):
        ### INPUT
        check_file([self.fq1, self.fq2])
        
        f1, f2 = mpp.fastx_read(self.fq1), mpp.fastx_read(self.fq2)
        barcode_dict = generate_seq_dict(self.barcode_list, self.n_mismatch)
        
        ### barcode statistics
        barcode_range = self.barcode_range.split(",")
        umi_range = self.umi_range.split(",")
        
        total_reads = 0  # total reads
        valid_barcode_reads = 0 # read with a valid barcode
        corrected_barcode_reads = 0 # read with a corrected barcode
        incorrect_barcode_reads = 0 # read with a incorrect barcode: 1. mismatch>1 or 
                                    # 2. mismatch=1 but mismatch base quality > min_qual

        ### make outdir
        check_dir([f'{self.outdir}'])
            
        if self.gzip:
            out_fq = gzip.open(f'{self.outdir}/{self.sample}.fq.gz', "wt")
        else:
            out_fq = open(f'{self.outdir}/{self.sample}.fq', "wt")
            
        for entry1, entry2 in zip(f1, f2):
            total_reads += 1
            if total_reads % 5000000 == 0:
                BARCODE.run.logger.info(f'Processed {total_reads} reads.')
            f1_seq = entry1[1]
            f1_qual = entry1[2]
            barcode = f1_seq[int(barcode_range[0])-1:int(barcode_range[1])]
            umi = f1_seq[int(umi_range[0])-1:int(umi_range[1])]
            umi_qual = f1_qual[int(umi_range[0])-1:int(umi_range[1])]
            
            ### At most one base mismatch is allowed, and the base must be a low-quality base.
            if barcode in barcode_dict:
                ### check barcode quality:
                bc_qual = f1_qual[int(barcode_range[0])-1:int(barcode_range[1])]
                bc_qual = [ord(i)-33 for i in bc_qual]
                diff_idx = [i for i in range(len(barcode)) if barcode[i]!=barcode_dict[barcode][i]]
                if diff_idx!=[]:
                    if diff_idx[0] < self.min_qual:
                        corrected_barcode_reads += 1
                        bc = barcode_dict[barcode]
                    else:
                        incorrect_barcode_reads += 1
                        continue
                else:
                    bc = barcode
                    valid_barcode_reads += 1
            else:
                incorrect_barcode_reads += 1
                continue
            
            new_head = f'@{bc}-{self.sample}_{umi}_{umi_qual}_{total_reads}'
            new_seq = f'{new_head}\n{entry2[1]}\n+\n{entry2[2]}\n'  
            out_fq.write(f'{new_seq}')      
        
        out_fq.close()
                
        
        ### sum barcode:
        barcode_summary = {
            "Total reads": total_reads,
            "Reads with a valid barcode": f'{valid_barcode_reads} ({round(valid_barcode_reads*100/total_reads, 2)}%)',
            "Reads with a corrected barcode": f'{corrected_barcode_reads} ({round(corrected_barcode_reads*100/total_reads, 2)}%)',
            "Reads with a incorrect barcode": f'{incorrect_barcode_reads} ({round(incorrect_barcode_reads*100/total_reads, 2)}%)',
        }
        barcode_summary = pd.DataFrame.from_dict(barcode_summary, orient="index")
        barcode_summary.to_csv(f'{self.outdir}/{self.sample}_summary.txt', sep='\t', header=False)
        
def barcode(args):
    step = "barcode"
    barcode_obj = BARCODE(args, step)
    barcode_obj.run()
    
    
def get_barcode_para(parser, optional=False):
    
    parser.add_argument("--fq1", help="R1 read.", required=True)
    parser.add_argument("--fq2", help="R2 read.", required=True)
    parser.add_argument("--barcode_list", help="Barcode list file.", 
                        required=True)
    parser.add_argument("--barcode_range", help="Barcode range in Read 1.",
                        default="1,10")
    parser.add_argument("--umi_range", help="UMI range in Read 1.", 
                        default="11,20")
    if optional:
        parser.add_argument("--gzip", help="Output gzip fastq file.", 
                            action="store_true")
        parser.add_argument("--min_qual", help="Min barcode base quality", 
                            default=20)
        parser = common_args(parser)
    
    return parser
