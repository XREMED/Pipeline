import argparse
import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from glob import glob

import pandas as pd

from pipeline.drug.__init__ import __STEPS__
from pipeline.pipeline import ArgFormatter
from pipeline.toolkits import utils


def get_fq(path, end='1'):
    suffixs = ['.gz', '.gzip', '']
    formats = ['.fq', '.fastq']
    links = [f'_{end}_', f'_R{end}_', f'_R{end}_001', f'_R{end}', f'_{end}']
    patterns = [f'{path}*{l}*{f}{s}' for s in suffixs for l in links for f in formats]
    fqs = set([glob(p)[0] for p in patterns if glob(p)])
    fqs = sorted(list(fqs))
    if len(fqs)==0:
        s = '\n'.join(patterns)
        raise IOError(f'No valid file path found in: {s}, \
                      please check your path and read name!')
    return(fqs)

def parse_SampleInfo(sample_info_file_path, ) -> dict:
    df = pd.read_csv(sample_info_file_path, sep='\t')
    df.sample_path = df['path'] + '/' + df['sample']
    samples = df.sample_path.tolist()
    sample_dict = defaultdict(dict)
    for s in samples:
        print(s)
        num = s.split('/')[-1]
        fq1 = get_fq(s, '1')
        s1 = ' '.join(fq1)

        fq2 = get_fq(s, '2')
        s2 = ' '.join(fq2)
        if len(fq1) != len(fq2):
            raise Exception('R1 and R2 file numbers are not equal.')
        if len(fq1) > 1:
            os.system(f'cat {s1} > {s}_R1.fq.gz')
            sample_dict[num]['R1'] = f'{s}_R1.fq.gz'
            os.system(f'cat {s2} > {s}_R2.fq.gz')
            sample_dict[num]['R2'] = f'{s}_R2.fq.gz'
            os.system(f'rm -rf {s1} {s2}')
        else:
            sample_dict[num]['R1'] = s1
            sample_dict[num]['R2'] = s2
    return(sample_dict)
        
class ALLINONE:
    def __init__(self, args):
        self.sample_info = args.sample_info
        self.assay = args.assay
        self.outdir = args.outdir
        # barcode
        self.barcode_list = args.barcode_list
        self.barcode_range = args.barcode_range
        self.umi_range = args.umi_range
        # trim
        self.adaptor_3 = args.adaptor_3
        self.min_qual = args.min_qual
        self.min_length = args.min_length
        self.er = args.error_rate
        # mapping
        self.genomeDir = args.genomeDir
        self.gtf = args.gtf
        self.mapping_para = args.mapping_para
        # featureCounts
        self.featureCounts_para = args.featureCounts_para
        
        # prepare
        # prepare sample_dict, sample outdir
        self.sample_dict = parse_SampleInfo(self.sample_info)
        self.outdir_dict = defaultdict()
        
        for s in __STEPS__[:-1]:
            idx = __STEPS__.index(s) + 1
            self.outdir_dict[s] = f'{self.outdir}/{"%02d" %idx}.{s}' 
             
        
    def get_parser(self, assay, step):
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers(required=True)
        step_module = utils.find_step_module(assay, step)
        func = getattr(step_module, step)
        func_opts = getattr(step_module, f"get_{step}_para")
        parser_step = subparsers.add_parser(step, formatter_class=ArgFormatter)
        func_opts(parser_step, optional=True)
        parser_step.set_defaults(func=func)
        return parser
      
        
    @utils.logit
    def barcode(self, fq1, fq2, sample):
        step = 'barcode'
        parser = self.get_parser(self.assay, step)
        cmd_str = f'{step} --fq1 {fq1} ' \
                  f'--fq2 {fq2} ' \
                  f'--barcode_list {self.barcode_list} ' \
                f'-o {self.outdir_dict[step]} -s {sample} '
        args = parser.parse_args(cmd_str.split())
        ALLINONE.barcode.logger.info(cmd_str)
        args.func(args)
     
    @utils.logit   
    def trim(self, fq, sample):
        step = 'trim'
        parser = self.get_parser(self.assay, step)
        cmd_str = f'{step} --fq {fq} --adaptor_3 {self.adaptor_3} ' \
                f'-o {self.outdir_dict[step]} -s {sample} -ml {self.min_length} ' \
                f'-mq {self.min_qual} -er {self.er} '
        args = parser.parse_args(cmd_str.split())
        ALLINONE.trim.logger.info(cmd_str)
        args.func(args)
      
    @utils.logit  
    def mapping(self, clean_fq, sample):
        step = 'mapping'
        parser = self.get_parser(self.assay, step)
        cmd_str = f'{step} --clean_fq {clean_fq} ' \
                f'--genomeDir {self.genomeDir} ' \
                f'-o {self.outdir_dict[step]} -s {sample} '
        if self.mapping_para != None:
            cmd_str += self.mapping_para
        args = parser.parse_args(cmd_str.split())
        ALLINONE.trim.logger.info(cmd_str)
        args.func(args)     
    
    @utils.logit   
    def featureCounts(self, input_bam, sample): 
        step = 'featureCounts'
        parser = self.get_parser(self.assay, step)
        cmd_str = f'{step} --gtf {self.gtf} --input_bam {input_bam} ' \
                f'-o {self.outdir_dict[step]} -s {sample} '
        if self.featureCounts_para != None:
            cmd_str += self.featureCounts_para
        args = parser.parse_args(cmd_str.split())
        ALLINONE.trim.logger.info(cmd_str)
        args.func(args)   
        
    def merge_bam(self):
        bams = glob(self.outdir_dict['featureCounts']+'/*_name_sorted.bam')
        s = ' '.join(bams)
        outs = self.outdir_dict['featureCounts']
        cmd = f'samtools merge -n {outs}/all_name_sorted.bam {s}'
        os.system(cmd)
        
    @utils.logit
    def count(self):
        self.merge_bam()
        step = 'count'
        parser = self.get_parser(self.assay, step)
        bam = self.outdir_dict['featureCounts'] + '/all_name_sorted.bam'
        outs = self.outdir_dict['count']
        cmd = f'{step} --bam {bam} ' \
            f'--gtf {self.gtf} -o {outs} -s all'
        args = parser.parse_args(cmd.split())
        ALLINONE.count.logger.info(cmd)
        args.func(args)    
    
    def one(self, fq1, fq2, sample):
        self.barcode(fq1, fq2, sample)
        self.trim(self.outdir_dict['barcode']+f'/{sample}.fq', sample)
        self.mapping(self.outdir_dict['trim']+f'/{sample}.clean.fq', sample)
        self.featureCounts(self.outdir_dict['mapping']+f'/{sample}_Aligned.sortedByCoord.out.bam', 
                           sample)
        
    def run(self):
        fq1s, fq2s, samples = [], [], []
        for i in self.sample_dict:
            fq1s.append(self.sample_dict[i]['R1'])
            fq2s.append(self.sample_dict[i]['R2'])
            samples.append(i)

        process_pools = []
        with ProcessPoolExecutor(3) as pool:
            for i in pool.map(self.one, fq1s, fq2s, samples):
                process_pools.append(i)
        self.count()
                
def allinone(args):
    allinone_obj = ALLINONE(args)
    allinone_obj.run()
        

def get_allinone_para(parser, optional=True):
    parser.add_argument('--sample_info', 
                        help='Required. Sample info file, `.txt` fromat.', 
                        required=True)
    parser.add_argument('--outdir', 
                        help='Default `./`. Output directory.',
                        default='./')
    parser.add_argument('--assay', help='Required. Assay name.', 
                        required=True)
    # barcode
    parser.add_argument('--barcode_list', 
                        help='Required. TXT file, recording barcode info.', 
                        required=True)
    parser.add_argument('--barcode_range', 
                        help='Required. Barcode range in the read 1.', 
                        default='1,10')
    parser.add_argument('--umi_range', 
                        help='Required. UMI range in the read 1.', 
                        default='11,20')
    # trim
    parser.add_argument('--adaptor_3', 
                        help='Required. Adaptor sequence to trim.', 
                        default='GCGGAAGCAGTGGTATCAACGCAGAGTACAACAAGGTAC')
    parser.add_argument('--min_length', 
                        help='Required. Mininum length after trimming.', 
                        default=50)
    parser.add_argument('--min_qual', 
                        help='Required. Mininum quality for trimming.', 
                        default=30)
    parser.add_argument('--error_rate', 
                        help='Required. Maximum allowed error rate (if 0 <= E < 1), \
                            or absolute number of errors for full-length adapter match (if E is an integer >= 1). \
                            Error rate = no. of errors divided by length of matching region (default: 0.1).', 
                        default=0.1)
    # mapping
    parser.add_argument('--gtf', 
                        help='Required. GTF path.', 
                        required=True)
    parser.add_argument('--genomeDir', 
                        help='Required. STAR genome index directory.', 
                        required=True)

    parser.add_argument('--mapping_para', 
                        help='Mapping parameters.', 
                        default=None)
    parser.add_argument('--featureCounts_para', 
                        help='FeatureCounts parameters.', 
                        default=None)
    
    return parser
    