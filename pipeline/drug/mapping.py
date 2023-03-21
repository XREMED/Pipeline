import re
import subprocess

import pandas as pd

from pipeline.__init__ import RUN_THREADS
from pipeline.toolkits import utils


class MAPPING():
    """
    Features:
    - Mapping reads to genome and sort bam file.
    
    Arguments:
    - `fq`: Clean fastq file after trimming. Required.
    - `genomeDir` Genome index directory. Required.
    - `out_unmapped` Ouput unmapped reads. 
    - `outFilterMatchNmin` Alignment will be output only if the number of matched bases is higher than or equal to this value. Default `0`.
    - `outFilterMultimapNmax` Maximum number of loci the read is allowed to map to. Alignments (all of them) will be output only if the read maps to no more loci than this value. Otherwise no alignments will be output, and the read will be counted as "mapped to too many loci" in the Log.final.out. Default `1`.
    - `STAR_param` Other star parameters. Default `None`.
    
    Outputs:
    - `{sample}_Aligned.out.bam` Unsorted bam file.
    - `{sample}_Aligned.sortedByCoord.out.bam` Sorted bam file.
    - `{sample}_Log.final.out` STAR summary file.
    - `{sample}_Log.out` STAR log.
    - `{sample}_Log.progress.out` STAR log.
    - `{sample}_SJ.out.tab`
    - `summary.txt` Mapping summary file.
    """
    def __init__(self, step, args):
        self.step = step
        self.outdir = args.outdir
        self.sample = args.sample
        self.fq = args.clean_fq
        self.genomeDir = args.genomeDir
        self.out_unmapped = args.out_unmapped
        self.outFilterMatchNmin = int(args.outFilterMatchNmin)
        self.outFilterMultimapNmax = int(args.outFilterMultimapNmax)
        self.STAR_param = args.STAR_param
        self.thread = args.thread

        # parse
        self.stat_prefix = 'Reads'

        # out
        self.outPrefix = f'{self.outdir}/{self.sample}_'

        self.STAR_map_log = f'{self.outPrefix}Log.final.out'
        self.unsort_STAR_bam = f'{self.outPrefix}Aligned.out.bam'
        self.STAR_bam = f'{self.outPrefix}Aligned.sortedByCoord.out.bam'
    
    @utils.logit
    def STAR(self):
        cmd = [
            'STAR',
            '--runThreadN', str(RUN_THREADS[self.step]),
            '--genomeDir', self.genomeDir,
            '--readFilesIn', self.fq,
            '--outFilterMultimapNmax', str(self.outFilterMultimapNmax),
            '--outFileNamePrefix', self.outPrefix,
            '--outSAMtype', 'BAM', 'Unsorted', # controls sort by Coordinate or not
            '--outFilterMatchNmin', str(self.outFilterMatchNmin)
        ]
        if self.out_unmapped:
            cmd += ['--outReadsUnmapped', 'Fastx']
        if self.fq[-3:] == ".gz":
            cmd += ['--readFilesCommand', 'zcat']
        cmd = ' '.join(cmd)
        if self.STAR_param != None:
            cmd += (" " + self.STAR_param)
        MAPPING.STAR.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    def run_star(self):
        self.STAR()
        self.gen_star_summary()
        self.sort_bam()
        self.index_bam()

    @utils.logit
    def sort_bam(self):
        cmd = (
            f'samtools sort {self.unsort_STAR_bam} '
            f'-o {self.STAR_bam} '
            f'--threads {self.thread} '
        )
        MAPPING.sort_bam.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.logit
    def index_bam(self):
        cmd = f"samtools index {self.STAR_bam}"
        MAPPING.index_bam.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)
    

    def gen_star_summary(self):
        """
        step metrics
        """

        dic = {}
        with open(self.STAR_map_log, 'r') as map_log:
            # number amd percent
            unique_reads_list = []
            multi_reads_list = []
            for line in map_log:
                if line.strip() == '':
                    continue
                if re.search(r'Uniquely mapped reads', line):
                    unique_reads_list.append(line.strip().split()[-1])
                if re.search(r'of reads mapped to too many loci', line):
                    multi_reads_list.append(line.strip().split()[-1])
        unique_reads = int(unique_reads_list[0])
        unique_reads_fraction = unique_reads_list[1]
        multi_reads = int(multi_reads_list[0])
        multi_reads_fraction = multi_reads_list[1]

        dic[f'Uniquely Mapped {self.stat_prefix}:'] = unique_reads
        dic[f'Uniquely Mapped {self.stat_prefix} fraction:'] = unique_reads_fraction
        dic[f'Multi-Mapped {self.stat_prefix}:'] = multi_reads
        dic[f'Multi-Mapped {self.stat_prefix} fraction:'] = multi_reads_fraction

        mappping_summary = pd.DataFrame.from_dict(dic, orient="index")
        mappping_summary.to_csv(f'{self.outdir}/{self.sample}_summary.txt', sep='\t', header=False)


def mapping(args):
    step = 'mapping'
    mapping_obj = MAPPING(step, args)
    mapping_obj.run_star()
    

def get_mapping_para(parser, optional=False):
    parser.add_argument('--clean_fq', help="Required. R2 fastq file.", required=True)
    parser.add_argument(
        '--genomeDir', 
        help='Required. Genome directory.'
    )
    if optional:
        parser.add_argument(
            '--outFilterMatchNmin', 
            help="""Default `0`. Alignment will be output only if the number of matched bases 
    is higher than or equal to this value.""", 
            default=0
        )
        parser.add_argument(
            '--out_unmapped', 
            help='Output unmapped reads', 
            action='store_true'
        )
        parser.add_argument('--STAR_param', help='Other STAR parameters', default=None)
        parser.add_argument(
            '--outFilterMultimapNmax', 
            help='Default `1`. How many places are allowed to match a read at most.', 
            default=1
        )
        # parser.add_argument(
        #     '--starMem', 
        #     help='Default `50`. Maximum memory that STAR can use.', 
        #     default=50
        # )
        
        parser = utils.common_args(parser)
    return parser