from pipeline.toolkits import utils
import subprocess
from pipeline.__init__ import RUN_THREADS
import pandas as pd

class TRIM():
    """
    Features:
    - Trim the adapter sequence in read.
    - Trim low-quality bases in reads.
    
    Arguments:
    - `fq` Fastq file path, required.
    - `adaptor_3`: Adaptor sequence in 3'end. Required.
    - `min_qual` Minimum base quality. Default `30`.
    - `min_length` Minimum seq length after trimming. Default `50`.
    - `error_rate` Maximum allowed error rate (if 0 <= E < 1), or absolute number of errors for full-length adapter match (if E is an integer >= 1). Error rate = no. of errors divided by length of matching region. Default: `0.1`.
    - `cutadapt_para` Other cutadapt parameters. Default `None`.
    
    Output:
    - `{sample}.clean.fq` Clean fastq file after trimming.
    - `stat.txt` Trim summary.
    """
    def __init__(self, args, step):
        self.step = step
        ## required parameters
        self.outdir = args.outdir
        self.sample = args.sample
        self.fq = args.fq
        self.adaptor_3 = args.adaptor_3
        ## default parameters
        self.n_cores = RUN_THREADS[self.step]
        self.min_qual = args.min_qual
        self.min_length = args.min_length
        self.error_rate = args.error_rate
        ## other cutadapt parameters
        self.cutadapt_para = args.cutadapt_para
        
    @utils.logit
    def process_trim(self):
        utils.check_dir(f'{self.outdir}')
        self.out_fq = f'{self.outdir}/{self.sample}.clean.fq'
        cmd = ['cutadapt', 
               '-a', f'{self.adaptor_3}', 
               '-o', self.out_fq,
               '-j', f'{self.n_cores}',
               '-m', f'{self.min_length}', 
               '-q', f'{self.min_qual}', 
               '-e', f'{self.error_rate}',
               f'{self.fq}']
        cmd = ' '.join(cmd)
        
        if self.cutadapt_para != None:
            cmd += f' {self.cutadapt_para}'
            
        cmd += f' > {self.outdir}/trim.log 2>&1'
        TRIM.process_trim.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)
        
    @utils.logit
    def fastqc(self):
        cmd = f'fastqc ' \
              f'-t 3 ' \
              f'-o {self.outdir} ' \
              f'{self.out_fq}'
        TRIM.fastqc.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)
    
    def gen_summary(self):
        dic = {}
        go = True
        for line in utils.iter_readline(f'{self.outdir}/trim.log'):
            if '=== Adapter 1 ===' in line:
                go = False
            if go:
                if 'Total reads processed:' in line:
                    dic['Total reads processed:'] = line.strip('Total reads processed:').strip(' ').strip('\n')
                if 'Reads with adapters:' in line:
                    dic['Reads with adapters:'] = line.strip('Reads with adapters:').strip(' ').strip('\n')
                if 'Reads that were too short:' in line:
                    dic['Reads that were too short:'] = line.strip('Reads that were too short:').strip(' ').strip('\n')
                if 'Reads written (passing filters)' in line:
                    dic['Reads written (passing filters)'] = line.strip('Reads written (passing filters)').strip(' ').strip('\n')
                if 'Total basepairs processed:' in line:
                    dic['Total basepairs processed:'] = line.strip('Total basepairs processed:').strip(' ').strip('\n')
                if 'Quality-trimmed:' in line:
                    dic['Quality-trimmed:'] = line.strip('Quality-trimmed:').strip(' ').strip('\n')
                if 'Total written (filtered):' in line:
                    dic['Total written (filtered):'] = line.strip('Total written (filtered):').strip(' ').strip('\n')
            else:
                break
        
        trim_summary = pd.DataFrame.from_dict(dic, orient="index")
        trim_summary.to_csv(f'{self.outdir}/{self.sample}_summary.txt', sep='\t', header=False)
        
    @utils.logit 
    def run(self):
        self.process_trim()
        self.fastqc()
        self.gen_summary()
        
        
def trim(args):
    step = 'trim'
    trim_obj = TRIM(args, step)
    trim_obj.run()
        
        
def get_trim_para(parser, optional=False):
    parser.add_argument("--fq", help="fastq data", required=True)
    parser.add_argument("--adaptor_3", 
                            help="3' end adaptor sequence.", 
                            default="GCGGAAGCAGTGGTATCAACGCAGAGTACAACAAGGTAC")
    if optional:
        parser.add_argument("-ml", "--min_length", 
                            help="Minum length after trimming.", default=50)
        parser.add_argument("-mq", "--min_qual", 
                            help="Minum base quality.", default=30) 
        parser.add_argument("-er", "--error_rate", 
                            help="Maximum allowed error rate (if 0 <= E < 1), \
                            or absolute number of errors for full-length adapter match (if E is an integer >= 1). \
                            Error rate = no. of errors divided by length of matching region",
                            default=0.1)
        parser.add_argument("--cutadapt_para", help="Cutadapt paramters.", default=None)   
        parser = utils.common_args(parser)
    return(parser)
    
    