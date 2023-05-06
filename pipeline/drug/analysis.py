from pipeline.toolkits import utils
from pipeline.__init__ import ROOT_DIR
import subprocess

ANA_TOOLS = f'{ROOT_DIR}/toolkits/analysis.R'

class ANALYSIS:
    def __init__(self, args, step):
        self.step = step
        self.outdir = args.outdir
        self.count_matrix = args.count_matrix
        self.group_info = args.group_info
        self.project = args.project
        self.group_by = args.group_by
        self.normalization_method = args.normalization_method
        self.dims = args.dims
        self.k_filter = args.k_filter
        
    @utils.logit
    def run(self):
        cmd = f'Rscript {ANA_TOOLS} ' \
            f'--outdir {self.outdir} ' \
            f'--count_matrix {self.count_matrix} ' \
            f'--group_info {self.group_info} ' \
            f'--project {self.project} ' \
            f'--group_by {self.group_by} ' \
            f'--normalization.method {self.normalization_method} ' \
            f'--k_filter {self.k_filter} '
        if self.dims != None:
            cmd+=f'--dims {self.dims}'
        ANALYSIS.run.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)
        

def analysis(args):
    step = 'analysis'
    obj = ANALYSIS(args, step)
    obj.run()
    
    
def get_analysis_para(parser, optional=False):
    parser.add_argument('--count_matrix', 
                        help='Required. Count matrix file, `.txt` fromat.', 
                        required=True)
    parser.add_argument('--group_info', 
                        help='Required. Group info file, `.txt` fromat.', 
                        required=True)
    parser.add_argument('--project', 
                        help='Required. Project name.', 
                        required=True)
    parser.add_argument('--group_by', 
                        help='Required. Group by.', 
                        default='tag')
    parser.add_argument('--normalization_method', 
                        help='Normalize method for integration.', 
                        default='LogNormalize',
                        choices=['LogNormalize', 'SCT'])
    parser.add_argument('--dims', 
                        help='PCA nums.', 
                        default=None)
    parser.add_argument('--k_filter', 
                        help='Mininum cell nums of sample.', 
                        default=30)
    if optional:
        parser = utils.common_args(parser)
    return parser
