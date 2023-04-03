from pipeline.toolkits import utils
import os
import subprocess

TOOLS_DIR = os.path.dirname(__file__)
ANA_TOOLS = f'{TOOLS_DIR}/analysis.R'

class ANALYSIS:
    def __init__(self, args, step):
        self.step = step
        self.outdir = args.outdir
        self.count_matrix = args.count_matrix
        self.group_info = args.group_info
        self.project = args.project
        self.group_by = args.group_by
        self.method = args.method
        self.min_cells = args.min_cells
        self.min_features = args.min_features
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
            f'--method {self.method} ' \
            f'--min_cells {self.min_cells} ' \
            f'--min_features {self.min_features} ' \
            f'--dims {self.dims} ' \
            f'--k_filter {self.k_filter}'
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
                        default='plate')
    parser.add_argument('--method', 
                        help='Normalize method for integration.', 
                        default='integrate',
                        choices=['integrate', 'SCTransform'])
    parser.add_argument('--min_cells', 
                        help='Include features detected in at least this many cells. \
                    Will subset the counts matrix as well. To reintroduce excluded features, create a new object with a lower cutoff.', 
                        default=4)
    parser.add_argument('--min_features', 
                        help='Include cells where at least this many features are detected.', 
                        default=200)
    parser.add_argument('--dims', 
                        help='PCA nums.', 
                        default=6)
    parser.add_argument('--k_filter', 
                        help='Mininum cell nums of sample.', 
                        default=80)
    if optional:
        parser = utils.common_args(parser)
    return parser
