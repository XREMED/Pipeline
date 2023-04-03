import importlib
import logging
import os
import re
import sys
import time
from collections import Counter, defaultdict
from datetime import timedelta
from functools import wraps


def logit(func):
    '''
    logging start and done.
    '''
    logging.basicConfig(level=logging.INFO, 
                        stream=sys.stdout,
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', 
                        datefmt='%Y/%m/%d %H:%M:%S')
    module = func.__module__
    name = func.__name__
    logger_name = f"{module}.{name}"
    logger = logging.getLogger(logger_name)

    @wraps(func)
    def wrapper(*args, **kwargs):
        if args and hasattr(args[0], 'debug') and args[0].debug:
            logger.setLevel(10) # debug

        logger.info('start...')
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        used = timedelta(seconds=end - start)
        logger.info('done. time used: %s', used)
        return result

    wrapper.logger = logger
    return wrapper 


@logit
def parse_sample(sample_file: str):
    """
    Args:
        sample_file (string): There are four columns in total, which are sample name, fastq file path, experimental processing, and remark.
        
    Return:
        sample dict (dict)
    """
    
    def split_line(line):
        return line.rstrip('\n').split('\t')
    
    dic = defaultdict(dict)
    
    # check sample file path and format
    try:
        sample_lines = open(sample_file) 
    except IOError:
        print(f'ERROR: No such file: "{sample_file}"!')
    else:
        if not sample_file.endswith('txt'):
            suffix_c = sample_file.split('.')[-1]
            raise Exception(f'ERROR: Invalid file format:.{suffix_c}! .txt file is required for sample file!')
        else:
            line = sample_lines.readline()
            s = split_line(line)
            if len(s) != 3 and len(s) != 4:
                raise Exception(f'ERROR: Invaild separation! Sample file should be separated by tab.')
            while line:
                line = sample_lines.readline()
                s = split_line(line)
                sample = s[0]
                fq1 = s[1].split(',')[0]
                fq2 = s[1].split(',')[1]
                if not os.path.exists(fq1):
                    raise IOError(f'ERROR: No such fastq file: "{fq1}"')
                elif not os.path.exists(fq2):
                    raise IOError(f'ERROR: No such fastq file: "{fq2}"')
                treat = s[2]
                if len(s) == 4:
                    remark = s[3]
                else:
                    remark = None
                dic[sample]['fq1'] = fq1
                dic[sample]['fq2'] = fq2
                dic[sample]['treat'] = treat
                dic[sample]['remark'] = remark
                
                return dic
            

def find_assay_init(assay):
    init_module = importlib.import_module(f"pipeline.{assay}.__init__")
    return init_module


def find_step_module(assay, step):
    init_module = find_assay_init(assay)
    try:
        step_module = importlib.import_module(f"pipeline.{assay}.{step}")
    except ModuleNotFoundError:
        try:
            step_module = importlib.import_module(f"pipeline.toolkits.{step}")
        except ModuleNotFoundError:
            module_path = init_module.IMPORT_DICT[step]
            step_module = importlib.import_module(f"{module_path}.{step}")

    return step_module

         
def common_args(parser):
    ### 这是一些通用的参数。
    parser.add_argument("-s", "--sample", help="Sample name")
    parser.add_argument("-o", "--outdir", help="Output directory.", default='./')
    parser.add_argument("-t", "--thread", help="Number of threads for each step.", default=5)
    return parser


def check_file(file_list):
    s = 'No such file(s): '
    if isinstance(file_list, list):
        for file in file_list:
            if not os.path.exists(file):
                s += f' {file}'
            else:
                continue
    else:
        if not os.path.exists(file_list):
            s += file_list

    if not s == 'No such file(s): ':
        raise FileNotFoundError(s)
    else:
        pass

    
def check_dir(dirs):
    if isinstance(dirs, list):
        for dir in dirs:
            if not os.path.exists(dir):
                os.system(f'mkdir -p {dir}')
            else:
                continue
    else:
        if not os.path.exists(dirs):
            os.system(f'mkdir -p {dirs}')
            
def iter_readline(file):
    with open(file, 'r') as f:
        while True:
            line = f.readline()
            if line:
                yield line
            else:
                return
            
            
def get_id_name_dict(gtf_file):
    """
    get gene_id:gene_name from gtf file
        - one gene_name with multiple gene_id: "_{count}" will be added to gene_name.
        - one gene_id with multiple gene_name: error.
        - duplicated (gene_name, gene_id): ignore duplicated records and print a warning.
    Returns:
        {gene_id: gene_name} dict
    """

    gene_id_pattern = re.compile(r'gene_id "(\S+)";')
    gene_name_pattern = re.compile(r'gene_name "(\S+)"')
    id_name = {}
    c = Counter()
    with open(gtf_file) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('#'):
                continue
            tabs = line.split('\t')
            gtf_type, attributes = tabs[2], tabs[-1]
            if gtf_type == 'gene':
                gene_id = gene_id_pattern.findall(attributes)[-1]
                gene_names = gene_name_pattern.findall(attributes)
                if not gene_names:
                    gene_name = gene_id 
                else:
                    gene_name = gene_names[-1]
                c[gene_name] += 1
                if c[gene_name] > 1:
                    if gene_id in id_name:
                        assert id_name[gene_id] == gene_name, (
                                'one gene_id with multiple gene_name '
                                f'gene_id: {gene_id}, '
                                f'gene_name this line: {gene_name}'
                                f'gene_name previous line: {id_name[gene_id]}'
                            )
                        get_id_name_dict.logger.warning(
                                'duplicated (gene_id, gene_name)'
                                f'gene_id: {gene_id}, '
                                f'gene_name {gene_name}'
                            )
                        c[gene_name] -= 1
                    else:
                        gene_name = f'{gene_name}_{c[gene_name]}'
                id_name[gene_id] = gene_name
    return id_name


def hamming_distance(string1, string2):
    distance = 0
    length = len(string1)
    length2 = len(string2)
    if (length != length2):
        raise Exception(f"string1({length}) and string2({length2}) do not have same length")
    for i in range(length):
        if string1[i] != string2[i]:
            distance += 1
    return distance