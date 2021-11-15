import requests
import json

from math import sqrt

class color:
   PURPLE = '\033[1;35;48m'
   CYAN = '\033[1;36;48m'
   BOLD = '\033[1;37;48m'
   BLUE = '\033[1;34;48m'
   GREEN = '\033[1;32;48m'
   YELLOW = '\033[1;33;48m'
   RED = '\033[1;31;48m'
   BLACK = '\033[1;30;48m'
   UNDERLINE = '\033[4;37;48m'
   END = '\033[1;37;0m'


def download_seq(chr: str, start: int, end: int, genome: str = "hg19", out: str = None):
    r = requests.get(f'https://api.genome.ucsc.edu/getData/sequence?genome={genome};chrom={chr};start={start};end={end};')
    if r.status_code == 200:
        data = r.json()['dna'].upper()
        if out:
            with open(f'{out}.fasta', 'w') as f:
                f.write(f'>{out}\n')
                f.write(data)
        return data
    
    print(f"{color.RED}[ERR]{color.END} {r.status_code}: {r.json()['error']}")



#download_seq('chr15', 25068795, 25223870, out='snrpn')

def compare_tools(gene: str):
    
    qmrlfs = open(f'qmrlfs_data/{gene}.out.bed', 'r').read().split('\n')[:-1]
    rloopt = open(f'rlooptracker_data/{gene}.bedgraph', 'r').read().split('\n')[2:-1]
    
    rloopt_coords = [(int(r[1]), int(r[2])) for r in [r.split(' ') for r in rloopt]]
    qmrlfs_coords = [(int(r[1]), int(r[2])) for r in [r.split('\t') for r in qmrlfs]]

    tp = len(set(rloopt_coords)) if ((set(rloopt_coords) & set(qmrlfs_coords)) == len(set(rloopt_coords))) else len(set(rloopt_coords) & set(qmrlfs_coords))
    tn = 1 if (len(rloopt_coords) == 0 and len(qmrlfs_coords) == 0) else 0
    fp = len(set(qmrlfs_coords) - set(rloopt_coords))
    fn = len(set(rloopt_coords) - set(qmrlfs_coords))

    if fp:
        print(f"False positives [{gene}]: {set(qmrlfs_coords) - set(rloopt_coords)}")
    if fn:
        print(f"False negatives [{gene}]: {set(rloopt_coords) - set(qmrlfs_coords)}")

    return {
        "TP": tp,
        "TN": tn,
        "FP": fp,
        "FN": fn
    }


_TRUE_POSITIVE = 0
_TRUE_NEGATIVE = 0
_FALSE_POSITIVE = 0
_FALSE_NEGATIVE = 0
genes = [
    'apoe',
    'bcl6',
    'fhit',
    'fmr1',
    'hk2',
    'jtb',
    'myc',
    'snrpn',
    'igha1',
    'rhoh',
    'pten',
    'cirh1a',
    'pkd2l1',
    'foxo4',
    'airn',
    'c9orf72'
]

for gene in genes:
    separate_stats = compare_tools(gene)
    
    _TRUE_POSITIVE += separate_stats["TP"]
    _TRUE_NEGATIVE += separate_stats["TN"]
    _FALSE_POSITIVE += separate_stats["FP"]
    _FALSE_NEGATIVE += separate_stats["FN"]

_ACCURACY = (_TRUE_POSITIVE + _TRUE_NEGATIVE) / (_TRUE_POSITIVE + _TRUE_NEGATIVE + _FALSE_POSITIVE + _FALSE_NEGATIVE)
_SENSITIVITY = _TRUE_POSITIVE / (_TRUE_POSITIVE + _FALSE_NEGATIVE)
_SPECIFICITY = _TRUE_NEGATIVE / (_TRUE_NEGATIVE + _FALSE_POSITIVE)
_PRECISION = _TRUE_POSITIVE / (_TRUE_POSITIVE + _FALSE_POSITIVE)
_MCC = (_TRUE_POSITIVE * _TRUE_NEGATIVE - _FALSE_POSITIVE * _FALSE_NEGATIVE) / sqrt((_TRUE_POSITIVE + _FALSE_POSITIVE) * (_TRUE_POSITIVE + _FALSE_NEGATIVE) * (_TRUE_NEGATIVE + _FALSE_POSITIVE) * (_TRUE_NEGATIVE + _FALSE_NEGATIVE))

print("\n")
print(f"Accuracy: {_ACCURACY*100:.2f} %")
print(f"Sensitivity: {_SENSITIVITY*100:.2f} %")
print(f"Specificity: {_SPECIFICITY*100:.2f} %")
print(f"Precision: {_PRECISION*100:.2f} %")
print(f"Matthews Correlation Coefficient: {_MCC:.2f}")
