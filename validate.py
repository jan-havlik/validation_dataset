import requests
import json
import os

from math import sqrt

JWT_TOKEN = "XXX"

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

def parse_genome_browser_position(position: str):

    # split original gene into two halves
    # return list of gene parts consiting of [first_half, second_half, whole_sequence]

    sp = position.split(":")

    chr = sp[0]
    start = int(sp[1].split("-")[0].replace(',', ''))
    end = int(sp[1].split("-")[1].replace(',', ''))
   
    return [chr, start, end]


def download_seq(pos: str, genome: str = "hg19", out: str = None):
    
    # downloads sequence via genome browser api

    splitted_gene = parse_genome_browser_position(pos)
    chr, start, end = splitted_gene[0], splitted_gene[1], splitted_gene[2]

    r = requests.get(f'https://api.genome.ucsc.edu/getData/sequence?genome={genome};chrom={chr};start={start};end={end};')
    if r.status_code == 200:
        data = r.json()['dna'].upper()
        if out:
            with open(f'gene_data/{out}.fasta', 'w') as f:
                f.write(f'>{out}\n')
                f.write(data)


def download_dripc(position: str, genome: str = "hg19", out: str = None):
    
    # downloads DRIPc analysis via genome browser api

    splitted_gene = parse_genome_browser_position(position)
    chr, start, end = splitted_gene[0], splitted_gene[1], splitted_gene[2]
    
    dripc = requests.get(f'https://api.genome.ucsc.edu/getData/track?genome={genome};chrom={chr};start={start};end={end};track=hub_124277_DRIPc_peak_NT2_Sanz_2016')
    if dripc.status_code == 200:
        if out:
            with open(f'experimental_data/{out}.txt', 'w') as f_p:
                f_p.write(dripc.text)


def upload_seq(position: str, gene: str):
    
    # uploads given sequence to DNA Analyser web server

    headers = {"Authorization": JWT_TOKEN, "Content-Type": "application/json"}
    api_url = "https://bioinformatics.ibp.cz"

    splitted_gene = parse_genome_browser_position(position)
    chr, start, end = splitted_gene[0], splitted_gene[1], splitted_gene[2]
    data_in = open(f'gene_data/{gene}.fasta', 'r').read()

    seq = requests.post(f"{api_url}/api/sequence/import/text", headers=headers, data=json.dumps({
        "circular": False,
        "data": data_in,
        "format": "FASTA",
        "name": f"{gene.capitalize()}_{position}",
        "tags": [
            f"review"
        ],
        "type": "DNA"
    }))
    sequence_id = seq.json()["payload"]["id"]


_GENE_MAP = {
    "JTB": "chr1:153946746-153950164",
    "Immunoglobulin": "chr14:106235440-106237742",
    "MYC": "chr8:128745680-128755674",
    "RHOH": "chr4:40192674-40248587",
    "ACTB": "chr7:5566783-5603415",
    "FMR1": "chrX:146993470-147032645",
    "SNRPN": "chr15:25068795-25223870",
    "HK2": "chr2:75061109-75120486",
    "CIRH1A": "chr16:69165195-69265033",
    "APOE": "chr19:45409012-45412650",
    "FHIT": "chr3:59735037-61237133",
    "PPM1D": "chr17:58677545-58741849",
    "TP53": "chr17:7565098-7590856",
    "PBX1": "chr1:164524822-164868533",
}


def compare_results(pos, gene):

    positive_rloop_match = 0
    positive_dripc_match = 0
    tp = tn = fp = fn = 0
    
    def true_positive(alg, exp): return alg > 0 and exp > 0
    def true_negative(alg, exp): return alg == 0 and exp == 0
    def false_positive(alg, exp): return alg > 0 and exp == 0
    def false_negative(alg, exp): return alg == 0 and exp > 0


    rloopt = open(f'rlooptracker_data/{gene}.bed', 'r').read().split('\n')[2:-1]
    dripc = json.load(open(f'experimental_data/{gene}.txt', 'r'))["hub_124277_DRIPc_peak_NT2_Sanz_2016"]


    # detecs positive strand r-loops only
    positive_rloop_match = len([(int(rloop.split(" ")[1]), int(rloop.split(" ")[2])) for rloop in rloopt if int(rloop.split(" ")[4]) >= 0])
    positive_dripc_match = len([(int(rloop["chromStart"]), int(rloop["chromEnd"])) for rloop in dripc if rloop["strand"] == "+"])

    # detects negative strand r-loops only
    negative_rloop_match = len([(int(rloop.split(" ")[1]), int(rloop.split(" ")[2])) for rloop in rloopt if int(rloop.split(" ")[4]) < 0])
    negative_dripc_match = len([(int(rloop["chromStart"]), int(rloop["chromEnd"])) for rloop in dripc if rloop["strand"] == "-"])

    tp = tp + 1 if true_positive(positive_rloop_match, positive_dripc_match) else tp
    tp = tp + 1 if true_positive(negative_rloop_match, negative_dripc_match) else tp
    
    tn = tn + 1 if true_negative(positive_rloop_match, positive_dripc_match) else tn
    tn = tn + 1 if true_negative(negative_rloop_match, negative_dripc_match) else tn
    
    fp = fp + 1 if false_positive(positive_rloop_match, positive_dripc_match) else fp
    fp = fp + 1 if false_positive(negative_rloop_match, negative_dripc_match) else fp

    fn = fn + 1 if false_negative(positive_rloop_match, positive_dripc_match) else fn
    fn = fn + 1 if false_negative(negative_rloop_match, negative_dripc_match) else fn

    return {
        "tp": tp,
        "tn": tn,
        "fp": fp,
        "fn": fn
    }


""" 
///////////////////////////////
        VALIDATION CALC         
///////////////////////////////
"""

TP = TN = FP = FN = 0
for gene, position in _GENE_MAP.items():
    download_seq(position, out=gene)
    download_dripc(position, out=gene)
    upload_seq(position, gene)
    res = compare_results(position, gene)


    TP += res["tp"]
    TN += res["tn"]
    FP += res["fp"]
    FN += res["fn"]

accuracy = (TP + TN) / (TP + TN + FP + FN)
specificity = TN / (TN + FP)
sensitivity = TP / (TP + FN)
precision = TP / (TP + FP)
mcc = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

print(f"Accuracy: {accuracy*100.0:.2f} %")
print(f"Specificity: {specificity*100.0:.2f} %")
print(f"Sensitivity: {sensitivity*100.0:.2f} %")
print(f"Precision: {precision*100.0:.2f} %")
print(f"Matthews Correlation Coefficient: {mcc:.2f}")
