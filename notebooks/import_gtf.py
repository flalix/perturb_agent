
from pathlib import Path
import numpy as np
import pandas as pd
import gzip

'''
# GTF (this is the one infercnvpy's genomic_position_from_gtf parses)
wget https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz

# GFF3 (only if you want it for something else — infercnv doesn't need it)
# wget https://ftp.ensembl.org/pub/release-110/gff3/homo_sapiens/Homo_sapiens.GRCh38.110.gff3.gz
'''

root_refseq = Path('/home/flavio/uv/perturb_agent/data/colab/refseq/')
fname = 'Homo_sapiens.GRCh38.110.gtf.gz'
filename = root_refseq / fname


def parse_gtf_attrs(s):
    d = {}
    for field in s.strip().split(";"):
        field = field.strip()
        if not field: continue
        key, _, val = field.partition(" ")      # GTF: key "value"
        d[key] = val.strip().strip('"')
    return d

fname2 = 'hsa_GRCh38_110_gtf_simple.tsv'
filename2 = root_refseq / fname2

if filename2.exists():
    dfgtf = pd.read_csv(filename2, sep='\t')
    print(f"Read: {filename2}")
else:
    rows = []
    with gzip.open(filename, "rt") as fh:
        for line in fh:
            if line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene": continue
            attr = parse_gtf_attrs(f[8])
            name = attr.get("gene_name")             # GTF uses gene_name
            if name:
                rows.append((name, f[0], int(f[3]), int(f[4])))   # 4 values -> 4 columns

    dfgtf = pd.DataFrame(rows, columns=["gene_name","chromosome","start","end"])
    dfgtf["chromosome"] = "chr" + dfgtf["chromosome"].astype(str)
    dfgtf = dfgtf.drop_duplicates("gene_name").set_index("gene_name")

    dfgtf.to_csv(filename2, sep='\t',index=True)
    print(f"Saved: {filename2}")

print(dfgtf.shape)          # expect ~30-40k genes for GTF release 110/116
dfgtf.head(3)

