#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
'''
# @Updated on 2026/03/20
# @Created on 2026/03/20
# @author: Flavio Lichtenstein
# @local: Home sweet home
'''

# import requests
# import json
import os, sys
import pandas as pd

from pathlib import Path

ROOT = Path().resolve().parent.parent
SRC = os.path.join(ROOT, "src")

if str(SRC) not in sys.path:
    sys.path.append(str(SRC))

print("ROOT:", ROOT)
print("SRC added:", SRC)

from libs.tcga_gdc_lib import *
from libs.Basic import *

root_data = "/media/flavio/36873e3e-7941-48d7-aecb-45900ef92659/colaboracoes/tcga"
case = "Breast Cancer"

def tool2(pathway: str) -> pd.DataFrame:
    gdc = GDC(case=case, root_data=root_data)

    # barcode case

    case_uuid2 = 'e3935ce4-64d3-4a66-ba11-d308b844b410'
    barcode = 'TCGA-E9-A5FL'

    case_uuid = gdc.get_case_uuid(barcode)

    ## expression table files
    response = gdc.get_files(case_uuid, data_type="Gene Expression Quantification")

    gdc.show_fnames()

    verbose=True
    force=False

    ret = gdc.download_file(force=force, verbose=verbose)

    if ret:
        filename = gdc.gdc_filename
    else:
        filename = 'not_found.tsv'
        return None

    df2 = gdc.clean_gdc_tsv()

    return df2




