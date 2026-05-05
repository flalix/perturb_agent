#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2024/04/11
# Udated  on 2024/05/03
# @author: Flavio Lichtenstein
# @local: Bioinformatics: CENTD/Molecular Biology; Instituto Butatan

from fileinput import filename
import json
import os

import numpy as np
import pybiopax
import pandas as pd
import requests

from pathlib import Path
from typing import List


from pybiopax import model_from_owl_file
from pybiopax.biopax import Pathway, BioPaxModel

from collections import deque

from libs.Basic import isfloat, pdreadcsv, pdwritecsv, create_dir


class Reactome(object):
    def __init__(self, ROOT0: Path):

        self.root_src = create_dir(ROOT0, "src")
        self.root_owl = create_dir(ROOT0, "owl")
        self.root_reactome = create_dir(ROOT0, "data/reactome")
        self.root_json = self.root_reactome / "pathway_json"

        self.fname_reactome = "reactome_pathways_human.tsv"
        self.fname_reactome_2024 = "Reactome_Pathways_2024.tsv"
        self.fname_reactome_abstract = "reactome_pathways_abstract.tsv"
        self.fname_reactome_gmt = "reactome_pathways_human_gmt.tsv"

        self.fname_human_owl = "Homo_sapiens.owl"
        self.hsa_model = None
        self.fname_owl = "%s_level3.owl"


        self.dfr, self.df_gmt, self.dfr_merge = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
        self.pathway_list = []
        self.n_pathways = 0

        self.script_curl = "curl -X GET --header 'Accept: application/json' 'https://reactome.org/ContentService/data/query/%s' --silent > %s"

        self.url_owl_level3 = "https://reactome.org/download/current/biopax/%s_level3.owl"

    def refresh_reactome_table(
        self, pathw_list: List, force: bool = False, verbose: bool = True
    ) -> bool:

        dfr = self.open_reactome_abstract(verbose=False)
        self.dfr = dfr

        if dfr is None:
            redo = True
        else:
            pathw_list = [x for x in pathw_list if x not in dfr.pathway_id]
            redo = len(pathw_list) > 0

        if redo:
            self.get_many_reactome_content_service_pathway(pathw_list, verbose=verbose)
            ret = self.build_reactome_table(dfr, force=force, verbose=verbose)

            if verbose:
                print("reactome table was updated.")
        else:
            if verbose:
                print("reactome table is Ok.")
            ret = False

        # return  if has reviewed
        return ret

    def get_many_reactome_content_service_pathway(
        self, pathw_list: List, force: bool = False, verbose: bool = False
    ):

        for pathway_id in pathw_list:
            # if not verbose: print(".", end='')
            self.get_reactome_content_service_pathway(pathway_id, force=force, verbose=verbose)

    def get_reactome_content_service_pathway(
        self, pathway_id: str, force: bool = False, verbose: bool = False
    ) -> bool:
        pathway_id_json = f"{pathway_id}.json"
        fullname = os.path.join(self.root_json, pathway_id_json)

        if os.path.exists(fullname) and not force:
            if verbose:
                print(f"File already exists: {fullname}")
            return True

        command = self.script_curl % (pathway_id, fullname)
        ret = os.system(command)
        if verbose:
            if ret == 0:
                print(f"File downloaded {ret}: {pathway_id} -> {fullname}")
            else:
                print(f"Error: could not download {ret}: {pathway_id} -> {fullname}")

        if os.path.exists(fullname) and os.path.getsize(fullname) == 0:
            os.remove(fullname)
            print(f"Error: size zero json file: {pathway_id} -> {fullname}")

        return ret == 0

    def open_reactome(self, verbose: bool = False) -> pd.DataFrame:

        fullname = os.path.join(self.root_reactome, self.fname_reactome)

        if not os.path.exists(fullname):
            self.pathway_list = []
            self.dfr = pd.DataFrame()
            self.n_pathways = 0

            print(f"Could not find: {fullname}")
            return pd.DataFrame()

        """
            bpx.dfr['pathway'] = [x.strip() for x in bpx.dfr.pathway]
            pdwritecsv(bpx.dfr, bpx.reactome.fname_reactome, bpx.reactome.root_reactome, verbose=True)
        """
        dfr = pdreadcsv(self.fname_reactome, self.root_reactome, verbose=verbose)

        self.dfr = dfr
        self.pathway_list = list(np.unique(dfr.pathway))
        self.n_pathways = len(self.pathway_list)

        return dfr

    def open_reactome_abstract(self, verbose: bool = False) -> pd.DataFrame:

        fullname = os.path.join(self.root_reactome, self.fname_reactome_abstract)

        if not os.path.exists(fullname):
            print(f"Could not find: {fullname}")
            return pd.DataFrame()

        df = pdreadcsv(self.fname_reactome_abstract, self.root_reactome, verbose=verbose)
        return df

    def open_reactome_gmt(self, verbose: bool = False) -> pd.DataFrame:
        fullname = os.path.join(self.root_reactome, self.fname_reactome_gmt)

        if not os.path.exists(fullname):
            print(f"There is no reactome gst file: '{fullname}'")
            return pd.DataFrame()

        df_gmt = pdreadcsv(self.fname_reactome_gmt, self.root_reactome, verbose=verbose)
        self.df_gmt = df_gmt

        return df_gmt

    def merge_reactome_table_gmt(self, verbose: bool = False) -> bool:
        self.dfr_merge = pd.DataFrame()

        dft = self.open_reactome_abstract(verbose=verbose)
        dfg = self.open_reactome_gmt(verbose=verbose)

        if dft is None or dft.empty:
            return False
        if dfg is None or dfg.empty:
            return False

        cols = ["pathway_id", "genes", "n"]
        dfg = dfg[cols].copy()
        dfg.columns = ["pathway_id", "genes_pathway", "ngenes_pathway"]

        self.dfr_merge = pd.merge(dft, dfg, how="inner", on="pathway_id")

        return True

    def build_reactome_table(self, dfr: pd.DataFrame, force: bool = False, verbose: bool = False) -> bool:

        if force:
            dfr = pd.DataFrame()
            files = [x for x in os.listdir(self.root_json) if x.endswith(".json")]

            if len(files) == 0:
                print(f"There are no json files in '{self.root_json}'")
                return False

        else:
            pathway_id_list = dfr.pathway_id.to_list()
            files = [
                x
                for x in os.listdir(self.root_json)
                if x.endswith(".json")
                if x + ".json" not in pathway_id_list
            ]

            if len(files) == 0:
                if verbose:
                    print("There are no json files to update")
                return False

        if verbose:
            print(f"Updating {len(files)} json files.")

        dicf = {}
        i = -1

        for fname in files:
            fullname = os.path.join(self.root_json, fname)

            if not os.path.exists(fullname):
                print(f"Could not find: {fullname}")
                continue

            try:
                f = open(fullname)
                dic = json.load(f)
            except ValueError:
                dic = None
                print(f"Could not open: {fullname}")
                continue

            try:
                ret_error = dic["code"]
                os.remove(fullname)
            except ValueError:
                ret_error = None
                pass

            if ret_error is not None:
                try:
                    ret_msg = dic["messages"]
                except ValueError:
                    ret_msg = "???"
                print(f"Error in json download: {fullname} - ret_error={ret_error} {ret_msg}")

                continue

            i += 1
            dicf[i] = {}
            dic2 = dicf[i]

            for key, val in dic.items():
                dic2[key] = val

                if key == "summation":
                    all_text = ""
                    summaries = val
                    if len(summaries) > 0:
                        for dic3 in summaries:
                            if all_text != "":
                                all_text += "\n\n"
                            all_text += dic3["text"]

                    dic2["abstract"] = all_text

        if len(dicf) == 0:
            print("Nothing updated.")
            return False

        df = pd.DataFrame(dicf).T
        cols = list(df.columns)

        # manually add R-HSA-71406
        special_pathwya_id = "R-HSA-71406"

        if dfr.empty:
            import_special = True
        else:
            import_special = False if special_pathwya_id in dfr.pathway_id else True

        if import_special:
            i += 1
            dicf[i] = {}
            dic2 = dicf[i]

            dic2[cols[1]] = "Pyruvate metabolism and Citric Acid (TCA) cycle"
            dic2[cols[2]] = special_pathwya_id
            dic2["abstract"] = (
                "In the citric acid or tricarboxylic acid (TCA) cycle, the acetyl group of acetyl CoA (derived primarily from oxidative decarboxylation of pyruvate, beta-oxidation of long-chain fatty acids, and catabolism of ketone bodies and several amino acids) can be completely oxidized to CO2 in reactions that also yield one high-energy phosphate bond (as GTP or ATP) and four reducing equivalents (three NADH + H+, and one FADH2). Then, the electron transport chain oxidizes NADH and FADH2 to yield nine more high-energy phosphate bonds (as ATP). All reactions of the citric acid cycle take place in the mitochondrion."
            )

            # again
            df = pd.DataFrame(dicf).T

        cols[1] = "pathway"
        cols[2] = "pathway_id"
        df.columns = cols

        if not dfr.empty:
            if verbose:
                print(f"There are {len(dfr)} dfr regs, incrementing {len(df)}.")

            cols = list(dfr.columns)
            # df with the same columns
            df = df[cols]

            df = pd.concat([dfr, df])

        goods = [
            True if isfloat(df.iloc[i]["dbId"]) or pd.isnull(df.iloc[i]["dbId"]) else False
            for i in range(len(df))
        ]
        df = df[goods]

        df = df.drop_duplicates("pathway_id")
        df.index = np.arange(0, len(df))

        pdwritecsv(df, self.fname_reactome_abstract, self.root_reactome, verbose=verbose)
        self.dfr_merge = df
        return True


    def download_reactome_owl(self, pathway_id: str,  force: bool = False,  timeout: int = 60) -> dict:

        pathway_id = str(pathway_id).strip()

        fname_owl = self.fname_owl%(pathway_id)
        filename = self.root_reactome / fname_owl

        if filename.exists() and filename.stat().st_size > 0 and not force:
            return {
                "pathway_id": pathway_id,
                "filename": str(filename),
                "status": "skipped_exists",
                "url": None,
                "http_status": None,
                "n_bytes": filename.stat().st_size,
                "msg": 'OWL file already exists',
            }

        url = self.url_owl_level3%(pathway_id)
        error_msg = ""
        try:
            r = requests.get(url, timeout=timeout)

            if r.status_code != 200:
                error_msg = f"{url} -> HTTP {r.status_code}"
            else:
                content = r.content

                if len(content) < 500:
                    error_msg = f"{url} -> too small: {len(content)} bytes"
                else:
                    with open(filename, "wb") as f:
                        f.write(content)

                    return {
                        "pathway_id": pathway_id,
                        "filename": str(filename),
                        "status": "downloaded",
                        "url": url,
                        "http_status": r.status_code,
                        "n_bytes": len(content),
                        "msg": 'OWL file downloaded successfully',
                    }

        except Exception as e:
            error_msg = f"Error occurred while downloading {pathway_id}: {repr(e)}"

        return {
            "pathway_id": pathway_id,
            "filename": str(filename),
            "status": "failed",
            "url": None,
            "http_status": None,
            "n_bytes": 0,
            "msg": error_msg
        }


    def open_human_owl(self, force: bool = False):
        if self.hsa_model is None or force:
            filename = self.root_reactome / self.fname_human_owl
            self.hsa_model = model_from_owl_file(str(filename))

        return self.hsa_model


    def get_reactome_ids_given_obj(self, obj):
        ids = set()

        # Sometimes stable IDs appear in uid
        uid = getattr(obj, "uid", None)
        if isinstance(uid, str) and uid.startswith("R-HSA-"):
            ids.add(uid)

        # Usually stable IDs are in xref
        for x in getattr(obj, "xref", []) or []:
            xid = getattr(x, "id", None)
            db = getattr(x, "db", None)

            if xid and str(xid).startswith("R-HSA-"):
                ids.add(str(xid))

            if db and "Reactome" in str(db) and xid:
                ids.add(str(xid))

        return ids

    def get_pathways_from_model(self):
        
        pathways = [ obj for obj in self.hsa_model.objects.values()  if type(obj).__name__ == "Pathway"  ]

        pathway_index = {}

        for p in pathways:
            ids = self.get_reactome_ids_given_obj(p)

            for rid in ids:
                pathway_index[rid] = p

        return pathways, pathway_index
    
    def pathways_to_df(self, pathways) -> pd.DataFrame:
        dic={}
        icount=-1

        for p in pathways:
            ids = self.get_reactome_ids_given_obj(p)

            for rid in ids:
                icount += 1
                dic[icount] = {}

                dic2 = dic[icount]

                dic2['pathway_object'] = p
                dic2['rid'] = rid
                dic2['is_ID'] = str(rid).startswith('R-HSA-')

        dfo = pd.DataFrame.from_dict(dic, orient='index')
        return dfo
                

    def is_biopax_object(self, x):
        return hasattr(x, "uid") and hasattr(x, "__dict__")


    def get_pathway_names(self, p, dfo: pd.DataFrame) -> tuple[str, str, str]:

        dfa = dfo[
            (dfo["pathway_object"].map(lambda x: x is p)) &
            (dfo["is_ID"])
        ]

        if dfa.empty:
            pathway_id = ""
        else:
            row = dfa.iloc[0]
            pathway_id = row.rid

        uid = p.uid
        pathway = p.display_name
    
        return uid, pathway_id, pathway

    def collect_referenced_objects(self, pthw_obj, include_inverse=False) -> dict:
        """
        Recursively collect all BioPAX objects reachable from pthw_obj.

        include_inverse=False is safer because inverse links such as:
            _participant_of
            _controller_of
        may pull in too much of the full model.

        For rebuilding one pathway OWL, usually use forward references only.
        """

        if pthw_obj is None:
            print("Pathway is None")
            return {}
        
        if pthw_obj.uid is None:
            print("Pathway has no stable ID")
            return {}
                    

        sub_objects = {}
        queue = deque([pthw_obj])

        while queue:
            obj = queue.popleft()

            uid = getattr(obj, "uid", None)
            if uid is None:
                continue

            if uid in sub_objects:
                continue

            sub_objects[uid] = obj

            for attr, value in obj.__dict__.items():

                # Skip inverse/back-reference fields unless explicitly requested
                if not include_inverse and attr.startswith("_"):
                    continue

                if value is None:
                    continue

                if self.is_biopax_object(value):
                    queue.append(value)

                elif isinstance(value, (list, tuple, set)):
                    for item in value:
                        if self.is_biopax_object(item):
                            queue.append(item)

                elif isinstance(value, dict):
                    for item in value.values():
                        if self.is_biopax_object(item):
                            queue.append(item)
        
        return sub_objects


    def save_owl_model(self, pathway_id:str, pathway:str, sub_objects: dict,
                       force:bool=False, verbose:bool=False) -> bool:
        
        fname_owl = self.fname_owl%(pathway_id)
        filename = self.root_owl / fname_owl

        if filename.exists() and not force:
            if verbose: print(f"File {filename} already exists.")
            return True

        try:
            sub_model = BioPaxModel(objects=sub_objects)
            pybiopax.model_to_owl_file(sub_model, filename)
            if verbose: print(f"Submodel '{pathway}' saved to {filename}")
            ret = True
        except Exception as e:
            print(f"Error saving submodel to {filename}: {e}")
            ret = False

        return ret
    
