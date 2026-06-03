#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2026/05/04
# Udated  on 2026/05/04
# @author: Flavio Lichtenstein
# @local: Home sweet home


import json
import math
import numpy as np
import os
import re
import subprocess
from os import path
from collections import Counter
import pandas as pd
import threading
import webbrowser
import socket
from pathlib import Path

import dash
from dash import html, dcc, dash_table, Input, Output, State, ctx, no_update
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
import networkx as nx
from rdflib import RDF, Graph, Namespace

# import py4cytoscape as p4c
from libs.Basic import create_dir, pdreadcsv, download_url_file
from libs.open_target_lib import OpenTarget


class DASH_CYTO(object):
    def __init__(self, root0: Path, root0_data: Path, dflfc_ori: pd.DataFrame, 
                 found_degs: list, pathway_genes: list, psi_id: str = "",
                 lfc_cutoff: float=1.0, fdr_cutoff: float=0.05):

        self.GENE_COLS = ["gene_id", "symbol", "gene_type"]

        self.root0 = root0
        self.root0_data = root0_data

        self.root_src = create_dir(root0, "src")
        self.root_styles = create_dir(self.root_src, "styles")

        self.root_colab= create_dir(root0_data, "colab")
        self.root_owl = create_dir(self.root_colab, "owl")
        self.root_ncbi = create_dir(self.root_colab, "ncbi")

        self.ot = OpenTarget(root_colab=self.root_colab)
        self.psi_id = psi_id
        self.disease = psi_id.replace("TCGA-", "")

        self.dflfc_ori = dflfc_ori
        self.lfc_cutoff = lfc_cutoff
        self.fdr_cutoff = fdr_cutoff
        self.found_degs = found_degs
        self.pathway_genes = pathway_genes
        self.settings_file = self.root_owl /'graph_settings.json'

        self.fname_pos = "positions_%s.json"

        self.url_reactome_image = "https://reactome.org/ContentService/exporter/diagram/%s.png?quality=7"
        self.url_reactome_owl   = "https://reactome.org/ReactomeRESTfulAPI/RESTfulWS/biopaxExporter/Level3/%s"        

        # has columns: ['ensembl_id', 'symbol', 'name', 'uniprot_id', 'ncbi_gene_id', 'synonyms', 'refseq_summary']
        self.fname_hugo       = "hugo_gene_table_refseq_uniprot.tsv"
        self.fname_gene_alias = "hugo_gene_alias_table.tsv"

        self.pathway_id = ""
        self.pathway = ""
        self.rdf = Graph()
        self.MAX_LEN_LABEL = 25

        BP = Namespace("http://www.biopax.org/release/biopax-level3.owl#")
        self.BP = BP

        self.G = nx.DiGraph()
        self.saved_positions = {}

        self.classes = [
            BP.Pathway,
            BP.BiochemicalReaction,
            BP.Catalysis,
            BP.Control,
            BP.Protein,
            BP.Complex,
            BP.SmallMolecule,
            BP.PhysicalEntity,
        ]

        self.relations = [
            BP.pathwayComponent,
            BP.left,
            BP.right,
            BP.controller,
            BP.controlled,
            BP.participant,
            BP.component,
        ]

        self.load_gene_annotation_table()
        self.load_gene_alias_table()

        self.START_PORT = 8051
        self.END_PORT = 8070

        if not self.is_render():
            self.kill_dash_ports()
        else:
            print("Render detected: skipping local Dash port cleanup")


        self.hub_list = None
        self.hub_top_n = None

        self.source_node_list = None
        self.sink_node_list = None

    def is_render(self) -> bool:
        return os.environ.get("RENDER") == "true" or os.environ.get("RENDER_SERVICE_ID") is not None

    def reset_graph(self):
        self.G = nx.DiGraph()
        self.saved_positions = {}
        self.layout_changed = False

        self.hub_list = None
        self.hub_top_n = None
        
        self.source_node_list = None
        self.sink_node_list = None        

    def read_owl(self, pathway_id: str, pathway: str, verbose: bool = False) -> bool:
                            
        self.reset_graph()

        self.rdf = Graph()
        self.pathway_id = pathway_id
        self.pathway = pathway

        #------------ first image ------------------------
        fname_png = f"{pathway_id}.png"
        filename = self.root_owl / fname_png

        if not filename.exists():
            url = self.url_reactome_image%(pathway_id)
            ret = download_url_file(url, fname=fname_png, root_file=self.root_owl, verbose=True)

        #------------ next, level3.owl ------------------
        fname_owl = f"{pathway_id}_level3.owl"
        filename = self.root_owl / fname_owl

        if not filename.exists():
            ID = pathway_id.replace("R-HSA-", "")
            url = self.url_reactome_owl%(ID)
            ret = download_url_file(url, fname=fname_owl, root_file=self.root_owl, verbose=True)
            if not ret: 
                return False

        if not filename.exists():
            print(f"File not found: {filename}")
            return False

        if verbose:
            print(f"Reading OWL file {filename} - pathway: '{pathway}'")

        #---------- read owl ------------------------
        try:
            self.rdf.parse(filename, format="xml")
        except Exception as e:
            print(f"Error parsing OWL file: {e}")

        #if verbose:
        #    print(f"Processing {len(list(self.rdf.subjects()))} RDF triples")

        for cls in self.classes:
            for node in self.rdf.subjects(RDF.type, cls):
                node_id = self.short(node)
                self.G.add_node(node_id, label=self.get_name(node, self.MAX_LEN_LABEL), biopax_type=self.short(cls))

        for rel in self.relations:
            for s, o in self.rdf.subject_objects(rel):
                s_id = self.short(s)
                o_id = self.short(o)

                if s_id in self.G.nodes and o_id in self.G.nodes:
                    self.G.add_edge(s_id, o_id, interaction=self.short(rel))

        return True

    def short(self, x) -> str:
        return str(x).split("#")[-1].split("/")[-1]

    def get_name(self, x, max_len: int | None = None) -> str:

        for prop in [self.BP.displayName, self.BP.standardName, self.BP.name]:
            value = next(self.rdf.objects(x, prop), None)
            if value:
                value = str(value)
                return value if max_len is None else value[:max_len]

        return self.short(x)

    def load_positions(self, pathway_id: str):

        fname = self.fname_pos % (pathway_id)
        filename = self.root_owl / fname

        dic_posi = {}
        try:
            if os.path.exists(filename):
                with open(filename) as fh:
                    print(f"Data loaded from {filename}")
                    dic_posi = json.load(fh)
        except Exception as e:
            print(f"Error loading positions from {filename}: {e}")

        return dic_posi


    def build_lfc_lookup(self) -> dict:
        if self.dflfc_ori is None or self.dflfc_ori.empty:
            return {}

        df = self.dflfc_ori.copy()

        df["symbol"] = df["symbol"].astype(str).str.upper()
        df["lfc"] = pd.to_numeric(df["lfc"], errors="coerce")
        df["fdr"] = pd.to_numeric(df["fdr"], errors="coerce")

        # If duplicated symbols, keep strongest absolute LFC
        df = (
            df.sort_values("abs_lfc", ascending=False)
            .drop_duplicates("symbol")
        )

        return {
            row["symbol"]: {
                "lfc": float(row["lfc"]),
                "fdr": float(row["fdr"]) if pd.notna(row["fdr"]) else None,
                "abs_lfc": float(row["abs_lfc"]),
            }
            for _, row in df.iterrows()
        }


    def get_lfc_bin(self, lfc: float | None) -> str:
        if lfc is None:
            return "none"

        try:
            lfc = float(lfc)
        except Exception:
            return "none"

        if lfc >= 3:
            return "up_3"
        elif lfc >= 2:
            return "up_2"
        elif lfc >= 1:
            return "up_1"
        elif lfc >= 0.4:
            return "+weak"
        elif lfc <= -3:
            return "down_3"
        elif lfc <= -2:
            return "down_2"
        elif lfc <= -1:
            return "down_1"
        elif lfc <= -0.4:
            return "-weak"

        return "~zero"
        

    def nx_to_cytoscape_elements(self, saved_positions=None) -> list:
        elements = []

        self.dic_lfc_lookup = self.build_lfc_lookup()

        if saved_positions is None:
            saved_positions = self.load_positions(self.pathway_id)

        self.saved_positions = saved_positions

        """
        missing = set(saved_positions.keys()) - set(self.G.nodes())
        if missing:
            print(f"⚠️ positions for {len(missing)} nodes not used (graph changed)")
            for node_id in missing:
                print(f"  - {node_id}")
        """

        for node_id, data in self.G.nodes(data=True):

            label = data.get("label", node_id)

            # --------------------------------------------------
            # NEW: annotate Reactome node label using alias table
            # --------------------------------------------------
            annot = self.get_gene_annotation_for_node({
                "id": node_id,
                "label": label,
                **data,
            })

            symbol = annot.get("symbol", "")
            ensembl_id = annot.get("ensembl_id", "NA")
            name = annot.get("name", "NA")
            uniprot_id = annot.get("uniprot_id", "NA")
            refseq_summary = annot.get("refseq_summary", "NA")
            n_matches = annot.get("n_matches", 0)
            matches = annot.get("matches", [])

            # Use official symbol for DEG lookup when available
            if symbol not in [None, ""]:
                lookup_symbol = str(symbol).upper()
            else:
                lookup_symbol = str(label).upper()

            dic_symbol = self.dic_lfc_lookup.get(lookup_symbol)

            if dic_symbol is not None:
                lfc = dic_symbol["lfc"]
                fdr = dic_symbol["fdr"]
            else:
                lfc = None
                fdr = None

            elem = {
                "data": {
                    "id": node_id,
                    "label": label,
                    "biopax_type": data.get("biopax_type", "Unknown"),

                    # --------------------------------------------------
                    # NEW: keep gene annotation inside Cytoscape data
                    # --------------------------------------------------
                    "symbol": symbol,
                    "ensembl_id": ensembl_id,
                    "name": name,
                    "uniprot_id": uniprot_id,
                    "refseq_summary": refseq_summary,
                    "n_matches": n_matches,
                    "matches": matches,

                    # DEG fields
                    "lfc": lfc,
                    "FDR": fdr,
                    "abs_lfc": abs(lfc) if lfc is not None else None,
                    "lfc_bin": self.get_lfc_bin(lfc),
                    "is_deg_gene": lfc is not None,
                }
            }

            if node_id in saved_positions:
                elem["position"] = saved_positions[node_id]

            elements.append(elem)

        for source, target, data in self.G.edges(data=True):
            elements.append(
                {
                    "data": {
                        "source": source,
                        "target": target,
                        "interaction": data.get("interaction", ""),
                    }
                }
            )

        return elements

    def extract_positions(self, elements):
        return {
            e["data"]["id"]: e["position"]
            for e in elements
            if "data" in e and "id" in e["data"] and "position" in e
        }

    def positions_changed(self, new: dict, tol_px: float = 2.0):
        old = self.saved_positions or {}

        if not old:
            return True

        if old.keys() != new.keys():
            print("Warning: different node sets.")
            return True

        for k in old:
            dx = new[k]["x"] - old[k]["x"]
            dy = new[k]["y"] - old[k]["y"]

            if math.hypot(dx, dy) > tol_px:
                return True

        return False

    def save_positions_if_changed(self, elements) -> bool:

        new_pos = self.extract_positions(elements)
            
        print(f"save_positions_if_changed() -> {self.positions_changed(new_pos)} {self.layout_changed}")

        if not self.positions_changed(new_pos) and not self.layout_changed:
            return False  # no change

        fname = self.fname_pos % (self.pathway_id)
        filename = self.root_owl / fname

        try:
            with open(filename, "w") as f:
                json.dump(new_pos, f, indent=2)

            self.saved_positions = new_pos  # update global here
            self.layout_changed = False
        except Exception as e:
            print(f"Error: saving Graph positions: {e}")
            return False

        return True
    
    #============ Hugo, Ensbembl, UNIPROT gene table ===========

    def load_gene_annotation_table(self, verbose: bool = False) -> pd.DataFrame:
        """
        Load HGNC/HUGO + RefSeq + UniProt annotation table.

        Expected columns:
            ensembl_id, symbol, name, uniprot_id,
            ncbi_gene_id, synonyms, refseq_summary
        """

        filename = self.root_ncbi / self.fname_hugo
        if not filename.exists():
            raise FileNotFoundError(f"Could not find gene annotation file: {filename}")

        df = pdreadcsv(self.fname_hugo, self.root_ncbi, verbose=verbose)

        required_cols = [
            "ensembl_id",
            "symbol",
            "name",
            "uniprot_id",
            "refseq_summary",
        ]

        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            raise ValueError(f"Missing columns in {filename}: {missing}")

        # Keep only useful columns for node annotation
        df = df[required_cols].copy()

        # Safety normalization
        df = df.fillna("")
        df["ensembl_id"] = df["ensembl_id"].astype(str).str.strip()
        df["symbol"] = df["symbol"].astype(str).str.strip()
        df["name"] = df["name"].astype(str).str.strip()
        df["uniprot_id"] = df["uniprot_id"].astype(str).str.strip()
        df["refseq_summary"] = df["refseq_summary"].astype(str).str.strip()

        df["ensembl_id_base"] = df["ensembl_id"].str.replace(
            r"\.\d+$",
            "",
            regex=True,
        )

        self.gene_annot_df = df
        if verbose:
            print(">>> df hugo", df.columns)

        # Fast dictionaries
        self.gene_annot_by_ensembl = (
            df[df["ensembl_id_base"] != ""]
            .drop_duplicates("ensembl_id_base")
            .set_index("ensembl_id_base")
            .to_dict(orient="index")
        )

        self.gene_annot_by_symbol = (
            df.dropna(subset=["symbol"])
            .assign(symbol=lambda x: x["symbol"].astype(str).str.strip())
            .drop_duplicates("symbol")
            .set_index("symbol")
            .to_dict(orient="index")
        )

        if verbose:
            print(f"Loaded gene annotation table: {filename}")
            print(f"Rows: {len(df):,}")
            print(f"Symbols: {len(self.gene_annot_by_symbol):,}")
            print(f"Ensembl IDs: {len(self.gene_annot_by_ensembl):,}")

        return df

    def load_gene_alias_table(self, verbose: bool = False) -> pd.DataFrame:
        df = pdreadcsv(self.fname_gene_alias, self.root_ncbi).fillna("")

        required_cols = [
            "alias",
            "alias_upper",
            "symbol",
            "ensembl_id",
            "name",
            "uniprot_id",
            "refseq_summary",
        ]

        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            raise ValueError(f"Missing columns in alias table: {missing}")

        df = df[required_cols].copy()

        df["alias"] = df["alias"].astype(str).str.strip()
        df["alias_upper"] = df["alias_upper"].astype(str).str.strip().str.upper()
        df["symbol"] = df["symbol"].astype(str).str.strip()

        # Important: do NOT drop duplicates only by alias_upper
        # because that would remove ties like TSP1 -> THBS1 / PRSS55.
        df = df.drop_duplicates(
            subset=["alias_upper", "symbol", "ensembl_id"]
        )

        self.gene_alias_df = df.set_index("alias_upper", drop=False).sort_index()

        if verbose:
            print(f"Loaded alias rows: {len(self.gene_alias_df):,}")
            print(f"Unique aliases: {self.gene_alias_df.index.nunique():,}")

        return self.gene_alias_df

    def lookup_gene_alias(self, label: str) -> list[dict]:
        if not hasattr(self, "gene_alias_df"):
            print("WARNING: gene_alias_df not loaded")
            return []

        if label in [None, "", "NA"]:
            return []

        key = str(label).strip().upper()

        try:
            hit = self.gene_alias_df.loc[key]
        except KeyError:
            return []

        # One match: pandas returns Series
        if isinstance(hit, pd.Series):
            return [
                hit.drop(labels=["alias_upper"], errors="ignore").to_dict()
            ]

        # Multiple matches: pandas returns DataFrame
        return hit.drop(columns=["alias_upper"], errors="ignore").to_dict(
            orient="records"
        )


    def get_aliases_for_symbol(self, symbol: str | None) -> str:
        if not hasattr(self, "gene_alias_df"):
            return "NA"

        if symbol in [None, "", "NA"]:
            return "NA"

        df = self.gene_alias_df

        # If gene_alias_df is indexed by alias_upper, symbol is still a column
        hits = df[df["symbol"].astype(str).str.upper() == str(symbol).upper()]

        if hits.empty:
            return "NA"

        aliases = (
            hits["alias"]
            .dropna()
            .astype(str)
            .str.strip()
        )

        aliases = sorted(set(a for a in aliases if a))

        return "; ".join(aliases) if aliases else "NA"

    def get_gene_annotation_for_node(self, node_data: dict) -> dict:
        empty = {
            "ensembl_id": "NA",
            "symbol": "NA",
            "name": "NA",
            "uniprot_id": "NA",
            "refseq_summary": "NA",
            "synonyms": "NA",
            "matches": [],
            "n_matches": 0,
        }

        label = (
            node_data.get("label")
            or node_data.get("symbol")
            or node_data.get("gene_symbol")
            or node_data.get("name")
            or ""
        )

        label = str(label).strip()

        if not label:
            return empty

        matches = self.lookup_gene_alias(label)

        if not matches:
            # print(f"No gene annotation match for label: {label}")
            return empty

        first = dict(matches[0])
        first["matches"] = matches
        first["n_matches"] = len(matches)

        return first



    def none_one_multiple_alias(self, info_data: dict):

        # print(">>> info_data", info_data)

        if info_data.get("symbol", "NA") == "NA":
            return None
        
        matches = info_data.get("matches", [])
        n_matches = info_data.get("n_matches", 0)

        # print(">>> n_matches", n_matches)

        if n_matches > 1:
            alias_note = html.Div(
                [
                    html.P(
                        [
                            html.B("Alias warning: "),
                            f'{info_data["label"]} maps to {n_matches} genes.',
                        ],
                        style={
                            "color": "#92400e",
                            "fontWeight": "700",
                            "marginBottom": "6px",
                        },
                    ),
                    html.Table(
                        [
                            html.Thead(
                                html.Tr(
                                    [
                                        html.Th("Symbol"),
                                        html.Th("Ensembl ID"),
                                        html.Th("Name"),
                                    ]
                                )
                            ),
                            html.Tbody(
                                [
                                    html.Tr(
                                        [
                                            html.Td(str(m.get("symbol", ""))),
                                            html.Td(str(m.get("ensembl_id", ""))),
                                            html.Td(str(m.get("name", ""))),
                                        ]
                                    )
                                    for m in matches
                                ]
                            ),
                        ],
                        style={
                            "width": "100%",
                            "fontSize": "12px",
                            "borderCollapse": "collapse",
                        },
                    ),
                ],
                style={
                    "backgroundColor": "#fef3c7",
                    "padding": "8px",
                    "borderRadius": "8px",
                    "marginBottom": "10px",
                },
            )
        elif info_data["label"] != info_data["symbol"]:
            alias_note = html.P(
                [
                    html.B("Alias mapping: "),
                    f'{info_data["label"]} → {info_data["symbol"]}',
                ],
                style={
                    "color": "#b45309",
                    "fontWeight": "600",
                    "backgroundColor": "#fef3c7",
                    "padding": "6px 8px",
                    "borderRadius": "8px",
                },
            )
        else:
            alias_note = None

        return alias_note

    #============ methods - expand-contract ==================
    def is_node(self, elem):
        return "source" not in elem.get("data", {}) and "target" not in elem.get("data", {})


    def is_edge(self, elem):
        return "source" in elem.get("data", {}) and "target" in elem.get("data", {})


    def element_id(elem):
        """
        Cytoscape elements should preferably have data["id"].
        For edges, create a stable fallback id if missing.
        """
        data = elem.get("data", {})

        if "id" in data:
            return data["id"]

        if "source" in data and "target" in data:
            return f"{data['source']}__{data['target']}"

        return None


    def get_visible_node_ids(self, elements):
        return {
            elem["data"]["id"]
            for elem in elements
            if self.is_node(elem) and "id" in elem.get("data", {})
        }


    def get_visible_edge_ids(self, elements):
        return {
            self.element_id(elem)
            for elem in elements
            if self.is_edge(elem)
        }


    def get_neighbor_node_ids(self, node_id, all_elements):
        neighbors = set()

        for elem in all_elements:
            if not self.is_edge(elem):
                continue

            src = elem["data"]["source"]
            tgt = elem["data"]["target"]

            if src == node_id:
                neighbors.add(tgt)
            elif tgt == node_id:
                neighbors.add(src)

        return neighbors


    def get_edges_touching_visible_nodes(self, all_elements, visible_node_ids):
        """
        Return edges where both source and target are currently visible.
        """
        edges = []

        for elem in all_elements:
            if not self.is_edge(elem):
                continue

            src = elem["data"]["source"]
            tgt = elem["data"]["target"]

            if src in visible_node_ids and tgt in visible_node_ids:
                edges.append(elem)

        return edges

    def apply_positions_to_elements(self, elements, positions):
        new_elements = []

        for elem in elements:
            data = elem.get("data", {})
            node_id = data.get("id")

            if node_id in positions and "source" not in data and "target" not in data:
                elem = dict(elem)
                elem["position"] = positions[node_id]

            new_elements.append(elem)

        return new_elements


    def toggle_expand_contract(self,
        selected_node,
        current_elements,
        expanded_nodes,
        all_elements,
    ):
        if selected_node is None:
            return current_elements, expanded_nodes

        node_id = selected_node["id"]

        if expanded_nodes is None:
            expanded_nodes = []

        expanded_nodes = set(expanded_nodes)

        all_nodes_by_id = {
            elem["data"]["id"]: elem
            for elem in all_elements
            if self.is_node(elem)
        }

        currently_visible_node_ids = self.get_visible_node_ids(current_elements)

        # --------------------------------------------------
        # CASE 1: node is already expanded -> CONTRACT
        # --------------------------------------------------
        if node_id in expanded_nodes:
            expanded_nodes.remove(node_id)

            # Start from root/current important nodes
            # Keep nodes that are expanded or connected to expanded nodes
            nodes_to_keep = set()

            # Always keep the clicked node itself
            nodes_to_keep.add(node_id)

            # Keep all nodes currently expanded
            nodes_to_keep.update(expanded_nodes)

            # Keep neighbors of all still-expanded nodes
            for expanded_node_id in expanded_nodes:
                nodes_to_keep.add(expanded_node_id)
                nodes_to_keep.update(
                    self.get_neighbor_node_ids(expanded_node_id, all_elements)
                )

            # Optional but useful:
            # keep nodes that were initially visible before expansion
            for elem in current_elements:
                data = elem.get("data", {})
                if self.is_node(elem) and data.get("root", False):
                    nodes_to_keep.add(data["id"])

            visible_node_ids = currently_visible_node_ids.intersection(nodes_to_keep)

        # --------------------------------------------------
        # CASE 2: node is not expanded -> EXPAND
        # --------------------------------------------------
        else:
            expanded_nodes.add(node_id)

            neighbor_ids = self.get_neighbor_node_ids(node_id, all_elements)

            visible_node_ids = set(currently_visible_node_ids)
            visible_node_ids.add(node_id)
            visible_node_ids.update(neighbor_ids)

        # Rebuild visible nodes
        new_nodes = [
            all_nodes_by_id[nid]
            for nid in visible_node_ids
            if nid in all_nodes_by_id
        ]

        # Rebuild visible edges
        new_edges = self.get_edges_touching_visible_nodes(
            all_elements=all_elements,
            visible_node_ids=visible_node_ids,
        )

        new_elements = new_nodes + new_edges

        return new_elements, sorted(expanded_nodes)

    def find_rdf_obj_from_node_id(self, node_id):
        for node in self.rdf.subjects():
            if self.short(node) == node_id:
                return node

        return None

    def extract_node_info(self, node_data: dict, ndigits: int = 4) -> dict:
        """
        Extract standardized node information from Cytoscape node data
        and enrich it with HGNC/HUGO + RefSeq + UniProt annotations.
        """

        def first_available(*keys, default="NA"):
            for key in keys:
                value = node_data.get(key)
                if value not in [None, "", "NA"]:
                    return value
            return default

        def fmt_number(x, ndigits=4):
            if x in [None, "NA", ""]:
                return "NA"
            try:
                return f"{float(x):.{ndigits}g}"
            except Exception:
                return str(x)

        matches = []
        n_matches = 0

        annot = self.get_gene_annotation_for_node(node_data)

        node_id = first_available("id")
        node = self.find_rdf_obj_from_node_id(node_id)
        label = "??" if node is None else self.get_name(node)

        biopax_type = first_available(
            "biopax_type",
            "type",
            "biopax_class",
        )

        # Prefer annotation table values when available
        ensembl_id = (
            annot.get("ensembl_id")
            or first_available("ensembl_id", "gene_id", "gene", default="NA")
        )

        symbol = (
            annot.get("symbol")
            or first_available("symbol", "gene_symbol", default="NA")
        )

        name = (
            annot.get("name")
            or first_available("name", default="NA")
        )

        uniprot_id = (
            annot.get("uniprot_id")
            or first_available("uniprot_id", "uniprot", "UniProt", default="NA")
        )

        synonyms = self.get_aliases_for_symbol(symbol)

        refseq_summary = annot.get("refseq_summary", "NA")

        gene_type = first_available(
            "gene_type",
            "biotype",
        )

        lfc = first_available(
            "log2FoldChange",
            "log2FC",
            "LFC",
            "lfc",
        )

        fdr = first_available(
            "padj",
            "FDR",
            "fdr",
            "qvalue",
        )

        matches = annot.get("matches", [])
        n_matches = annot.get("n_matches", 0)

        # print("ANNOT:", annot)
        # print("matches:", matches)
        # print("n_matches:", n_matches)

        return {
            "id": node_id,
            "label": label,
            "biopax_type": biopax_type,
            "ensembl_id": ensembl_id,
            "symbol": symbol,
            "name": name,
            "gene_type": gene_type,
            "uniprot_id": uniprot_id,
            "synonyms": synonyms,
            "refseq_summary": refseq_summary,
            "lfc": fmt_number(lfc, ndigits),
            "fdr": fmt_number(fdr, ndigits),

            # keep alias/synonym ambiguity information
            "matches": matches,
            "n_matches": n_matches,
        }
    
    def make_node_info_panel(self, node_data: dict):
        info_data = self.extract_node_info(node_data)

        alias_note = self.none_one_multiple_alias(info_data)
        n_matches = info_data.get("n_matches", 0)

        if n_matches < 2:
            if str(info_data["ensembl_id"]) != 'NA':
                middle_gene_info = html.Div([
                    html.P([html.B("BioPAX type: "), str(info_data["biopax_type"])]),
                    html.P([html.B("Ensembl ID: "), str(info_data["ensembl_id"])]),
                    html.P([html.B("Symbol: "), str(info_data["symbol"])]),
                    html.P([html.B("Name: "), str(info_data["name"])]),
                    html.P([html.B("Gene type: "), str(info_data["gene_type"])]),
                    html.P([html.B("UniProt ID: "), str(info_data["uniprot_id"])]),
                    html.Hr(),
                    html.P([html.B("Synonyms: "), str(info_data["synonyms"])]),
                    html.Hr(),
                    html.P([html.B("LFC: "), str(info_data["lfc"])]),
                    html.P([html.B("FDR: "), str(info_data["fdr"])]),

                    html.Hr(),
                    html.P([html.B("RefSeq summary:")]), 

                    html.Div(
                        str(info_data["refseq_summary"]),
                        style={
                            "maxHeight": "220px",
                            "overflowY": "auto",
                            "fontSize": "13px",
                            "lineHeight": "1.45",
                            "backgroundColor": "#f9fafb",
                            "padding": "8px",
                            "borderRadius": "8px",
                            "border": "1px solid #e5e7eb",
                        },
                    ),                    
                ])
            else:
                middle_gene_info = html.Div([
                    html.P([html.B("BioPAX type: "), str(info_data["biopax_type"])]),
                ])                
        else:
            middle_gene_info = None

        return html.Div(
            [
                html.P([html.B("ID: "), str(info_data["id"])]),
                html.P([html.B("Label: "), str(info_data["label"])]),

                html.Hr(),
                alias_note,
                middle_gene_info,
            ]
        )

    def load_cyto_settings(self, default_font_size: int = 10) -> dict:
        if not self.settings_file.exists():
            return {"font_size": default_font_size}

        try:
            with open(self.settings_file, "r", encoding="utf-8") as f:
                settings = json.load(f)

            font_size = int(settings.get("font_size", default_font_size))
            font_size = max(6, min(font_size, 30))

            return {"font_size": font_size}

        except Exception as e:
            print(f"Could not load Cytoscape settings: {e}")
            return {"font_size": default_font_size}


    def save_cyto_settings(self, font_size: int) -> None:
        font_size = int(font_size)
        font_size = max(6, min(font_size, 30))

        settings = {
            "font_size": font_size,
        }

        with open(self.settings_file, "w", encoding="utf-8") as f:
            json.dump(settings, f, indent=2)    

    def get_main_hub_node_ids_cached(self, all_elements, top_n=10):
        """
        Return cached main hubs if already calculated.
        Main hubs = nodes with highest total degree.
        """

        if self.hub_list is not None and self.hub_top_n == top_n:
            return set(self.hub_list)

        node_ids = set()
        degree = Counter()

        for elem in all_elements:
            data = elem.get("data", {})

            is_node = (
                "id" in data
                and "source" not in data
                and "target" not in data
            )

            is_edge = (
                "source" in data
                and "target" in data
            )

            if is_node:
                node_ids.add(data["id"])

            elif is_edge:
                src = data["source"]
                tgt = data["target"]

                degree[src] += 1
                degree[tgt] += 1

        hubs = [
            node_id
            for node_id, _ in degree.most_common(top_n)
            if node_id in node_ids
        ]

        self.hub_list = hubs
        self.hub_top_n = top_n

        print(f">>>> calculated hubs once: {self.hub_list}")

        return set(self.hub_list)


    def get_source_sink_node_ids_cached(self, all_elements):
        """
        Source nodes = nodes with outgoing edges but no incoming edges.
        Sink nodes   = nodes with incoming edges but no outgoing edges.
        """

        if self.source_node_list is not None and self.sink_node_list is not None:
            return set(self.source_node_list), set(self.sink_node_list)

        node_ids = set()
        incoming = set()
        outgoing = set()

        for elem in all_elements:
            data = elem.get("data", {})

            is_node = (
                "id" in data
                and "source" not in data
                and "target" not in data
            )

            is_edge = (
                "source" in data
                and "target" in data
            )

            if is_node:
                node_ids.add(data["id"])

            elif is_edge:
                src = data["source"]
                tgt = data["target"]

                outgoing.add(src)
                incoming.add(tgt)

        source_nodes = sorted((node_ids & outgoing) - incoming)
        sink_nodes = sorted((node_ids & incoming) - outgoing)

        self.source_node_list = source_nodes
        self.sink_node_list = sink_nodes

        print(f">>>> calculated source nodes once: {len(source_nodes)}")
        print(f">>>> calculated sink nodes once: {len(sink_nodes)}")

        return set(source_nodes), set(sink_nodes)


    def create_cytoscape_app(self, height: str = "95%", width: str = "100%", marginTop: str = "20px", port: int = 8050):

        self.G = self.G.subgraph([n for n in self.G.nodes() if self.G.degree(n) > 0]).copy()

        elements = self.nx_to_cytoscape_elements()
        self.elements = elements

        app = dash.Dash(__name__, 
                        assets_folder=str(self.root_styles),
                        external_stylesheets=[dbc.themes.BOOTSTRAP],
                        suppress_callback_exceptions=True,
                        )

        title = f"{self.pathway} - {self.pathway_id}"

        print(">>> create_cytoscape_app() ", title)

        """
        save_and_notify          -> updates toast + saved-output
        show_node_info           -> updates node-info + selected-node-store
        update_graph_visibility  -> updates elements + expanded-nodes-store
        update_layout            -> updates layout
        show_selected_nodes      -> updates selected-output
        """

        base_font_size = 10
        cyto_settings = self.load_cyto_settings(default_font_size=base_font_size)
        initial_font_size = cyto_settings["font_size"]

        def make_stylesheet(font_size=10):
            edge_font_size = max(font_size - 3, 6)
            return [
                {
                    "selector": "node",
                    "style": {
                        "label": "data(label)",
                        "font-size": f"{font_size}px",
                        "text-valign": "center",
                        "text-halign": "center",
                        "background-color": "#8ecae6",
                        "width": 30,
                        "height": 30,
                    },
                },

                {
                    "selector": '[lfc_bin = "up_1"]',
                    "style": {
                        "background-color": "#fcae91",
                        "border-color": "#cb181d",
                        "border-width": 2,
                        "shape": "ellipse",
                    },
                },
                {
                    "selector": '[lfc_bin = "up_2"]',
                    "style": {
                        "background-color": "#fb6a4a",
                        "border-color": "#a50f15",
                        "border-width": 3,
                        "shape": "ellipse",
                    },
                },
                {
                    "selector": '[lfc_bin = "up_3"]',
                    "style": {
                        "background-color": "#cb181d",
                        "border-color": "#67000d",
                        "border-width": 4,
                        "shape": "ellipse",
                    },
                },
                {
                    "selector": '[lfc_bin = "down_1"]',
                    "style": {
                        "background-color": "#bdd7e7",
                        "border-color": "#2171b5",
                        "border-width": 2,
                        "shape": "ellipse",
                    },
                },
                {
                    "selector": '[lfc_bin = "down_2"]',
                    "style": {
                        "background-color": "#6baed6",
                        "border-color": "#08519c",
                        "border-width": 3,
                        "shape": "ellipse",
                    },
                },
                {
                    "selector": '[lfc_bin = "down_3"]',
                    "style": {
                        "background-color": "#2171b5",
                        "border-color": "#08306b",
                        "border-width": 4,
                        "shape": "ellipse",
                    },
                },
                {
                    "selector": '[lfc_bin = "weak"]',
                    "style": {
                        "background-color": "#eeeeee",
                        "border-color": "#999999",
                        "border-width": 1,
                    },
                },
                {
                    "selector": '[lfc_bin = "none"]',
                    "style": {
                        "background-color": "#dddddd",
                        "border-color": "#aaaaaa",
                        "border-width": 1,
                    },
                },


                {
                    "selector": "node:selected",
                    "style": {
                        "border-width": 4,
                        "border-color": "black",
                        "background-color": "#ffcc00",
                    },
                },
                {
                    "selector": "edge",
                    "style": {
                        "curve-style": "bezier",
                        "target-arrow-shape": "triangle",
                        "label": "data(interaction)",
                        "font-size": f"{edge_font_size}px",
                    },
                },
                {
                    "selector": '[biopax_type = "BiochemicalReaction"]',
                    "style": {
                        "background-color": "#f94144",
                        "shape": "hexagon",
                    },
                },
                {
                    "selector": '[biopax_type = "Protein"]',
                    "style": {
                        "background-color": "#8ecae6",
                        "shape": "ellipse",
                    },
                },
                {
                    "selector": '[biopax_type = "Complex"]',
                    "style": {
                        "background-color": "#ffb703",
                        "shape": "round-rectangle",
                    },
                },
                {
                    "selector": '[biopax_type = "SmallMolecule"]',
                    "style": {
                        "background-color": "#90be6d",
                        "shape": "diamond",
                    },
                },
                {
                    "selector": '[biopax_type = "Pathway"]',
                    "style": {
                        "background-color": "#cdb4db",
                        "shape": "rectangle",
                    },
                },
            ]

        app.layout = html.Div(
            [
                html.H2(
                    title,
                    style={
                        "textAlign": "center",
                        "color": "#333",
                        "marginTop": "5px",
                        "marginBottom": "5px",
                    },
                ),

                # Main horizontal layout: dropdown / show/hide button
                html.Div(
                    [
                        html.Div(
                            dcc.Dropdown(
                                id="layout-dropdown",
                                options=[
                                    {"label": "COSE force-directed", "value": "cose"},
                                    {"label": "Breadthfirst / hierarchical", "value": "breadthfirst"},
                                    {"label": "Circle", "value": "circle"},
                                    {"label": "Grid", "value": "grid"},
                                    {"label": "Preset", "value": "preset"},
                                ],
                                value="preset",
                                clearable=False,
                                style={
                                    "width": "350px",
                                },
                            ),
                            style={
                                "minWidth": "0",
                            },
                        ),

                        html.Div(
                            html.Button(
                                "Hide node info",
                                id="toggle-node-panel-button",
                                n_clicks=0,
                                className="cyto-button cyto-button-hide",
                                style={
                                    "width": "100%",
                                    "fontSize": "12px",
                                    "padding": "8px 12px",
                                    "backgroundColor": "#ddc3a5",
                                    "border": "1px solid #ddc3a5",
                                },
                            ),
                            style={
                                "width": "100%",
                            },
                        ),
                    ],
                    id="top-cyto-controls",
                    style={
                        "display": "grid",
                        "gridTemplateColumns": "5fr 1fr",
                        "gap": "12px",
                        "alignItems": "center",
                        "width": "100%",
                        "marginBottom": "10px",
                    },
                ),
                
                # Graph + Sidebar
                html.Div(
                    [
                        html.Div(
                            [
                                cyto.Cytoscape(
                                            id="reactome-network",
                                            elements=elements,
                                            layout={"name": "preset"},
                                            boxSelectionEnabled=True,
                                            autoungrabify=False,
                                            autounselectify=False,
                                            stylesheet=make_stylesheet(initial_font_size),
                                            zoom=1.0,
                                            minZoom=0.2,
                                            maxZoom=3.0,
                                            style={
                                                "width": "100%",
                                                "height": "80vh",
                                                "backgroundColor": "#fff9c4",
                                                "marginTop": marginTop,
                                                "border": "1px solid #ddd",
                                                "borderRadius": "12px",
                                            },
                                        ),

                                html.Div(
                                    id="cyto-popup-menu",
                                    children=[
                                        html.Button("1️⃣ Select 1st neighbors", id="btn-select-neighbors", className="popup-button"),
                                        html.Button("⬆️ Select upstream", id="btn-select-upstream", className="popup-button"),
                                        html.Button("⬇️ Select downstream", id="btn-select-downstream", className="popup-button"),
                                        html.Button("⭐ Select main hubs", id="btn-select-hubs", className="popup-button"),
                                        html.Button("🔼 Select most upstream nodes", id="btn-select-sources", className="popup-button"),
                                        html.Button("🔽 Select most downstream nodes", id="btn-select-sinks", className="popup-button"),
                                        
                                        html.Hr(style={"margin": "8px 0"}),

                                        html.Button("🎯 Open Targets: pathways", id="btn-ot-pathways", className="popup-button popup-button-ot", style={"display": "none"},),
                                        html.Button("🧬 Open Targets: diseases", id="btn-ot-diseases", className="popup-button popup-button-ot", style={"display": "none"},),
                                        html.Button("💊 Open Targets: drugs for TCGA disease", id="btn-ot-drugs", className="popup-button popup-button-ot",style={"display": "none"},),

                                        html.Hr(style={"margin": "8px 0"}),

                                        html.Button("✖ Close", id="btn-close-popup", className="popup-button popup-close"),
                                    ],
                                    style={
                                        "display": "none",
                                        "position": "fixed",
                                        "zIndex": 9999,
                                        "backgroundColor": "white",
                                        "border": "1px solid #ddd",
                                        "borderRadius": "12px",
                                        "boxShadow": "0 8px 24px rgba(0,0,0,0.18)",
                                        "padding": "10px",
                                        "width": "230px",
                                    },
                                ),
                            ],
                            id="graph-panel",
                             style={
                                "position": "relative",
                                "minWidth": "0",
                            },
                        ),
   
                        # Sidebar
                        html.Div(
                            [
                                html.H4(
                                    "Node information",
                                    className="node-panel-title",
                                    style={"margin": "0"},
                                ),

                                html.Div(id="node-info", className="node-info-box"),
                                dcc.Store(id="cyto-font-size-store", data=initial_font_size),
                                dcc.Store(id="cyto-zoom-store", data=1.0),

                                # Search for genes
                                html.Div(
                                    html.Div(
                                        [
                                            html.Span(
                                                "Search for a gene:",
                                                style={
                                                    "fontWeight": "600",
                                                    "fontSize": "13px",
                                                    "marginRight": "8px",
                                                    "display": "block",
                                                },
                                            ),


                                            html.Div(
                                                [
                                                    dcc.Input(
                                                        id="gene-search-input",
                                                        type="text",
                                                        placeholder="e.g. TP, MAPK, CCN...",
                                                        debounce=True,
                                                        style={
                                                            "width": "80px",
                                                            "padding": "6px 8px",
                                                            "border": "1px solid #d0d7de",
                                                            "borderRadius": "8px",
                                                            "fontSize": "13px",
                                                        },
                                                    ),

                                                    html.Button(
                                                        "🔎",
                                                        id="gene-search-button",
                                                        n_clicks=0,
                                                        title="Find genes starting with this text",
                                                        className="cyto-button",
                                                        style={
                                                            "marginLeft": "6px",
                                                            "padding": "6px 10px",
                                                        },
                                                    ),
                                                ],
                                                style={
                                                    "display": "flex",
                                                    "alignItems": "center",
                                                    "gap": "4px",
                                                },
                                            ),
                                        ],
                                        style={
                                            "display": "flex",
                                            "flexDirection": "column",
                                            "alignItems": "flex-start",
                                            "gap": "4px",
                                        },
                                    ),
                                ),

                                # Search for pathways between genes
                                html.Div(
                                    [
                                        html.Span(
                                            "Find all paths between:",
                                            style={
                                                "fontWeight": "600",
                                                "fontSize": "13px",
                                                "marginBottom": "6px",
                                                "display": "block",
                                            },
                                        ),

                                        html.Div(
                                            [
                                                dcc.Input(
                                                    id="path-source-input",
                                                    type="text",
                                                    placeholder="source gene",
                                                    debounce=True,
                                                    style={
                                                        "width": "80px",
                                                        "padding": "6px 8px",
                                                        "border": "1px solid #d0d7de",
                                                        "borderRadius": "8px",
                                                        "fontSize": "13px",
                                                    },
                                                ),

                                                html.Span(
                                                    "and",
                                                    style={
                                                        "fontSize": "13px",
                                                        "margin": "0 4px",
                                                    },
                                                ),

                                                dcc.Input(
                                                    id="path-target-input",
                                                    type="text",
                                                    placeholder="target gene",
                                                    debounce=True,
                                                    style={
                                                        "width": "80px",
                                                        "padding": "6px 8px",
                                                        "border": "1px solid #d0d7de",
                                                        "borderRadius": "8px",
                                                        "fontSize": "13px",
                                                    },
                                                ),

                                                html.Button(
                                                    "🔎",
                                                    id="find-paths-button",
                                                    n_clicks=0,
                                                    title="Find paths between both genes",
                                                    className="cyto-button",
                                                    style={
                                                        "marginLeft": "6px",
                                                        "padding": "6px 10px",
                                                    },
                                                ),
                                            ],
                                            style={
                                                "display": "flex",
                                                "alignItems": "center",
                                                "gap": "4px",
                                            },
                                        ),
                                    ],
                                    style={
                                        "display": "flex",
                                        "flexDirection": "column",
                                        "alignItems": "flex-start",
                                        "gap": "4px",
                                    },
                                ),

                                html.Div(
                                    html.Button(
                                        "🧬 Show DEGs",
                                        id="show-degs-button",
                                        n_clicks=0,
                                        className="cyto-button",
                                        style={
                                            "width": "100%",
                                            "marginTop": "8px",
                                            "padding": "8px 12px",
                                            "fontWeight": "700",
                                            "backgroundColor": "#7bd398",
                                            "border": "1px solid #7bd398",
                                        },
                                    ),
                                ),

                                # font controls
                                html.Div(
                                    [
                                        html.Span(
                                            "Font size",
                                            style={
                                                "fontWeight": "600",
                                                "fontSize": "13px",
                                                "marginRight": "8px",
                                                "whiteSpace": "nowrap",
                                            },
                                        ),

                                        html.Div(
                                            [
                                                html.Button(
                                                    "A−",
                                                    id="decrease-font-btn",
                                                    n_clicks=0,
                                                    title="Decrease font size",
                                                    className="cyto-font-button",
                                                    style={
                                                        "borderTopRightRadius": "0px",
                                                        "borderBottomRightRadius": "0px",
                                                        "borderRight": "0px",
                                                    },
                                                ),
                                                html.Button(
                                                    "A+",
                                                    id="increase-font-btn",
                                                    n_clicks=0,
                                                    title="Increase font size",
                                                    className="cyto-font-button",
                                                    style={
                                                        "borderTopLeftRadius": "0px",
                                                        "borderBottomLeftRadius": "0px",
                                                    },
                                                ),
                                            ],
                                            style={
                                                "display": "flex",
                                                "flexDirection": "row",
                                                "alignItems": "center",
                                            },
                                        ),

                                        html.Span(
                                            id="font-size-label",
                                            children=f"{initial_font_size}px",
                                            style={
                                                "fontSize": "12px",
                                                "color": "#555",
                                                "marginLeft": "8px",
                                                "minWidth": "36px",
                                            },
                                        ),
                                    ],
                                    style={
                                        "display": "flex",
                                        "flexDirection": "row",
                                        "alignItems": "center",
                                        "justifyContent": "flex-start",
                                        "gap": "0px",
                                        "marginTop": "10px",
                                        "marginBottom": "12px",
                                        "padding": "8px",
                                        "border": "1px solid #ddd",
                                        "borderRadius": "10px",
                                        "backgroundColor": "#fafafa",
                                    },
                                ),

                                # zoom controls
                                html.Div(
                                    [
                                        html.Span(
                                            "Zoom",
                                            style={
                                                "fontWeight": "600",
                                                "fontSize": "13px",
                                                "marginRight": "8px",
                                                "whiteSpace": "nowrap",
                                            },
                                        ),

                                        dcc.Slider(
                                            id="cyto-zoom-slider",
                                            min=0.2,
                                            max=3.0,
                                            step=0.1,
                                            value=1.0,
                                            marks={
                                                0.5: "0.5x",
                                                1.0: "1x",
                                                2.0: "2x",
                                                3.0: "3x",
                                            },
                                            tooltip={
                                                "placement": "bottom",
                                                "always_visible": False,
                                            },
                                        ),
                                    ],
                                    style={
                                        "marginTop": "10px",
                                        "marginBottom": "12px",
                                        "padding": "8px",
                                        "border": "1px solid #ddd",
                                        "borderRadius": "10px",
                                        "backgroundColor": "#fafafa",
                                    },
                                ),

                                html.Button(
                                    "💾 Save node positions",
                                    id="save-graph-button",
                                    n_clicks=0,
                                    className="cyto-button cyto-button-save",
                                ),

                                # reset/save buttons
                                html.Button(
                                    "🔄 Reset graph",
                                    id="reset-graph-button",
                                    n_clicks=0,
                                    className="cyto-button cyto-button-reset",
                                ),

                            ],
                            id="node-sidebar",
                            className="node-sidebar",
                        ),
                    ],
                    id="main-cyto-layout",
                    style={
                        "width": "100%",
                        "display": "grid",
                        "gridTemplateColumns": "5fr 1fr",
                        "gap": "12px",
                        "alignItems": "start",
                    },
                ),

                dcc.Store(id="selected-node-store"),
                dcc.Store(id="hover-node-store"),
                dcc.Store(id="right-click-event-store"),
                dcc.Store(id="right-click-node-store"),
                dcc.Store(id="expanded-nodes-store", data=[]),
                dcc.Store(id="all-elements-store", data=elements),

                dcc.Store(id="opentarget-action-store"),

                html.Pre(id="saved-output"),
                html.Pre(id="selected-output"),

                dbc.Alert(
                    id="path-alert",
                    children="",
                    color="warning",
                    is_open=False,
                    dismissable=True,
                    style={
                        "marginTop": "8px",
                        "fontSize": "13px",
                        "padding": "8px 10px",
                        "borderRadius": "10px",
                    },
                ),

                dbc.Toast(
                    id="save-toast",
                    header="Success",
                    is_open=False,
                    dismissable=True,
                    duration=2500,
                    icon="success",
                    style={
                        "position": "fixed",
                        "top": 20,
                        "right": 20,
                        "width": 360,
                        "zIndex": 9999,
                    },
                ),

                dbc.Modal(
                    [
                        dbc.ModalHeader( dbc.ModalTitle("DEGs in Reactome graph") ),
                        dbc.ModalBody( html.Div(id="deg-summary-window") ),

                        dbc.ModalFooter(
                            dbc.Button(
                                "Close",
                                id="close-deg-summary-button",
                                className="ms-auto",
                                n_clicks=0,
                            )
                        ),
                    ],
                    id="deg-summary-modal",
                    is_open=False,
                    size="xl",
                    scrollable=True,
                ),

                dbc.Modal(
                    [
                        dbc.ModalHeader(
                            dbc.ModalTitle(id="opentarget-modal-title")
                        ),

                        dbc.ModalBody(
                            [
                                html.Div(
                                    id="opentarget-modal-message",
                                    style={
                                        "marginBottom": "12px",
                                        "fontSize": "13px",
                                        "color": "#444",
                                    },
                                ),

                                html.Div(id="opentarget-modal-table"),
                            ]
                        ),

                        dbc.ModalFooter(
                            dbc.Button(
                                "Close",
                                id="btn-close-opentarget-modal",
                                className="ms-auto",
                                n_clicks=0,
                            )
                        ),
                    ],
                    id="opentarget-modal",
                    is_open=False,
                    size="xl",
                    scrollable=True,
                ),                

            ],
        )

        def classify_degs_in_graph(all_elements):
            """
            Classify DEGs as:
            - found once in graph
            - found more than once in graph
            - not found in graph

            Uses self.dflfc_ori symbols as DEG reference.
            Uses Cytoscape node labels as graph symbols.
            """

            if self.dflfc_ori is None or self.dflfc_ori.empty:
                return 0, [], [], [], set()
            
            dflfc = self.dflfc_ori[
                (self.dflfc_ori["abs_lfc"] >= self.lfc_cutoff) &
                (self.dflfc_ori["fdr"] < self.fdr_cutoff)
            ].copy()

            if dflfc.empty:
                return 0, [], [], [], set()

            dflfc.reset_index(drop=True, inplace=True)

            # Reference DEG symbols from the DEG table
            deg_symbols = np.unique(dflfc["symbol"].to_list())

            # Count graph node labels
            graph_symbol_counts = Counter()
            deg_node_ids = set()

            for elem in all_elements:
                data = elem.get("data", {})

                is_node = (
                    "id" in data
                    and "source" not in data
                    and "target" not in data
                )

                if not is_node:
                    continue

                node_id = data.get("id")
                label = str(data.get("label", node_id)).strip().upper()

                if not label:
                    continue

                graph_symbol_counts[label] += 1

            found_once = []
            found_multiple = []

            for symbol in deg_symbols:
                n = graph_symbol_counts.get(symbol, 0)

                if n == 1:
                    found_once.append(symbol)
                elif n > 1:
                    found_multiple.append(f"{symbol} ({n}x)")


            not_found = []
            # self.pathway_genes:
            for enr_symbol in self.found_degs:
                if enr_symbol not in deg_symbols:
                    not_found.append(enr_symbol)

            # Select all graph nodes whose labels are DEG symbols
            deg_symbol_set = set(deg_symbols)

            for elem in all_elements:
                data = elem.get("data", {})

                is_node = (
                    "id" in data
                    and "source" not in data
                    and "target" not in data
                )

                if not is_node:
                    continue

                node_id = data.get("id")
                label = str(data.get("label", node_id)).strip().upper()

                if label in deg_symbol_set:
                    deg_node_ids.add(node_id)

            return len(deg_symbols), found_once, found_multiple, not_found, deg_node_ids


        def make_deg_list_panel(title, items, color, background):
            if not items:
                children = html.Div(
                    "None",
                    style={
                        "fontSize": "13px",
                        "color": "#666",
                        "fontStyle": "italic",
                    },
                )
            else:
                children = html.Div(
                    ", ".join(items),
                    style={
                        "fontSize": "13px",
                        "lineHeight": "1.6",
                        "maxHeight": "220px",
                        "overflowY": "auto",
                        "whiteSpace": "normal",
                    },
                )

            return html.Div(
                [
                    html.H5(
                        f"{title} ({len(items)})",
                        style={
                            "color": color,
                            "fontWeight": "800",
                            "marginBottom": "8px",
                        },
                    ),
                    children,
                ],
                style={
                    "backgroundColor": background,
                    "border": f"1px solid {color}",
                    "borderRadius": "12px",
                    "padding": "12px",
                    "marginBottom": "12px",
                },
            )

        def select_nodes_in_elements(elements, selected_ids):
            """
            Deselect all elements, then select only nodes in selected_ids.
            """
            selected_ids = set(selected_ids or [])
            new_elements = []

            for elem in elements:
                elem2 = dict(elem)
                data = elem2.get("data", {})

                is_node = (
                    "id" in data
                    and "source" not in data
                    and "target" not in data
                )

                elem2["selected"] = False

                if is_node and data.get("id") in selected_ids:
                    elem2["selected"] = True

                new_elements.append(elem2)

            return new_elements


        def force_select_single_node(elements, node_id):
            """
            Deselect all elements and select only node_id.
            """
            if not elements or not node_id:
                return dash.no_update

            new_elements = []

            for elem in elements:
                elem2 = dict(elem)
                data = elem2.get("data", {})

                is_node = (
                    "id" in data
                    and "source" not in data
                    and "target" not in data
                )

                # remove old selection from everything
                elem2["selected"] = False

                # select only the right-clicked node
                if is_node and data.get("id") == node_id:
                    elem2["selected"] = True

                new_elements.append(elem2)

            return new_elements


        def get_upstream_node_ids(node_id, all_elements):
            """
            Direct upstream = nodes with edges source -> node_id.
            """
            upstream = set()

            for elem in all_elements:
                data = elem.get("data", {})

                if "source" not in data or "target" not in data:
                    continue

                if data["target"] == node_id:
                    upstream.add(data["source"])

            return upstream


        def get_downstream_node_ids(node_id, all_elements):
            """
            Direct downstream = nodes with edges node_id -> target.
            """
            downstream = set()

            for elem in all_elements:
                data = elem.get("data", {})

                if "source" not in data or "target" not in data:
                    continue

                if data["source"] == node_id:
                    downstream.add(data["target"])

            return downstream


        def normalize_gene_text(x):
            if x is None:
                return ""
            return str(x).strip().upper()


        def get_node_label_from_elements(node_id, elements):
            for elem in elements:
                data = elem.get("data", {})

                is_node = (
                    "id" in data
                    and "source" not in data
                    and "target" not in data
                )

                if is_node and data.get("id") == node_id:
                    return data.get("label", node_id)

            return node_id


        def find_nodes_starting_with(prefix, all_elements):
            """
            Find nodes whose label or id starts with prefix.
            """
            prefix = normalize_gene_text(prefix)

            if not prefix:
                return set()

            hits = set()

            for elem in all_elements:
                data = elem.get("data", {})

                is_node = (
                    "id" in data
                    and "source" not in data
                    and "target" not in data
                )

                if not is_node:
                    continue

                node_id = str(data.get("id", ""))
                label = str(data.get("label", ""))

                if (
                    normalize_gene_text(label).startswith(prefix)
                    or normalize_gene_text(node_id).startswith(prefix)
                ):
                    hits.add(node_id)

            return hits


        def resolve_gene_query_to_node_ids(query, all_elements):
            """
            Resolve user text to matching node ids.
            First tries exact label/id match.
            Then tries startswith.
            """
            query_norm = normalize_gene_text(query)

            if not query_norm:
                return []

            exact_hits = []
            prefix_hits = []

            for elem in all_elements:
                data = elem.get("data", {})

                is_node = (
                    "id" in data
                    and "source" not in data
                    and "target" not in data
                )

                if not is_node:
                    continue

                node_id = str(data.get("id", ""))
                label = str(data.get("label", ""))

                node_id_norm = normalize_gene_text(node_id)
                label_norm = normalize_gene_text(label)

                if query_norm in {node_id_norm, label_norm}:
                    exact_hits.append(node_id)

                elif label_norm.startswith(query_norm) or node_id_norm.startswith(query_norm):
                    prefix_hits.append(node_id)

            if exact_hits:
                return exact_hits

            return prefix_hits


        def select_nodes_and_edges_in_elements(elements: list, selected_node_ids: set=set(), selected_edge_pairs: set=set()) -> list:
            """
            Select nodes and optionally edges.
            selected_edge_pairs should contain tuples: (source, target)
            """
            selected_node_ids = set(selected_node_ids)
            selected_edge_pairs = set(selected_edge_pairs)

            new_elements = []

            for elem in elements:
                elem2 = dict(elem)
                data = elem2.get("data", {})

                is_node = (
                    "id" in data
                    and "source" not in data
                    and "target" not in data
                )

                is_edge = (
                    "source" in data
                    and "target" in data
                )

                elem2["selected"] = False

                if is_node and data.get("id") in selected_node_ids:
                    elem2["selected"] = True

                elif is_edge:
                    edge_pair = (data.get("source"), data.get("target"))

                    if edge_pair in selected_edge_pairs:
                        elem2["selected"] = True

                new_elements.append(elem2)

            return new_elements

        def find_paths_between_gene_queries(
            source_query,
            target_query,
            all_elements,
            cutoff=8,
            max_paths=100,
        ):
            source_nodes = resolve_gene_query_to_node_ids(source_query, all_elements)
            target_nodes = resolve_gene_query_to_node_ids(target_query, all_elements)

            if not source_nodes or not target_nodes:
                return set(), set(), source_nodes, target_nodes, 0

            selected_nodes = set()
            selected_edges = set()
            n_paths = 0

            for source_id in source_nodes:
                for target_id in target_nodes:
                    if source_id not in self.G.nodes or target_id not in self.G.nodes:
                        continue

                    try:
                        paths = nx.all_simple_paths(
                            self.G,
                            source=source_id,
                            target=target_id,
                            cutoff=cutoff,
                        )

                        for path in paths:
                            n_paths += 1

                            selected_nodes.update(path)

                            for a, b in zip(path[:-1], path[1:]):
                                selected_edges.add((a, b))

                            if n_paths >= max_paths:
                                return (
                                    selected_nodes,
                                    selected_edges,
                                    source_nodes,
                                    target_nodes,
                                    n_paths,
                                )

                    except nx.NetworkXNoPath:
                        continue

                    except nx.NodeNotFound:
                        continue

            return selected_nodes, selected_edges, source_nodes, target_nodes, n_paths


        def clean_and_select_elements(elements, selected_node_ids=None, selected_edge_pairs=None):
            """
            Remove all previous selection states, then select requested nodes/edges.
            Keeps positions and other element fields.
            """
            selected_node_ids = set(selected_node_ids or [])
            selected_edge_pairs = set(selected_edge_pairs or [])

            new_elements = []

            for elem in elements:
                elem2 = dict(elem)
                data = elem2.get("data", {})

                is_node = (
                    "id" in data
                    and "source" not in data
                    and "target" not in data
                )

                is_edge = (
                    "source" in data
                    and "target" in data
                )

                # very important: remove previous selection
                elem2["selected"] = False

                if is_node:
                    if data.get("id") in selected_node_ids:
                        elem2["selected"] = True

                elif is_edge:
                    pair = (data.get("source"), data.get("target"))
                    if pair in selected_edge_pairs:
                        elem2["selected"] = True

                new_elements.append(elem2)

            return new_elements

        def get_node_gene_symbol(node_data: dict) -> str:
            if not node_data:
                return ""

            return (
                node_data.get("symbol")
                or node_data.get("gene_symbol")
                or node_data.get("approvedSymbol")
                or node_data.get("label")
                or ""
            )

        def is_gene_node(node_data: dict) -> bool:
            if not node_data:
                return False

            symbol = str(
                node_data.get("symbol")
                or node_data.get("label")
                or ""
            ).strip()

            ensembl_id = str(
                node_data.get("ensembl_id")
                or ""
            ).strip()

            has_symbol = symbol not in ["", "-", "NA", "None"]
            has_ensembl = ensembl_id not in ["", "-", "NA", "None"] and ensembl_id.startswith("ENSG")

            # For Open Targets I recommend requiring Ensembl,
            # because symbol-only aliases can be ambiguous.
            return has_symbol and has_ensembl


        def get_gene_symbol_from_node(node_data: dict) -> str:
            if not is_gene_node(node_data):
                return ""

            symbol = node_data.get("symbol")

            if symbol in [None, "", "NA", "None"]:
                symbol = node_data.get("label", "")

            return str(symbol).strip()


        def df_to_dash_table(df, max_rows: int = 500):
            if df is None or df.empty:
                return html.Div(
                    "No records found.",
                    style={
                        "padding": "12px",
                        "backgroundColor": "#f8f9fa",
                        "borderRadius": "8px",
                        "color": "#555",
                    },
                )

            df = df.copy().head(max_rows)

            return dash_table.DataTable(
                data=df.to_dict("records"),
                columns=[{"name": str(c), "id": str(c)} for c in df.columns],
                page_size=15,
                sort_action="native",
                filter_action="native",
                style_table={
                    "overflowX": "auto",
                    "maxHeight": "650px",
                    "overflowY": "auto",
                    "border": "1px solid #eee",
                    "borderRadius": "8px",
                },
                style_cell={
                    "fontFamily": "Arial",
                    "fontSize": "12px",
                    "padding": "7px",
                    "textAlign": "left",
                    "maxWidth": "320px",
                    "whiteSpace": "normal",
                    "height": "auto",
                },
                style_header={
                    "fontWeight": "bold",
                    "backgroundColor": "#f6f8fa",
                    "borderBottom": "1px solid #ddd",
                },
            )

        @app.callback(
            Output("reactome-network", "elements", allow_duplicate=True),
            Output("deg-summary-window", "children"),
            Output("deg-summary-modal", "is_open"),

            Input("show-degs-button", "n_clicks"),
            Input("close-deg-summary-button", "n_clicks"),

            State("reactome-network", "elements"),
            State("all-elements-store", "data"),
            State("deg-summary-modal", "is_open"),

            prevent_initial_call=True,
        )
        def show_degs_in_graph(
            show_clicks,
            close_clicks,
            current_elements,
            all_elements,
            is_open,
        ):
            trigger = dash.callback_context.triggered[0]["prop_id"].split(".")[0]

            if trigger == "close-deg-summary-button":
                return dash.no_update, dash.no_update, False

            if trigger != "show-degs-button" or not show_clicks:
                raise dash.exceptions.PreventUpdate

            n_degs, found_once, found_multiple, not_found, deg_node_ids = classify_degs_in_graph(
                all_elements=all_elements,
            )

            new_elements = clean_and_select_elements(
                elements=current_elements,
                selected_node_ids=deg_node_ids,
                selected_edge_pairs=[],
            )

            summary_children = html.Div(
                [
                    html.P(
                        [
                            html.B(f"{n_degs} DEGs; {len(deg_node_ids)} DEGs found as nodes.")
                        ],
                        style={
                            "fontSize": "18px",
                            "marginBottom": "12px",
                        },
                    ),

                    make_deg_list_panel(
                        title="DEGs found once in graph",
                        items=found_once,
                        color="#166534",
                        background="#dcfce7",
                    ),

                    make_deg_list_panel(
                        title="DEGs found more than once in graph",
                        items=found_multiple,
                        color="#0f172a",
                        background="#dbeafe",
                    ),

                    make_deg_list_panel(
                        title="DEGs not found in graph",
                        items=not_found,
                        color="#b91c1c",
                        background="#fee2e2",
                    ),

                    make_deg_list_panel(
                        title="DEGs in the pathway",
                        items=self.pathway_genes,
                        color="#b91c1c",
                        background="#fee2e2",
                    ),
                ]
            )

            return new_elements, summary_children, True

        @app.callback(
            Output("reactome-network", "elements", allow_duplicate=True),
            Output("selected-output", "children", allow_duplicate=True),

            Input("gene-search-button", "n_clicks"),
            State("gene-search-input", "value"),
            State("reactome-network", "elements"),
            State("all-elements-store", "data"),

            prevent_initial_call=True,
        )
        def search_gene_by_prefix(
            n_clicks,
            query,
            current_elements,
            all_elements,
        ):
            if not n_clicks:
                raise dash.exceptions.PreventUpdate

            selected_ids = find_nodes_starting_with(
                prefix=query,
                all_elements=all_elements,
            )

            new_elements = clean_and_select_elements(
                elements=current_elements,
                selected_node_ids=selected_ids,
            )

            msg = f"Found {len(selected_ids)} node(s) starting with: {query}"

            return new_elements, msg

        @app.callback(
            Output("reactome-network", "elements", allow_duplicate=True),
            Output("selected-output", "children", allow_duplicate=True),
            Output("path-alert", "children"),
            Output("path-alert", "is_open"),

            Input("find-paths-button", "n_clicks"),

            State("path-source-input", "value"),
            State("path-target-input", "value"),
            State("reactome-network", "elements"),
            State("all-elements-store", "data"),

            prevent_initial_call=True,
        )
        def find_paths_between_genes(
            n_clicks,
            source_query,
            target_query,
            current_elements,
            all_elements,
        ):
            if not n_clicks:
                raise dash.exceptions.PreventUpdate

            if not source_query or not target_query:
                return (
                    dash.no_update,
                    "Please provide both source and target genes.",
                    "Please provide both source and target genes.",
                    True,
                )

            selected_nodes, selected_edges, source_hits, target_hits, n_paths = (
                find_paths_between_gene_queries(
                    source_query=source_query,
                    target_query=target_query,
                    all_elements=all_elements,
                    cutoff=8,
                    max_paths=100,
                )
            )

            # Case 1: source or target not found
            if not source_hits or not target_hits:
                msg = (
                    f"Could not find source or target node. "
                    f"Source hits: {len(source_hits)} | "
                    f"Target hits: {len(target_hits)}"
                )

                new_elements = clean_and_select_elements(
                    elements=current_elements,
                    selected_node_ids=[],
                    selected_edge_pairs=[],
                )

                return new_elements, msg, msg, True            

            # Case 2: nodes exist, but no path was found
            if n_paths == 0:
                msg = f"There is no pathway between {source_query} and {target_query}."

                new_elements = clean_and_select_elements(
                    elements=current_elements,
                    selected_node_ids=[],
                    selected_edge_pairs=[],
                )

                return new_elements, msg, msg, True
    
            # Case 3: paths found
            new_elements = clean_and_select_elements(
                elements=current_elements,
                selected_node_ids=selected_nodes,
                selected_edge_pairs=selected_edges,
            )

            msg = (
                f"Source hits: {len(source_hits)} | "
                f"Target hits: {len(target_hits)} | "
                f"Paths found: {n_paths} | "
                f"Selected nodes: {len(selected_nodes)} | "
                f"Selected edges: {len(selected_edges)}"
            )

            return new_elements, msg, "", False

        @app.callback(
            Output("reactome-network", "elements", allow_duplicate=True),
            Output("cyto-popup-menu", "style", allow_duplicate=True),

            Input("btn-select-upstream", "n_clicks"),

            State("right-click-node-store", "data"),
            State("reactome-network", "elements"),
            State("all-elements-store", "data"),
            State("cyto-popup-menu", "style"),

            prevent_initial_call=True,
        )
        def select_upstream_and_close_popup(
            n_clicks,
            right_click_node,
            current_elements,
            all_elements,
            current_style,
        ):
            if not n_clicks or not right_click_node:
                raise dash.exceptions.PreventUpdate

            node_id = right_click_node.get("id")

            if not node_id:
                raise dash.exceptions.PreventUpdate

            upstream_ids = get_upstream_node_ids(
                node_id=node_id,
                all_elements=all_elements,
            )

            selected_ids = set(upstream_ids)
            selected_ids.add(node_id)

            new_elements = clean_and_select_elements(
                elements=current_elements,
                selected_node_ids=selected_ids,
            )

            style = dict(current_style or {})
            style["display"] = "none"

            return new_elements, style


        @app.callback(
            Output("reactome-network", "elements", allow_duplicate=True),
            Output("cyto-popup-menu", "style", allow_duplicate=True),

            Input("btn-select-downstream", "n_clicks"),

            State("right-click-node-store", "data"),
            State("reactome-network", "elements"),
            State("all-elements-store", "data"),
            State("cyto-popup-menu", "style"),

            prevent_initial_call=True,
        )
        def select_downstream_and_close_popup(
            n_clicks,
            right_click_node,
            current_elements,
            all_elements,
            current_style,
        ):
            if not n_clicks or not right_click_node:
                raise dash.exceptions.PreventUpdate

            node_id = right_click_node.get("id")

            if not node_id:
                raise dash.exceptions.PreventUpdate

            downstream_ids = get_downstream_node_ids(
                node_id=node_id,
                all_elements=all_elements,
            )

            selected_ids = set(downstream_ids)
            selected_ids.add(node_id)

            new_elements = select_nodes_in_elements(
                current_elements,
                selected_ids=selected_ids,
            )

            style = dict(current_style or {})
            style["display"] = "none"

            return new_elements, style


        @app.callback(
            Output("reactome-network", "elements", allow_duplicate=True),
            Output("cyto-popup-menu", "style", allow_duplicate=True),

            Input("btn-select-hubs", "n_clicks"),

            State("reactome-network", "elements"),
            State("all-elements-store", "data"),
            State("cyto-popup-menu", "style"),

            prevent_initial_call=True,
        )
        def select_main_hubs_and_close_popup(
            n_clicks,
            current_elements,
            all_elements,
            current_style,
        ):
            if not n_clicks:
                raise dash.exceptions.PreventUpdate

            hub_ids = self.get_main_hub_node_ids_cached(
                all_elements=all_elements,
                top_n=10,
            )

            new_elements = select_nodes_in_elements(
                current_elements,
                selected_ids=hub_ids,
            )

            style = dict(current_style or {})
            style["display"] = "none"

            return new_elements, style

        @app.callback(
            Output("hover-node-store", "data"),
            Input("reactome-network", "mouseoverNodeData"),
            prevent_initial_call=True,
        )
        def store_hovered_node(node_data):
            if node_data is None:
                raise dash.exceptions.PreventUpdate

            return node_data
        
        """
            Right-click node A
                ↓
            All nodes deselected
                ↓
            Node A selected
                ↓
            Popup opens

            Click “Select 1st neighbors”
                ↓
            Node A + first neighbors selected
                ↓
            Popup closes        
        """

        @app.callback(
            Output("reactome-network", "elements", allow_duplicate=True),
            Output("cyto-popup-menu", "style", allow_duplicate=True),

            Input("btn-select-sources", "n_clicks"),

            State("reactome-network", "elements"),
            State("all-elements-store", "data"),
            State("cyto-popup-menu", "style"),

            prevent_initial_call=True,
        )
        def select_most_upstream_and_close_popup(
            n_clicks,
            current_elements,
            all_elements,
            current_style,
        ):
            if not n_clicks:
                raise dash.exceptions.PreventUpdate

            source_ids, _ = self.get_source_sink_node_ids_cached(
                all_elements=all_elements,
            )

            new_elements = select_nodes_in_elements(
                current_elements,
                selected_ids=source_ids,
            )

            style = dict(current_style or {})
            style["display"] = "none"

            return new_elements, style

        @app.callback(
            Output("reactome-network", "elements", allow_duplicate=True),
            Output("cyto-popup-menu", "style", allow_duplicate=True),

            Input("btn-select-sinks", "n_clicks"),

            State("reactome-network", "elements"),
            State("all-elements-store", "data"),
            State("cyto-popup-menu", "style"),

            prevent_initial_call=True,
        )
        def select_most_downstream_and_close_popup(
            n_clicks,
            current_elements,
            all_elements,
            current_style,
        ):
            if not n_clicks:
                raise dash.exceptions.PreventUpdate

            _, sink_ids = self.get_source_sink_node_ids_cached(
                all_elements=all_elements,
            )

            new_elements = select_nodes_in_elements(
                current_elements,
                selected_ids=sink_ids,
            )

            style = dict(current_style or {})
            style["display"] = "none"

            return new_elements, style


        @app.callback(
            Output("cyto-popup-menu", "style"),
            Output("right-click-node-store", "data"),
            Output("reactome-network", "elements", allow_duplicate=True),
            Output("node-info", "children", allow_duplicate=True),
            Output("selected-node-store", "data", allow_duplicate=True),

            Output("btn-ot-pathways", "style"),
            Output("btn-ot-diseases", "style"),
            Output("btn-ot-drugs", "style"),

            Input("right-click-event-store", "data"),
            Input("btn-close-popup", "n_clicks"),

            State("hover-node-store", "data"),
            State("cyto-popup-menu", "style"),
            State("reactome-network", "elements"),

            prevent_initial_call=True,
        )
        def show_or_hide_popup(
            right_click_event,
            close_clicks,
            hover_node_data,
            current_style,
            current_elements,
        ):
            trigger = dash.callback_context.triggered[0]["prop_id"].split(".")[0]

            style = dict(current_style or {})

            hidden_ot_style = {"display": "none"}

            visible_ot_style = {
                "display": "block",
                "width": "100%",
                "textAlign": "left",
            }

            if trigger == "btn-close-popup":
                style["display"] = "none"
                return (
                    style,
                    dash.no_update,
                    dash.no_update,
                    dash.no_update,
                    dash.no_update,
                    hidden_ot_style,
                    hidden_ot_style,
                    hidden_ot_style,
                )

            if trigger == "right-click-event-store":
                if not right_click_event or not hover_node_data:
                    raise dash.exceptions.PreventUpdate

                node_id = hover_node_data.get("id")

                if not node_id:
                    raise dash.exceptions.PreventUpdate

                x = right_click_event.get("x", 20)
                y = right_click_event.get("y", 70)

                style.update(
                    {
                        "display": "block",
                        "position": "fixed",
                        "left": f"{x}px",
                        "top": f"{y}px",
                        "zIndex": 9999,
                        "backgroundColor": "white",
                        "border": "1px solid #ddd",
                        "borderRadius": "12px",
                        "boxShadow": "0 8px 24px rgba(0,0,0,0.18)",
                        "padding": "10px",
                        "width": "230px",
                    }
                )

                new_elements = force_select_single_node(
                    current_elements,
                    node_id=node_id,
                )

                node_info_panel = self.make_node_info_panel(hover_node_data)

                if is_gene_node(hover_node_data):
                    ot_pathways_style = visible_ot_style
                    ot_diseases_style = visible_ot_style
                    ot_drugs_style = visible_ot_style
                else:
                    ot_pathways_style = hidden_ot_style
                    ot_diseases_style = hidden_ot_style
                    ot_drugs_style = hidden_ot_style

                return (
                    style,
                    hover_node_data,
                    new_elements,
                    node_info_panel,
                    hover_node_data,
                    ot_pathways_style,
                    ot_diseases_style,
                    ot_drugs_style,
                )

            raise dash.exceptions.PreventUpdate

        @app.callback(
            Output("opentarget-action-store", "data"),

            Input("btn-ot-pathways", "n_clicks"),
            Input("btn-ot-diseases", "n_clicks"),
            Input("btn-ot-drugs", "n_clicks"),

            prevent_initial_call=True,
        )
        def choose_opentarget_action(n_pathways, n_diseases, n_drugs):
            trigger = dash.callback_context.triggered[0]["prop_id"].split(".")[0]

            if trigger == "btn-ot-pathways":
                return "pathways"

            if trigger == "btn-ot-diseases":
                return "diseases"

            if trigger == "btn-ot-drugs":
                return "drugs"

            raise dash.exceptions.PreventUpdate


        @app.callback(
            Output("opentarget-modal", "is_open"),
            Output("opentarget-modal-title", "children"),
            Output("opentarget-modal-message", "children"),
            Output("opentarget-modal-table", "children"),
            Output("cyto-popup-menu", "style", allow_duplicate=True),

            Input("opentarget-action-store", "data"),
            Input("btn-close-opentarget-modal", "n_clicks"),

            State("right-click-node-store", "data"),
            State("opentarget-modal", "is_open"),
            State("cyto-popup-menu", "style"),

            prevent_initial_call=True,
        )
        def open_opentarget_modal(
            action,
            close_clicks,
            right_click_node,
            is_open,
            popup_style,
        ):
            trigger = dash.callback_context.triggered[0]["prop_id"].split(".")[0]

            if trigger == "btn-close-opentarget-modal":
                return False, dash.no_update, dash.no_update, dash.no_update, dash.no_update

            if not action:
                raise dash.exceptions.PreventUpdate

            popup_style = dict(popup_style or {})
            popup_style["display"] = "none"

            if not is_gene_node(right_click_node):
                return (
                    True,
                    "Open Targets",
                    "Selected node is not a valid gene node.",
                    html.Div("A valid gene node must have both symbol and Ensembl gene ID."),
                    popup_style,
                )

            symbol = get_gene_symbol_from_node(right_click_node)

            try:
                if action == "pathways":
                    target_id, df = self.ot.get_reactome_pathways_for_target(symbol)

                    title = f"Open Targets: Reactome pathways — {symbol}"
                    message = f"Target ID: {target_id}" if target_id else "Target ID not found."

                    return (
                        True,
                        title,
                        message,
                        df_to_dash_table(df),
                        popup_style,
                    )

                if action == "diseases":
                    target_id, df = self.ot.get_reactome_disease_evidence_for_target(symbol)

                    title = f"Open Targets: disease evidence — {symbol}"
                    message = f"Target ID: {target_id}" if target_id else "Target ID not found."

                    return (
                        True,
                        title,
                        message,
                        df_to_dash_table(df),
                        popup_style,
                    )

                if action == "drugs":
                    df = self.ot.get_drugs_for_disease(
                        disease=self.disease,
                        limit=None,
                    )

                    title = f"Open Targets: drugs for TCGA disease — {self.disease}"
                    message = f"Current project: {self.psi_id}; disease query: {self.disease}"

                    return (
                        True,
                        title,
                        message,
                        df_to_dash_table(df),
                        popup_style,
                    )

            except Exception as e:
                return (
                    True,
                    "Open Targets error",
                    f"Action: {action}; selected gene: {symbol}",
                    html.Pre(
                        str(e),
                        style={
                            "backgroundColor": "#fff5f5",
                            "color": "#b00020",
                            "padding": "12px",
                            "borderRadius": "8px",
                            "whiteSpace": "pre-wrap",
                        },
                    ),
                    popup_style,
                )

            raise dash.exceptions.PreventUpdate


        @app.callback(
            Output("reactome-network", "elements", allow_duplicate=True),
            Output("cyto-popup-menu", "style", allow_duplicate=True),

            Input("btn-select-neighbors", "n_clicks"),

            State("right-click-node-store", "data"),
            State("reactome-network", "elements"),
            State("all-elements-store", "data"),
            State("cyto-popup-menu", "style"),

            prevent_initial_call=True,
        )
        def select_first_neighbors_and_close_popup(
            n_clicks,
            right_click_node,
            current_elements,
            all_elements,
            current_style,
        ):
            if not n_clicks or not right_click_node:
                raise dash.exceptions.PreventUpdate

            node_id = right_click_node.get("id")

            if not node_id:
                raise dash.exceptions.PreventUpdate

            neighbors = self.get_neighbor_node_ids(
                node_id=node_id,
                all_elements=all_elements,
            )

            selected_ids = set(neighbors)
            selected_ids.add(node_id)

            new_elements = select_nodes_in_elements(
                current_elements,
                selected_ids=selected_ids,
            )

            style = dict(current_style or {})
            style["display"] = "none"

            return new_elements, style


        @app.callback(
            Output("node-info", "children"),
            Output("selected-node-store", "data"),
            Input("reactome-network", "tapNodeData"),
        )
        def show_node_info(node_data):
            if node_data is None:
                return "Click a node to see details.", None
            
            return self.make_node_info_panel(node_data), node_data


        @app.callback(
            Output("main-cyto-layout", "style"),
            Output("node-sidebar", "style"),
            Output("toggle-node-panel-button", "children"),
            Input("toggle-node-panel-button", "n_clicks"),
        )
        def toggle_node_sidebar(n_clicks):
            if n_clicks is None:
                n_clicks = 0

            hidden = n_clicks % 2 == 1

            if hidden:
                return (
                    {
                        "width": "100%",
                        "display": "grid",
                        "gridTemplateColumns": "1fr",
                        "gap": "0px",
                        "alignItems": "start",
                    },
                    {"display": "none"},
                    "Show node info",
                )

            return (
                {
                    "width": "100%",
                    "display": "grid",
                    "gridTemplateColumns": "5fr 1fr",
                    "gap": "12px",
                    "alignItems": "start",
                },
                {},
                "Hide node info",
            )


        @app.callback(
            Output("reactome-network", "elements", allow_duplicate=True),
            Output("reactome-network", "layout", allow_duplicate=True),
            Output("reactome-network", "zoom", allow_duplicate=True),
            Output("expanded-nodes-store", "data", allow_duplicate=True),
            Output("selected-node-store", "data", allow_duplicate=True),
            Output("layout-dropdown", "value", allow_duplicate=True),
            Input("reset-graph-button", "n_clicks"),
            State("all-elements-store", "data"),
            prevent_initial_call=True,
        )
        def reset_graph(n_clicks, all_elements):
            if not n_clicks:
                return no_update, no_update, no_update, no_update, no_update, no_update

            return (
                all_elements,          # restore original graph
                {"name": "preset"},    # restore saved positions
                1.0,                   # reset zoom
                [],                    # no expanded nodes
                None,                  # no selected node
                "preset",              # reset dropdown
            )

        @app.callback(
            Output("reactome-network", "layout"),
            Input("layout-dropdown", "value"),
        )
        def update_layout(layout_name):
            if layout_name is None:
                layout_name = "preset"

            if layout_name == "preset":
                return {"name": "preset"}

            if layout_name == "cose":
                return {
                    "name": "cose",
                    "animate": True,
                    "fit": True,
                    "padding": 30,
                    "nodeRepulsion": 8000,
                    "idealEdgeLength": 80,
                }

            if layout_name == "breadthfirst":
                return {
                    "name": "breadthfirst",
                    "directed": True,
                    "fit": True,
                    "padding": 30,
                    "spacingFactor": 1.2,
                }

            if layout_name == "circle":
                return {
                    "name": "circle",
                    "fit": True,
                    "padding": 30,
                }

            if layout_name == "grid":
                return {
                    "name": "grid",
                    "fit": True,
                    "padding": 30,
                }

            return {"name": "preset"}


        @app.callback(
            Output("save-toast", "is_open"),
            Output("save-toast", "children"),
            Output("saved-output", "children"),
            Input("save-graph-button", "n_clicks"),
            State("reactome-network", "elements"),
            prevent_initial_call=True,
        )
        def save_and_notify(n_clicks, elements):
            """
            save JSON
            ↓
            update self.saved_positions
            ↓
            rebuild elements with positions
            ↓
            return elements to Cytoscape
            """
            print("SAVE BUTTON CLICKED:", n_clicks)
            print("Number of elements:", len(elements) if elements else 0)
            
            changed = self.save_positions_if_changed(elements)

            if changed:
                message = f"Graph saved for: {self.pathway} - {self.pathway_id}"
                # new_elements = self.apply_positions_to_elements(elements, self.saved_positions,)
            else:
                message = "No changes detected"
                # new_elements = elements

            return True, message, message

        @app.callback(
            Output("cyto-font-size-store", "data"),
            Input("increase-font-btn", "n_clicks"),
            Input("decrease-font-btn", "n_clicks"),
            State("cyto-font-size-store", "data"),
            prevent_initial_call=True,
        )
        def update_font_size(n_inc, n_dec, current_size):

            if current_size is None:
                current_size = base_font_size

            button_id = ctx.triggered_id

            if button_id == "increase-font-btn":
                current_size += 1

            elif button_id == "decrease-font-btn":
                current_size -= 1

            print(">>> current_size", current_size)

            current_size = max(6, min(current_size, 30))

            return current_size

        @app.callback(
            Output("reactome-network", "stylesheet"),
            Output("font-size-label", "children"),
            Input("cyto-font-size-store", "data"),
        )
        def refresh_graph_font_size(font_size):
            if font_size is None:
                font_size = base_font_size

            font_size = int(font_size)
            font_size = max(6, min(font_size, 30))

            self.save_cyto_settings(font_size)

            return make_stylesheet(font_size), f"{font_size}px"

       
        @app.callback(
            Output("cyto-zoom-store", "data"),
            Input("cyto-zoom-slider", "value"),
        )
        def update_graph_zoom1(zoom_value):
            if zoom_value is None:
                zoom_value = 1.0

            return float(zoom_value)

        @app.callback(
            Output("reactome-network", "zoom"),
            Input("cyto-zoom-store", "data"),
        )
        def update_graph_zoom2(zoom_value):
            if zoom_value is None:
                zoom_value = 1.0

            return float(zoom_value)

        return app
    

    def kill_dash_ports(self, dry_run: bool = False):
        """
        Kill processes occupying Dash ports, but never kill Streamlit.
        Intended for local development.
        """

        for port in range(self.START_PORT, self.END_PORT + 1):
            result = subprocess.run(
                ["lsof", "-ti", f"tcp:{port}"],
                capture_output=True,
                text=True,
            )

            pids = result.stdout.strip().splitlines()

            if not pids:
                continue

            for pid in pids:
                info = subprocess.run(
                    ["ps", "-p", pid, "-o", "args="],
                    capture_output=True,
                    text=True,
                )

                cmdline = info.stdout.strip()
                cmdline_lower = cmdline.lower()

                print(f">>> PID {pid} on port {port}: {cmdline}")

                if "streamlit" in cmdline_lower:
                    print(f"Skipping PID {pid}: looks like Streamlit")
                    continue

                looks_like_dash = (
                    "dash" in cmdline_lower
                    or "flask" in cmdline_lower
                    or "werkzeug" in cmdline_lower
                    or "dashcyto_lib" in cmdline_lower
                    or "libs.dashcyto_lib" in cmdline_lower
                )

                if looks_like_dash:
                    if dry_run:
                        print(f"Would kill Dash PID {pid} on port {port}")
                    else:
                        print(f"Killing Dash PID {pid} on port {port}")
                        subprocess.run(["kill", "-9", pid])
                else:
                    print(f"Skipping PID {pid}: not clearly Dash")
                    
    def is_port_free(self, port: int, host: str = "127.0.0.1") -> bool:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            return s.connect_ex((host, port)) != 0


    def find_free_dash_port(self) -> int:

        for port in range(self.START_PORT, self.END_PORT + 1):
            if self.is_port_free(port):
                return port

        raise RuntimeError(
            f"No free Dash port found between {self.START_PORT} and {self.END_PORT}"
        )

    def run_app(self, height: str = "95%", width: str = "100%", marginTop: str = "20px", port: int | None = None,):

        if port is None:
            port = self.find_free_dash_port()
            
        app = self.create_cytoscape_app(height=height, width=width, marginTop=marginTop, port=port)
        url = f"http://localhost:{port}"

        threading.Timer(1.0, lambda: webbrowser.open(url)).start()
        print(f"Open in browser: {url}")

        app.run(
            debug=False,
            port=port,
            use_reloader=False,
            jupyter_mode="external",
        )

    def build_gene_alias_table(self, force:bool=False, verbose:bool=False) -> pd.DataFrame:
        """
        Build helper table mapping:
            official symbol -> official symbol
            synonym/alias   -> official symbol

        Save symbol table -> to synonyms lookup table
        """

        filename = self.root_ncbi / self.fname_gene_alias

        if filename.exists() and not force:
            return pdreadcsv(self.fname_gene_alias, self.root_ncbi, verbose=verbose)

        #--- read hugo table
        df = pdreadcsv(self.fname_hugo, self.root_ncbi, verbose=verbose)

        rows = []
        for _, row in df.iterrows():
            symbol = row.symbol.strip()

            if not symbol:
                continue

            base_record = {
                "symbol": symbol,
                "ensembl_id": None if pd.isnull(row.ensembl_id) else row.ensembl_id.strip(),
                "name": None if pd.isnull(row['name']) else row['name'].strip(),
                "uniprot_id": row.uniprot_id,
                "synonyms": row.synonyms,
                "refseq_summary": None if pd.isnull(row.refseq_summary) else row.refseq_summary.strip(),
            }

            # Include official symbol as its own alias
            rows.append(
                {
                    "alias": symbol,
                    "alias_upper": symbol.upper(),
                    **base_record,
                }
            )

            synonyms = str(row.get("synonyms", "")).strip()

            if synonyms:
                # Accept ;, |, or comma as separators
                for alias in re.split(r"[;|,]", synonyms):
                    alias = alias.strip()

                    if not alias:
                        continue

                    rows.append(
                        {
                            "alias": alias,
                            "alias_upper": alias.upper(),
                            **base_record,
                        }
                    )

        df_alias = pd.DataFrame(rows)

        if df_alias.empty:
            return df_alias

        df_alias = df_alias.drop_duplicates(
            subset=["alias_upper", "symbol", "ensembl_id"]
        ).reset_index(drop=True)

        print(f"Gene alias table rows: {len(df_alias):,}")

        _ = pdwritecsv(df_alias, self.fname_gene_alias, self.root_ncbi, verbose=verbose)

        return df_alias

