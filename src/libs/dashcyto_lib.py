#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2026/05/04
# Udated  on 2026/05/04
# @author: Flavio Lichtenstein
# @local: Home sweet home


import json
import math
import os
import re
import subprocess
from os import path
import pandas as pd
import threading
import webbrowser
import socket
from pathlib import Path

import dash
from dash import html, dcc, Input, Output, State, ctx, no_update
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
import networkx as nx
from rdflib import RDF, Graph, Namespace

# import py4cytoscape as p4c
from libs.Basic import create_dir, pdreadcsv, download_url_file


class DASH_CYTO(object):
    def __init__(self, root0: Path, root0_data: Path, dflfc_ori: pd.DataFrame):

        self.GENE_COLS = ["gene_id", "symbol", "gene_type"]

        self.root0 = root0
        self.root0_data = root0_data

        self.root_src = create_dir(root0, "src")
        self.root_styles = create_dir(self.root_src, "styles")

        self.root_colab= create_dir(root0_data, "colab")
        self.root_owl = create_dir(self.root_colab, "owl")
        self.root_ncbi = create_dir(self.root_colab, "ncbi")

        self.dflfc_ori = dflfc_ori
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
        self.kill_dash_ports()

    def reset_graph(self):
        self.G = nx.DiGraph()
        self.saved_positions = {}
        self.layout_changed = False

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
                self.G.add_node(node_id, label=self.get_name(node), biopax_type=self.short(cls))

        for rel in self.relations:
            for s, o in self.rdf.subject_objects(rel):
                s_id = self.short(s)
                o_id = self.short(o)

                if s_id in self.G.nodes and o_id in self.G.nodes:
                    self.G.add_edge(s_id, o_id, interaction=self.short(rel))

        return True

    def short(self, x) -> str:
        return str(x).split("#")[-1].split("/")[-1]

    def get_name(self, x) -> str:

        for prop in [self.BP.displayName, self.BP.standardName, self.BP.name]:
            value = next(self.rdf.objects(x, prop), None)
            if value:
                return str(value)

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


    def get_lfc_bin(self, log2FC: float | None) -> str:
        if log2FC is None:
            return "none"

        try:
            log2FC = float(log2FC)
        except Exception:
            return "none"

        if log2FC >= 3:
            return "up_3"
        elif log2FC >= 2:
            return "up_2"
        elif log2FC >= 1:
            return "up_1"
        elif log2FC >= 0.4:
            return "+weak"
        elif log2FC <= -3:
            return "down_3"
        elif log2FC <= -2:
            return "down_2"
        elif log2FC <= -1:
            return "down_1"
        elif log2FC <= -0.4:
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
            symbol = str(label).upper()

            dic_symbol = self.dic_lfc_lookup.get(symbol)

            #if symbol == 'CCNB2':
            #    print(f"Found CCNB2: {dic_symbol}")

            if dic_symbol is not None:
                log2FC = dic_symbol["lfc"]
                fdr = dic_symbol["fdr"]
            else:
                log2FC = None
                fdr = None


            elem = {
                "data": {
                    "id": node_id,
                    "label": data.get("label", node_id),
                    "biopax_type": data.get("biopax_type", "Unknown"),
                    "log2FC": log2FC,
                    "FDR": fdr,
                    "abs_log2FC": abs(log2FC) if log2FC is not None else None,
                    "lfc_bin": self.get_lfc_bin(log2FC),
                    "is_deg_gene": log2FC is not None,                    
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
        label = first_available("label", default=node_id)

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
                    html.P([html.B("LFC / log2FC: "), str(info_data["lfc"])]),
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

                        html.Button(
                            "Hide node info",
                            id="toggle-node-panel-button",
                            n_clicks=0,
                            className="cyto-button",
                            style={
                                "fontSize": "12px",
                                "padding": "6px 12px",
                            },
                        ),
                    ],
                    style={
                        "display": "flex",
                        "alignItems": "center",
                        "justifyContent": "space-between",
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
                            ],
                            id="graph-panel",
                            style={"minWidth": "0"},
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

                                # reset/save buttons
                                html.Button(
                                    "🔄 Reset graph",
                                    id="reset-graph-button",
                                    n_clicks=0,
                                    className="cyto-button cyto-button-reset",
                                ),

                                html.Button(
                                    "💾 Save node positions",
                                    id="save-graph-button",
                                    n_clicks=0,
                                    className="cyto-button cyto-button-save",
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
                dcc.Store(id="expanded-nodes-store", data=[]),
                dcc.Store(id="all-elements-store", data=elements),

                html.Pre(id="saved-output"),
                html.Pre(id="selected-output"),

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
            ]
        )

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
            Output("node-info", "children"),
            Output("selected-node-store", "data"),
            Input("reactome-network", "tapNodeData"),
        )
        def show_node_info(node_data):
            if node_data is None:
                return "Click a node to see details.", None

            return self.make_node_info_panel(node_data), node_data


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

