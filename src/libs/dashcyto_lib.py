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
from os import path
import pandas as pd
import threading
import webbrowser
from pathlib import Path

import dash
from dash import html, dcc, Input, Output, State, ctx
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
import networkx as nx
from rdflib import RDF, Graph, Namespace

# import py4cytoscape as p4c
from libs.Basic import create_dir, pdreadcsv


class DASH_CYTO(object):
    def __init__(self, ROOT0: Path):

        self.GENE_COLS = ["gene_id", "symbol", "gene_type"]

        self.ROOT0 = ROOT0

        self.root_src = create_dir(ROOT0, "src")
        self.root_styles = create_dir(self.root_src, "styles")
        self.root_owl = create_dir(ROOT0, "data/owl")
        self.root_ncbi = create_dir(ROOT0, "data/ncbi")

        self.fname_pos = "positions_%s.json"

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


    def read_owl(self, pathway_id: str, pathway: str) -> Graph:
        self.pathway_id = pathway_id
        self.pathway = pathway

        fname_owl = f"{pathway_id}_level3.owl"
        filename = self.root_owl / fname_owl

        rdf = Graph()
        try:
            rdf.parse(filename, format="xml")
        except Exception as e:
            print(f"Error parsing OWL file: {e}")

        self.rdf = rdf

        for cls in self.classes:
            for node in rdf.subjects(RDF.type, cls):
                node_id = self.short(node)
                self.G.add_node(node_id, label=self.get_name(node), biopax_type=self.short(cls))

        for rel in self.relations:
            for s, o in rdf.subject_objects(rel):
                s_id = self.short(s)
                o_id = self.short(o)

                if s_id in self.G.nodes and o_id in self.G.nodes:
                    self.G.add_edge(s_id, o_id, interaction=self.short(rel))

        return rdf

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

    def nx_to_cytoscape_elements(self, saved_positions=None) -> list:
        elements = []

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
            elem = {
                "data": {
                    "id": node_id,
                    "label": data.get("label", node_id),
                    "biopax_type": data.get("biopax_type", "Unknown"),
                    "log2FC": data.get("log2FC", None),
                    "FDR": data.get("FDR", None),
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

        if not self.positions_changed(new_pos):
            return False  # no change

        fname = self.fname_pos % (self.pathway_id)
        filename = self.root_owl / fname

        try:
            with open(filename, "w") as f:
                json.dump(new_pos, f, indent=2)

            self.saved_positions = new_pos  # update global here
        except Exception as e:
            print(f"Error: saving Graph positions: {e}")
            return False

        return True
    
    #============ Hugo, Ensbembl, UNIPROT gene table ===========

    def load_gene_annotation_table(self) -> pd.DataFrame:
        """
        Load HGNC/HUGO + RefSeq + UniProt annotation table.

        Expected columns:
            ensembl_id, symbol, name, uniprot_id,
            ncbi_gene_id, synonyms, refseq_summary
        """

        filename = self.root_ncbi / self.fname_hugo
        if not filename.exists():
            raise FileNotFoundError(f"Could not find gene annotation file: {filename}")

        df = pdreadcsv(self.fname_hugo, self.root_ncbi)

        required_cols = [
            "ensembl_id",
            "symbol",
            "name",
            "uniprot_id",
            "synonyms",
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
        df["synonyms"] = df["synonyms"].astype(str).str.strip()
        df["refseq_summary"] = df["refseq_summary"].astype(str).str.strip()

        df["ensembl_id_base"] = df["ensembl_id"].str.replace(
            r"\.\d+$",
            "",
            regex=True,
        )

        self.gene_annot_df = df

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

        print(f"Loaded gene annotation table: {filename}")
        print(f"Rows: {len(df):,}")
        print(f"Symbols: {len(self.gene_annot_by_symbol):,}")
        print(f"Ensembl IDs: {len(self.gene_annot_by_ensembl):,}")

        return df

    def load_gene_alias_table(self) -> pd.DataFrame:
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

        print(f"Loaded alias rows: {len(self.gene_alias_df):,}")
        print(f"Unique aliases: {self.gene_alias_df.index.nunique():,}")

        print("DEBUG TSP1:")
        print(self.gene_alias_df.loc["TSP1"])

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


    def get_gene_annotation_for_node(self, node_data: dict) -> dict:
        empty = {
            "ensembl_id": "NA",
            "symbol": "NA",
            "name": "NA",
            "uniprot_id": "NA",
            "refseq_summary": "NA",
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
            print(f"No gene annotation match for label: {label}")
            return empty

        first = dict(matches[0])
        first["matches"] = matches
        first["n_matches"] = len(matches)

        return first


    def get_gene_annotation_for_node2(self, node_data: dict) -> dict:
        """
        Find gene annotation for a Cytoscape node.

        Lookup priority:
            1. Ensembl ID, if the node has a real ENSG ID
            2. Official symbol / synonym using node label
            3. Symbol-like fields

        Returns the first match plus all matches if the alias is ambiguous.
        """

        empty = {
            "ensembl_id": "",
            "symbol": "",
            "name": "",
            "uniprot_id": "",
            "synonyms": "",
            "refseq_summary": "",
            "matches": [],
            "n_matches": 0,
        }

        # -----------------------------
        # 1. Try Ensembl ID only if real
        # -----------------------------
        ensembl_id = (
            node_data.get("ensembl_id")
            or node_data.get("gene_id")
            or ""
        )

        ensembl_id = str(ensembl_id).strip()
        ensembl_id_base = ensembl_id.split(".")[0]

        if (
            ensembl_id_base.startswith("ENSG")
            and hasattr(self, "gene_annot_by_ensembl")
            and ensembl_id_base in self.gene_annot_by_ensembl
        ):
            record = self.gene_annot_by_ensembl[ensembl_id_base]
            record = dict(record)
            record["matches"] = [record]
            record["n_matches"] = 1
            return record

        # -----------------------------
        # 2. Try label/symbol via alias table
        # -----------------------------
        label = (
            node_data.get("label")
            or node_data.get("symbol")
            or node_data.get("gene_symbol")
            or node_data.get("name")
            or ""
        )

        label = str(label).strip()

        if label and hasattr(self, "gene_alias_lookup"):
            matches = self.lookup_gene_alias(label)

            if matches:
                first = dict(matches[0])
                first["matches"] = matches
                first["n_matches"] = len(matches)
                return first

        # -----------------------------
        # 3. Fallback: official symbol table
        # -----------------------------
        if label and hasattr(self, "gene_annot_by_symbol"):
            if label in self.gene_annot_by_symbol:
                record = dict(self.gene_annot_by_symbol[label])
                record["matches"] = [record]
                record["n_matches"] = 1
                return record

            label_upper = label.upper()
            for symbol, record in self.gene_annot_by_symbol.items():
                if str(symbol).upper() == label_upper:
                    record = dict(record)
                    record["matches"] = [record]
                    record["n_matches"] = 1
                    return record

        return empty


    def none_one_multiple_alias(self, info_data: dict):

        print(">>> info_data", info_data)

        if info_data.get("symbol", "NA") == "NA":
            return None
        
        matches = info_data.get("matches", [])
        n_matches = info_data.get("n_matches", 0)

        print(">>> n_matches", n_matches)

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

        synonyms = (
            annot.get("synonyms")
            or first_available("synonyms", "alias", "aliases", default="NA")
        )

        refseq_summary = (
            annot.get("refseq_summary", "NA")
        )

        gene_type = first_available(
            "gene_type",
            "biotype",
        )

        lfc = first_available(
            "log2FoldChange",
            "log2fc",
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


        print("ANNOT:", annot)
        print("matches:", matches)
        print("n_matches:", n_matches)

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

        return html.Div(
            [
                html.P([html.B("ID: "), str(info_data["id"])]),
                html.P([html.B("Label: "), str(info_data["label"])]),

                html.Hr(),
                alias_note,
        
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
            ]
        )

    def create_cytoscape_app(self, height: str = "95%", width: str = "100%", marginTop: str = "20px"):

        self.G = self.G.subgraph([n for n in self.G.nodes() if self.G.degree(n) > 0]).copy()

        elements = self.nx_to_cytoscape_elements()
        self.elements = elements

        app = dash.Dash(__name__, 
                        assets_folder=str(self.root_styles),
                        external_stylesheets=[dbc.themes.BOOTSTRAP],
                        suppress_callback_exceptions=True,
                        )

        title = f"{self.pathway} - {self.pathway_id}"

        """
        save_and_notify          -> updates toast + saved-output
        show_node_info           -> updates node-info + selected-node-store
        update_graph_visibility  -> updates elements + expanded-nodes-store
        update_layout            -> updates layout
        show_selected_nodes      -> updates selected-output
        """

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

                html.Div(
                    [
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
                                        "marginBottom": "10px",
                                    },
                                ),

                                cyto.Cytoscape(
                                    id="reactome-network",
                                    elements=elements,
                                    layout={"name": "preset"},
                                    boxSelectionEnabled=True,
                                    autoungrabify=False,
                                    autounselectify=False,
                                    style={
                                        "width": "100%",
                                        "height": "80vh",   # important
                                        "backgroundColor": "#fff9c4",
                                        "marginTop": marginTop,
                                        "border": "1px solid #ddd",
                                        "borderRadius": "12px",
                                    },
                                    stylesheet=[
                                        {
                                            "selector": "node",
                                            "style": {
                                                "label": "data(label)",
                                                "font-size": "10px",
                                                "background-color": "#8ecae6",
                                                "width": 30,
                                                "height": 30,
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
                                                "font-size": "7px",
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
                                    ],
                                ),
                            ],
                            style={
                                "flex": "1 1 auto",
                                "minWidth": "0",
                            },
                        ),

                        html.Div(
                            [
                                html.H4("Node information", className="node-panel-title"),

                                html.Div(id="node-info", className="node-info-box"),

                                html.Button(
                                    "⭐ Mark as root",
                                    id="mark-root-button",
                                    n_clicks=0,
                                    className="cyto-button cyto-button-root",
                                ),

                                html.Button(
                                    "🏁 Mark as terminal",
                                    id="mark-terminal-button",
                                    n_clicks=0,
                                    className="cyto-button cyto-button-terminal",
                                ),

                                html.Button(
                                    "🕸️ Expand neighbors",
                                    id="expand-neighbors-button",
                                    n_clicks=0,
                                    className="cyto-button cyto-button-expand",
                                ),

                                html.Button(
                                    "🔄 Reset graph",
                                    id="reset-graph-button",
                                    n_clicks=0,
                                    className="cyto-button cyto-button-reset",
                                ),
                                                                
                            ],
                            className="node-sidebar",
                        ),
                    ],
                    style={
                        "width": "100%",
                        "display": "flex",
                        "gap": "12px",
                        "alignItems": "flex-start",
                    },
                ),

                dcc.Store(id="selected-node-store"),
                dcc.Store(id="expanded-nodes-store", data=[]),
                dcc.Store(id="all-elements-store", data=elements),

                html.Button("💾 Save node positions", id="save-button", n_clicks=0),

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
            Output("save-toast", "is_open"),
            Output("save-toast", "children"),
            Output("saved-output", "children"),
            Input("save-button", "n_clicks"),
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
            Output("reactome-network", "elements"),
            Output("expanded-nodes-store", "data"),
            Input("expand-neighbors-button", "n_clicks"),
            Input("reset-graph-button", "n_clicks"),
            State("reactome-network", "selectedNodeData"),
            State("reactome-network", "elements"),
            State("expanded-nodes-store", "data"),
            State("all-elements-store", "data"),
            prevent_initial_call=True,
        )
        def update_graph_visibility(
            expand_clicks,
            reset_clicks,
            selected_node_data,
            current_elements,
            expanded_nodes,
            all_elements,
        ):
            trigger = ctx.triggered_id

            if trigger == "reset-graph-button":
                print("RESET GRAPH CLICKED")
                return all_elements, []

            if trigger == "expand-neighbors-button":
                if not selected_node_data:
                    return current_elements, expanded_nodes

                selected_node = selected_node_data[0]

                new_elements, new_expanded_nodes = self.toggle_expand_contract(
                    selected_node=selected_node,
                    current_elements=current_elements,
                    expanded_nodes=expanded_nodes,
                    all_elements=all_elements,
                )

                return new_elements, new_expanded_nodes

            return current_elements, expanded_nodes

        @app.callback(
            Output("reactome-network", "layout"),
            Input("layout-dropdown", "value"),
        )
        def update_layout(layout_name):

            if layout_name == "cose":
                return {
                    "name": "cose",
                    "animate": True,
                    "fit": True,
                    "padding": 100,
                    "nodeRepulsion": 15000,
                    "idealEdgeLength": 180,
                    "edgeElasticity": 100,
                    "gravity": 40,
                    "numIter": 3000,
                }

            if layout_name == "breadthfirst":
                return {
                    "name": "breadthfirst",
                    "directed": True,
                    "fit": True,
                    "padding": 100,
                    "spacingFactor": 1.8,
                }

            if layout_name == "circle":
                return {
                    "name": "circle",
                    "fit": True,
                    "padding": 100,
                }

            if layout_name == "grid":
                return {
                    "name": "grid",
                    "fit": True,
                    "padding": 100,
                }

            return {"name": "preset"}

        @app.callback(
            Output("selected-output", "children"),
            Input("reactome-network", "selectedNodeData"),
        )
        def show_selected_nodes(selected_nodes):
            if not selected_nodes:
                return "No nodes selected"

            ids = [n["id"] for n in selected_nodes]
            return f"Selected nodes: {', '.join(ids)}"

        return app

    def run_app(
        self, height: str = "95%", width: str = "100%", marginTop: str = "20px", port: int = 8050
    ):
        app = self.create_cytoscape_app(height=height, width=width, marginTop=marginTop)
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

