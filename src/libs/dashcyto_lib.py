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
import threading
import webbrowser
from pathlib import Path

import dash
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
import networkx as nx
from dash import Input, Output, State, html
from rdflib import RDF, Graph, Namespace

# import py4cytoscape as p4c
from libs.Basic import create_dir


class DASH_CYTO(object):
    def __init__(self, ROOT0: Path):

        self.GENE_COLS = ["gene_id", "symbol", "gene_type"]

        self.ROOT0 = ROOT0

        self.root_src = create_dir(ROOT0, "src")
        self.root_owl = create_dir(ROOT0, "owl")

        self.fname_pos = "positions_%s.json"

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

    def create_cytoscape_app(
        self, height: str = "900px", width: str = "100%", marginTop: str = "20px"
    ):

        self.G = self.G.subgraph([n for n in self.G.nodes() if self.G.degree(n) > 0]).copy()

        elements = self.nx_to_cytoscape_elements()
        self.elements = elements

        app = dash.Dash(__name__)

        title = f"{self.pathway} - {self.pathway_id}"

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
                cyto.Cytoscape(
                    id="reactome-network",
                    elements=elements,
                    layout={"name": "preset"},
                    boxSelectionEnabled=True,   # drag box to select many nodes
                    autoungrabify=False,        # nodes can be moved
                    autounselectify=False,      # nodes can be selected
                    style={
                        "width": width,
                        "height": height,
                        "backgroundColor": "#fff9c4",  # light yellow
                        "marginTop": marginTop,
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
                            "style": {"background-color": "#f94144", "shape": "hexagon"},
                        },
                        {
                            "selector": '[biopax_type = "Protein"]',
                            "style": {"background-color": "#8ecae6", "shape": "ellipse"},
                        },
                        {
                            "selector": '[biopax_type = "Complex"]',
                            "style": {"background-color": "#ffb703", "shape": "round-rectangle"},
                        },
                        {
                            "selector": '[biopax_type = "SmallMolecule"]',
                            "style": {"background-color": "#90be6d", "shape": "diamond"},
                        },
                        {
                            "selector": '[biopax_type = "Pathway"]',
                            "style": {"background-color": "#cdb4db", "shape": "rectangle"},
                        },
                    ],
                ),
                html.Button("Save node positions", id="save-button"),
                html.Pre(id="saved-output"),
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
            Output("reactome-network", "elements"),
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

            changed = self.save_positions_if_changed(elements)

            if changed:
                message = f"Graph saved for: {self.pathway} - {self.pathway_id}"
                new_elements = self.nx_to_cytoscape_elements(self.saved_positions)
            else:
                message = "No changes detected"
                new_elements = elements

            return True, message, new_elements
        
        @app.callback(
            Output("saved-output", "children"),
            Input("reactome-network", "selectedNodeData"),
        )
        def show_selected_nodes(selected_nodes):
            """
            Streamlit itself is not doing the selection. 
            The selection happens inside the Dash Cytoscape browser app, and Streamlit can present/open that app.            
            """
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
