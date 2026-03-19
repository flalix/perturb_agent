#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2023/06/25
# Udated  on 2024/10/11; 2024/05/06; 2024/03/29; 2023/08/28; 2023/08/16
# @author: Flavio Lichtenstein
# @local: Bioinformatics: CENTD/Molecular Biology; Instituto Butatan

# import numpy as np
import requests
# import os, sys, shutil
from os.path import join as osjoin
from os.path import exists as exists
# import pandas as pd
# import pickle, sys, gzip
# from collections import Counter #, OrderedDict, defaultdict
from typing import  Tuple, List  # Optional, Iterable, Set, Any
# from datetime import datetime
# import psutil

# import matplotlib.pyplot as plt
# import matplotlib.colors as mpl_colors
# import matplotlib.gridspec as gridspec

# import plotly.express as px
# import plotly.graph_objects as go
# import plotly.figure_factory as ff
# import plotly.io as pio
# from   plotly.subplots import make_subplots

# from scipy.stats import shapiro
# from scipy import stats
# from scipy.stats import norm
# from scipy.stats import spearmanr

import py3Dmol

from pdbfixer import PDBFixer
from openmm.app import PDBFile
from openmm.app import ForceField
from openmm import VerletIntegrator
from openmm.app import Simulation
from openmm import unit as openmm_unit
from openmm import Platform

import MDAnalysis as mda

from docking_gnina.Basic import *


class PDB(object):
	def __init__(self, root_data:str='../data/'):
		
		self.root_data = root_data
		self.root_pdb	  = create_dir(root_data, 'pdb')
		self.root_docking = create_dir(root_data, 'docking_files')

		self.pdb_url	= "https://files.rcsb.org/download/%s.pdb"
		self.molssi_url = "https://raw.githubusercontent.com/MolSSI-Education/iqb-2025/refs/heads/main/data/%s.pdb"
		
		self.pdb_id = ''
		self.fixer = None


	def set_pdb_id(self, pdb_id:str, verbose:bool=False):
		self.pdb_id = pdb_id
		self.fixer = PDBFixer()
		
		if verbose: print(f"Setted: {pdb_id}")

		
	def get_pdb_by_id(self,  force:bool=False, verbose:bool=False) -> requests.Response:
		if self.pdb_id is None or self.pdb_id == '':
			print("PDB ID is not setted. Use set_pdb_id() method.")
			return requests.Response()
		
		protein_url = self.pdb_url%(self.pdb_id)
		
		try:
			protein_request = requests.get(protein_url)
		
			protein_request.raise_for_status()
		except:
			print(f"Bad request for {self.pdb_id}")
			return requests.Response()
	
		protein_filename = f"{self.pdb_id}.pdb"
		filename = osjoin(self.root_pdb, protein_filename)
	
		if exists(filename) and not force:
			if verbose: print(f"File already saved to '{filename}'")
		else:
			try:
				with open(filename, "w") as f:
					f.write(protein_request.text)
				if verbose: print(f"Saved PDB protein file to {filename}")
			except:
				print("Could not write the pdb file to '{filename}'")
				return requests.Response()
		
		return protein_request

	def get_pdb_text(self) -> str:
		protein_request = self.get_pdb_by_id()

		if protein_request is None or protein_request.status_code != 200:
			return "???"

		s_html = f"""<div style="height:400px; overflow:auto;"><pre>{protein_request.text}</pre></div>"""

		return s_html

		
	def py3Dmol_view(self, pdb_id:str, chains:List, colors:List,
					 setStyle_dic_list:List=[], zoomTo_model:int=0,
					 rotate_y:int=-25, rotate_z:int=90) -> py3Dmol.view:

		view = py3Dmol.view()
		
		try:
			view.addModel(open(f"{self.root_pdb}/{pdb_id}.pdb").read())
		except:
			print(f"Could not read {pdb_id}")
			return py3Dmol.view()
		
		for chain, color in zip(chains, colors):
			view.setStyle({'chain':chain}, {'cartoon': {'color': color}})
			
		for setStyle_dic in setStyle_dic_list:
			try:
				view.setStyle(setStyle_dic['resi'], setStyle_dic['stick'])
			except:
				pass
			
		if zoomTo_model >= 0:
			view.zoomTo({'model':zoomTo_model})

		view.rotate(rotate_y, "y")
		view.rotate(rotate_z, "z")

		return view
	
	def run_pdbfixer(self, removeHeterogens:bool=True, removeWater:bool=True,
					 force:bool=False, verbose:bool=False) -> Tuple[bool, object, object, object, object]:
		
		# Load the PDB into the PDBFixer class
		filename = f"{self.root_pdb}/{self.pdb_id}.pdb"
		fixer = PDBFixer(filename=filename)

		fixer.findMissingResidues()
		missing_residues = fixer.missingResidues

		fixer.findNonstandardResidues()
		nonstandardResidues = fixer.nonstandardResidues
		fixer.replaceNonstandardResidues()

		if removeHeterogens:
			fixer.removeHeterogens(keepWater=not removeWater)

		fixer.findMissingAtoms()
		missingAtoms = fixer.missingAtoms
		missingTerminals = fixer.missingTerminals

		fixer.addMissingAtoms()

		filename = f"{self.root_pdb}/{self.pdb_id}_fix_heavy.pdb"
		ret = True

		if exists(filename) and not force:
			if verbose: print("Fixed PDB already exists for {filename}")
		else:
			# Toplology, Positions, file stream, and keep chain ID's
			try:
				with open(filename, 'w') as fhandle:
					PDBFile.writeFile(fixer.topology, fixer.positions, fhandle, True)
			except:
				ret = False
				
		self.fixer = fixer

		return ret, missing_residues, nonstandardResidues, missingAtoms, missingTerminals


	
	def pdb_addMissingHydrogens(self, forcefield_xml:str="amber/protein.ff14SB.xml",
								force:bool=False, verbose:bool=False):
	
		forcefield = ForceField(forcefield_xml)
		# also accepts a "pH" argument, default is pH=7.0
		self.fixer.addMissingHydrogens(forcefield=forcefield)
		
		filename = f"{self.root_pdb}/{self.pdb_id}_all_atoms.pdb"
		ret = True

		if exists(filename) and not force:
			if verbose: print(f"Fixed PDB already exists for {filename}")
		else:
			try:
				with open(filename, 'w') as fhandle:
					# Toplology, Positions, file stream, and keep chain ID's
					PDBFile.writeFile(self.fixer.topology, self.fixer.positions, fhandle, True)
					if verbose: print(f"File saved: '{filename}'")
			except:
				print(f"Could not saved all atoms: {filename}")
				ret = False
				
		return ret


	def energy_minimization_simulation(self, forcefield_xml:str="amber/protein.ff14SB.xml", platform=None,
								unit=None, verbose:bool=False) -> object:
	
		if unit is None:
			unit = openmm_unit.picosecond
		
		forcefield = ForceField(forcefield_xml)
		system = forcefield.createSystem(self.fixer.topology)
		
		if verbose: print(system.getForces())
		
		# Use a generic VerletIntegrator
		integrator = VerletIntegrator(0.001 * unit)

		# Create simulation for minimization
		# Optional, if you have access to a CUDA GPU, comment out the next line and uncomment the one after it
		
		if platform is None:
			platform = None
		else:
			platform = Platform.getPlatformByName('CUDA')

		# Define the OpenMM Simulation object, which serves as a convenince wrapper for the OpenMM Context and Reporter objects
		simulation = Simulation(self.fixer.topology, system, integrator, platform)
		# Set the position of our atoms
		simulation.context.setPositions(self.fixer.positions)
		
		# Minimize energy
		print('Minimizing energy...')
		simulation.minimizeEnergy()
		print("------------- end -------------")
		
		# Get minimized positions. We have to copy the positions of the compiled object back into Python's memory
		minimized_positions = simulation.context.getState(getPositions=True).getPositions()

		# Write minimized structure to a PDB file
		filename = f"{self.root_pdb}/{self.pdb_id}_fixed.pdb"
		
		try:
			with open(filename, 'w') as output:
				PDBFile.writeFile(self.fixer.topology, minimized_positions, output)
			print(f'Minimization complete. Minimized structure saved to {filename}_fixed_simple.pdb')
		except:
			print("Could not save file at '{filename}'")
						
		return simulation

	s_ligand_select="segid A and record_type HETATM and not resname HOH"

	def calc_residues_ligands_MDAnalysis(self, verbose:bool=False,
									  ligand_select:str=s_ligand_select) -> Tuple[object, object]:
		"""
		Determine residues nearby the ligand with MDAnalysis
		"""

		# Load the original PDB
		u = mda.Universe(f"{self.root_pdb}/{self.pdb_id}.pdb")

		# Select atoms using the MDAnalysis selection language
		
		ligand = u.select_atoms(ligand_select)

		# Find and residues within a certain distance from the ligand
		active_site = u.select_atoms("around 3.5 group ligand and segid A",
									 periodic=False, ligand=ligand)  # Uses generic select_name=object as kwargs
		
		if verbose:
			print(active_site.residues.resids)
			print(ligand.residues.resids)
		
		return active_site, ligand
	

	
	def visualize_poses(self, protein_file:str, pose_file:str, receptor_radius:float=0.1,
		cognate_file:str='', cognate_radius:float=0.25,
		animate:bool=True, receptor_color:str='navy', cognate_color:str="yellow", 
		pose_color:str="green", highlight_residues:List=[], rotate:int=270,
		width:int=800, height:int=600):

		"""
		Creates a 3D visualization of a protein with multiple ligand poses from a single SDF file.

		Parameters:
		-----------
		protein_file : str
			Path to the protein PDB file to visualize
		pose_file : str
			Path to the SDF file containing multiple ligand poses
		cognate_file : str, optional
			Path to the cognate (static structure) file if available
		animate : bool, optional
			Whether to animate through the poses (default: True)
		cognate_color: str, optional
			Color for cognate. Default is yellow.
		pose_color: str, optional
			Color for ligand poses. Default is green.
		highlight_residues: list, optional
			List of residue numbers to highlight with VDW representation

		Returns:
		--------
		py3Dmol.view
			The 3D visualization view object
		"""
		# Create the viewer
		view = py3Dmol.view(width=width, height=height)

		# Add protein structure
		try:
			view.addModel(open(protein_file).read())
			# Set protein style
			view.setStyle({"cartoon": {}, "stick": {"colorscheme": receptor_color+"Carbon","radius": receptor_radius}})
		except:
			print(f"Could not open {protein_file}")
			return None


		# Add VDW representation for highlighted residues if specified
		for res in highlight_residues:
			view.addStyle(
				{"model": 0, "resi": str(res)},
				{"sphere": {"color": "magenta", "opacity": 0.7, "scale": 0.7}},
			)

		# Add cognate ligand if provided
		if cognate_file and cognate_file != '':
			try:
				view.addModel(open(cognate_file).read())
				view.setStyle({"model": 1}, {"stick": {"colorscheme": cognate_color+"Carbon", "radius": cognate_radius}},)
			except:
				print(f"Could not open {cognate_file}")
				return None


		# Add all poses from the SDF file
		model_offset = 1 if cognate_file else 0

		if pose_file is not None and pose_file != '':
			try:
				with open(pose_file) as f:
					pose_content = f.read()
			except:
				print(f"Could not open {pose_file}")
				return None

			if animate:
				# Add all poses as animation frames
				view.addModelsAsFrames(pose_content)
				view.setStyle(
					{"model": model_offset + 1},
					{"stick": {"colorscheme": pose_color+"Carbon"}},
				)
				view.animate({"interval": 1000})
			else:
				# Add as separate models
				view.addModel(pose_content)

				view.setStyle(
					{"model": model_offset + 1},
					{"stick": {"colorscheme": pose_color+"Carbon"}},
				)

			# Set view - zoom to the correct structure
			if cognate_file:
				# Zoom to cognate if provided
				zoom_model = 1
			else:
				# Zoom to docked poses if no cognate
				zoom_model = model_offset + 1
		else:
			zoom_model = model_offset

		view.zoomTo({"model": zoom_model})
		view.rotate(rotate)

		return view
			

	def visualize_receptor_ligand(self, receptor_file:str, ligand_file:str, 
		carbon_color:str='navy', receptor_radius:float=0.1, ligand_color:str='green', ligand_radius:float=0.15,
		cognate_file:str='', cognate_color:str="yellow", cognate_radius:float=.25,
		model:int=1, zoom_to_model:int=1, rotate:int=270, width:int=800, height:int=600):			

		view = py3Dmol.view(width=width, height=height)

		if receptor_file and receptor_file != '':
			try:
				view.addModel(open(receptor_file).read())
				view.setStyle({'cartoon':{},'stick':{"colorscheme": carbon_color+'Carbon', 'radius':receptor_radius}})
			except:
				print(f"Could not open {receptor_file}")
				return None

	
		if ligand_file and ligand_file != '':
			try:
				view.addModel(open(ligand_file).read())

				if receptor_file != '':
					view.setStyle({"model": model}, {"stick": {"color": ligand_color, "radius": ligand_radius}})
				else:
					view.setStyle({"stick": {"color": ligand_color}})  #  "radius": ligand_radius
			except:
				print(f"Could not open {ligand_file}")
				return None


		# Add cognate ligand if provided
		if cognate_file and cognate_file != '':
			try:
				view.addModel(open(cognate_file).read())
				view.setStyle({"model": model}, {"stick": {"color": cognate_color, "radius": cognate_radius}})
			except:
				print(f"Could not open {cognate_file}")
				return None
			
		if zoom_to_model >= 0:
			view.zoomTo({'model':zoom_to_model})
			view.rotate(rotate)
	
		return view



	def visualize_one_pose(self, receptor_file:str, carbon_color:str='navy', receptor_radius:float=0.1, 
		ligand_file:str='', ligand_color_scheme:str='lime', ligand_radius:float=0.15, the_content:int=1,
		zoom_to_model:bool=True, rotate:int=270, width:int=800, height:int=600):			

		view = py3Dmol.view(width=width, height=height)
		model = 1

		if receptor_file and receptor_file != '':
			try:
				view.addModel(open(receptor_file).read())
				view.setStyle({'cartoon':{},'stick':{"colorscheme": carbon_color+'Carbon', 'radius':receptor_radius}})
			except:
				print(f"Could not open {receptor_file}")
				return None


		# Add ligand ligand if provided
		if ligand_file and ligand_file != '':
			try:
				with open(ligand_file) as f:
					pose_content = f.read()
			except:
				print(f"Could not open pose file: '{ligand_file}'")
				return None

			poses = pose_content.split("ENDMDL")
			print(f">>> pose {the_content} / {len(poses)}")

			chosen_pose = poses[the_content] + "ENDMDL"

			try:
				view.addModel(chosen_pose)
				view.setStyle({"model": model}, {"stick": {"colorscheme": ligand_color_scheme+'Carbon', "radius": ligand_radius}})
			except:
				print(f"Could not open {ligand_file}")
				return None
			
		if zoom_to_model:
			view.zoomTo({'model': model})
			view.rotate(rotate)

		return view



	def visualize_a_ligand(self, ligand_file:str, 
		ligand_color:str='green', ligand_radius:float=0.15,
		width:int=800, height:int=600):

		view = py3Dmol.view(width=width, height=height)
		view.addModel(open(ligand_file).read(), "sdf")

		if ligand_color:
			if ligand_radius:
				view.setStyle({"stick": {"color": ligand_color, "radius": ligand_radius}})
		else:
			view.setStyle({"stick": {}})

		view.zoomTo()
		return view		
