#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2011/09/09
# Udated  on 2014/05/13, 2015/10/02, 2026/05/07
# @author: Flavio Lichtenstein


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 

import numpy as np
import math, os, sys
import pandas as pd
import plotly.graph_objects as go
# import plotly.express as px

from scipy.cluster.hierarchy import linkage, dendrogram

from libs.Basic import title_replace # create_dir, pdwritecsv, pdreadcsv

plotly_colors = ["aliceblue", "antiquewhite", "aqua", "aquamarine", "azure", "beige", "bisque", "black", "blanchedalmond", "blue", "blueviolet",
				 "brown", "burlywood", "cadetblue", "chartreuse", "chocolate", "coral", "cornflowerblue", "cornsilk", "crimson", "cyan",
				 "darkblue", "darkcyan", "darkgoldenrod", "darkgray", "darkgrey", "darkgreen", "darkkhaki", "darkmagenta", "darkolivegreen",
				 "darkorange", "darkorchid", "darkred", "darksalmon", "darkseagreen", "darkslateblue", "darkslategray", "darkslategrey",
				 "darkturquoise", "darkviolet", "deeppink", "deepskyblue", "dimgray", "dimgrey", "dodgerblue", "firebrick", "floralwhite",
				 "forestgreen", "fuchsia", "gainsboro", "ghostwhite", "gold", "goldenrod", "gray", "grey", "green", "greenyellow", "honeydew",
				 "hotpink", "indianred", "indigo", "ivory", "khaki", "lavender", "lavenderblush", "lawngreen", "lemonchiffon", "lightblue",
				 "lightcoral", "lightcyan", "lightgoldenrodyellow", "lightgray", "lightgrey", "lightgreen", "lightpink", "lightsalmon",
				 "lightseagreen", "lightskyblue", "lightslategray", "lightslategrey", "lightsteelblue", "lightyellow", "lime", "limegreen",
				 "linen", "magenta", "maroon", "mediumaquamarine", "mediumblue", "mediumorchid", "mediumpurple", "mediumseagreen",
				 "mediumslateblue", "mediumspringgreen", "mediumturquoise", "mediumvioletred", "midnightblue", "mintcream", "mistyrose",
				 "moccasin", "navajowhite", "navy", "oldlace", "olive", "olivedrab", "orange", "orangered", "orchid", "palegoldenrod",
				 "palegreen", "paleturquoise", "palevioletred", "papayawhip", "peachpuff", "peru", "pink", "plum", "powderblue", "purple",
				 "red", "rosybrown", "royalblue", "rebeccapurple", "saddlebrown", "salmon", "sandybrown", "seagreen", "seashell", "sienna",
				 "silver", "skyblue", "slateblue", "slategray", "slategrey", "snow", "springgreen", "steelblue", "tan", "teal", "thistle",
				 "tomato", "turquoise", "violet", "wheat", "white", "whitesmoke", "yellow", "yellowgreen"]

plotly_colors_proteins = ['olivedrab', 'navy', 'red', 'darkcyan', 'gold', 'mediumvioletred',
						  'green', 'darkgreen', 'orange', 'brown', 'gray', 'pink',
						  'mediumvioletred', 'darksalmon', 'lightgray', 'yellowgreen', 'magenta', 'darkturquoise',
						  'lime', 'orange', 'hotpink', 'indigo', "magenta", "maroon", "gold", "black"] + \
						  plotly_colors

plotly_my_colors = ['olivedrab', 'navy', 'red', 'darkcyan', 'gold', 'orange', 'pink',
					'mediumvioletred', 'green', 'darkgreen',  'brown', 'gray', 
					'mediumvioletred', 'darksalmon', 'gray', 'yellowgreen', 'magenta', 'darkturquoise',
					'lime', 'orange', 'hotpink', 'indigo', "magenta", "maroon", "gold", "black"] + \
					 plotly_colors

colorscales =   ['aggrnyl', 'agsunset', 'algae', 'amp', 'armyrose', 'balance',
				 'blackbody', 'bluered', 'blues', 'blugrn', 'bluyl', 'brbg',
				 'brwnyl', 'bugn', 'bupu', 'burg', 'burgyl', 'cividis', 'curl',
				 'darkmint', 'deep', 'delta', 'dense', 'earth', 'edge', 'electric',
				 'emrld', 'fall', 'geyser', 'gnbu', 'gray', 'greens', 'greys',
				 'haline', 'hot', 'hsv', 'ice', 'icefire', 'inferno', 'jet',
				 'magenta', 'magma', 'matter', 'mint', 'mrybm', 'mygbm', 'oranges',
				 'orrd', 'oryel', 'oxy', 'peach', 'phase', 'picnic', 'pinkyl',
				 'piyg', 'plasma', 'plotly3', 'portland', 'prgn', 'pubu', 'pubugn',
				 'puor', 'purd', 'purp', 'purples', 'purpor', 'rainbow', 'rdbu',
				 'rdgy', 'rdpu', 'rdylbu', 'rdylgn', 'redor', 'reds', 'solar',
				 'spectral', 'speed', 'sunset', 'sunsetdark', 'teal', 'tealgrn',
				 'tealrose', 'tempo', 'temps', 'thermal', 'tropic', 'turbid',
				 'turbo', 'twilight', 'viridis', 'ylgn', 'ylgnbu', 'ylorbr',
				 'ylorrd']

def define_delta_y(val_max_min:float, general_max:float):

	if val_max_min < 0:
		is_negative = True
		val_max_min = - val_max_min
	else:
		is_negative = False

	if val_max_min >= 100000:
		delta_y = 2000
	elif val_max_min >= 50000:
		delta_y = 2000
	elif val_max_min >= 10000:
		delta_y = 1000
	elif val_max_min >= 5000:
		delta_y = 200
	elif val_max_min >= 1000:
		delta_y = 100
	elif val_max_min >= 500:
		if general_max >=500: 
			delta_y = 80
		else:
			delta_y = 40
	elif val_max_min >= 100:
		if general_max >=500: 
			delta_y = 40
		else:
			delta_y = 20
	elif val_max_min >= 50:
		if general_max >=500: 
			delta_y = 100
		else:
			delta_y = 10
	elif val_max_min >= 10:
		if general_max >=500: 
			delta_y = 50
		else:
			delta_y = 2
	elif val_max_min >= 5:
		if general_max >=500: 
			delta_y = 50
		else:
			delta_y = 1
	elif val_max_min >= 1:
		if general_max >=500: 
			delta_y = 50
		else:
			delta_y = 0.2
	else:
		if general_max >=500: 
			delta_y = 50
		else:
			delta_y = 0.1

	if is_negative:
		delta_y = - delta_y

	return delta_y

def save_a_plot(fig, title, root='.'):
	try:
		title2 = title_replace(title)
		filefig = os.path.join(root, title2 + ".html")
		fig.write_html(filefig)

		filefig = filefig.replace(".html", ".png")
		fig.write_image(filefig)
	except:
		print("Could not save figure", title)

class BarGraphic:
	def __init__(self, dpi=300):
		self.dpi = dpi
		self.fontsize = 16- 2* int(round( (dpi-100)/ 100.,0))


	def gHist(self, lines, columns, numOfFig, havexSeq, xSeq, ySeq, yLine, parSDV, maxEntropy, par_xmin, par_xmax, par_yMin, par_ymax, par_xlabel, par_ylabel, par_title, par_width, par_color, wantTicks=True, seqSdv=None):
		plt.subplot(lines, columns, numOfFig)

		plt.subplots_adjust(left=0.10, right=0.95, bottom=0.05, top=0.95)

		wid = par_width / 35

		if (len(yLine) == 0):
			if (parSDV == 0):
				plt.bar(xSeq, ySeq, width=wid)
			else:
				plt.bar(xSeq, ySeq, yerr=parSDV, width=wid, error_kw=dict(elinewidth=wid / 3, color=par_color, ecolor='red'))
		else:
			plt.bar(xSeq, ySeq, yerr=parSDV, width=wid, error_kw=dict(elinewidth=wid / 3, color=par_color, ecolor='red'))
			plt.plot(xSeq, yLine, color='red')

			seqSup = []
			seqInf = []
			seqMean = []

			for _ in range(len(yLine)):
				seqSup.append(maxEntropy)
				seqMean.append(maxEntropy - parSDV)
				seqInf.append(maxEntropy - 2 * parSDV)

			plt.plot(xSeq, seqSup, color='black')
			plt.plot(xSeq, seqMean, color='red')
			plt.plot(xSeq, seqInf, color='black')

		if wantTicks and havexSeq:
			plt.xticks(xSeq)

		self.yMax = par_ymax
		self.yMin = par_yMin

		self.xMax = par_xmax
		self.xMin = par_xmin

		plt.xlim(par_xmin, par_xmax)
		plt.ylim(par_yMin, par_ymax)

		plt.title(par_title, fontsize=self.fontsize)
		plt.xlabel(par_xlabel, fontsize=self.fontsize)
		plt.ylabel(par_ylabel, fontsize=self.fontsize)


	# typeGraph par_width, yLine,parSDV,  maxEntropy,par_xmax
	def gBar(self, lines, columns, numOfFig, havexSeq, xSeq, ySeq, par_xmin, par_yMin, par_ymax, par_xlabel, par_ylabel, par_title, wantTicks=True, seqSdv=None):
		plt.subplot(lines, columns, numOfFig)

		plt.subplots_adjust(left=0.10, right=0.95, bottom=0.05, top=0.95)

		par_xmax = np.max(xSeq)
		plt.plot(xSeq, ySeq, 'r--', color='blue')
		if seqSdv:
			plt.errorbar(x=xSeq, y=ySeq, yerr=seqSdv, ecolor='red')

		if wantTicks and havexSeq:
			plt.xticks(xSeq)

		plt.xlim(par_xmin, par_xmax)
		plt.ylim(par_yMin, par_ymax)

		plt.title(par_title, fontsize=self.fontsize)
		plt.xlabel(par_xlabel, fontsize=self.fontsize)
		plt.ylabel(par_ylabel, fontsize=self.fontsize)


	def sameBar(self, lines, columns, numOfFig, xSeq, ySeq, linestyleCode, color):
		plt.subplot(lines, columns, numOfFig)

		yMax = np.max(ySeq)
		yMin = np.min(ySeq)

		if yMax > self.yMax:
			self.yMax = yMax

			if self.yMax >= 0:
				self.yMax *= 1.1
			else:
				self.yMax *= .9

		if yMin < self.yMin:
			self.yMin = yMin

			if self.yMin >= 0:
				self.yMin *= .9
			else:
				self.yMin *= 1.1

		plt.ylim(self.yMin, self.yMax)
		plt.plot(xSeq, ySeq, linestyleCode, color=color)

	def gNestHist(self, lines, columns, numOfImage, xSeq, ySeq, yLine, 
			   	  parError, par_xmin, par_xmax, par_yMin, par_ymax, par_xlabel, par_ylabel, par_title, par_width, par_color, 
				  meanF:float=None, stdF:float=None, medianF:float=None):
		ax = plt.subplot(lines, columns, numOfImage)

		if (len(yLine) == 0):
			ax.bar(xSeq, ySeq, yerr=parError, width=par_width, error_kw=dict(elinewidth=par_width/3, color=par_color, ecolor='red'))
		else:
			ax.bar(xSeq, ySeq, yerr=parError, width=par_width, error_kw=dict(elinewidth=par_width/3, color=par_color, ecolor='red'))

			plt.plot(xSeq, yLine, color='red')
			seqSup = []
			seqInf = []
			for i in range(len(yLine)):
				seqSup.append(yLine[i] + parError)
				seqInf.append(yLine[i] - parError)
			plt.plot(xSeq, seqSup, color='black')
			plt.plot(xSeq, seqInf, color='black')

		# vertical line from (70,100) to (70, 250)
		if meanF:
			plt.plot([meanF, meanF], [par_yMin, par_ymax], 'k-', lw=2, color='black')
			plt.plot([meanF+stdF, meanF+stdF], [par_yMin, par_ymax], '--', lw=2, color='red')
			plt.plot([meanF-stdF, meanF-stdF], [par_yMin, par_ymax], '--', lw=2, color='red')
			plt.plot([medianF, medianF], [par_yMin, par_ymax], 'k-', lw=2, color='yellow')

			ax.annotate(r'$1\sigma$', xy=((meanF+stdF)*1.05, (par_ymax-par_yMin)*.75), color='red')

		if ((par_yMin < 0.5) and (par_yMin > 0)):
			par_yMin = 0

		plt.xlim(par_xmin, par_xmax)
		plt.ylim(par_yMin, par_ymax)
		plt.title(par_title, fontsize=self.fontsize)

		plt.xlabel(par_xlabel, fontsize=self.fontsize)
		plt.ylabel(par_ylabel, fontsize=self.fontsize)

	def gNextBar(self, lines, columns, numOfImage, xSeq, ySeq, par_xmin, par_xmax, par_yMin, par_ymax, par_xlabel, par_ylabel, par_title, par_width, par_color, meanF=None, stdF=None, medianF=None):
		plt.subplot(lines, columns, numOfImage)

		plt.plot(xSeq, ySeq, 'r--', color='blue')

		if ((par_yMin < 0.5) and (par_yMin > 0)):
			par_yMin = 0

		plt.xlim(par_xmin, par_xmax)
		plt.ylim(par_yMin, par_ymax)
		plt.title(par_title, fontsize=self.fontsize)

		plt.xlabel(par_xlabel, fontsize=self.fontsize)
		plt.ylabel(par_ylabel, fontsize=self.fontsize)


	def printBar(self):
		plt.show()



class MultiLineWithTitle:
	def __init__(self, numLines, numCols, left, top, legColumns, legendTitle='', title='',dpi=120):

		self.myPlot = Plot()

		self.fig = plt.figure(1, dpi=dpi)

		self.numLines = numLines
		self.numCols = numCols
		self.left = left
		self.top = top
		self.legColumns = legColumns

		self.legendTitle = legendTitle
		''' 100 dpi = fs 10, 200 = fs 9, 300 = fs 8 '''
		self.fontsize = 10- 2* int(round( (dpi-100)/ 100.,0))

		self.fig.text(.5, .95, title, ha='center', fontsize=self.fontsize, color="blue")

		self.colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k',
					   'aqua', 'teal', 'lightslategray',
					   'maroon', 'chartreuse', 'fuchsia',
					   'navy', 'purple', 'silver', 'violet',
					   'b', 'g', 'r', 'c', 'm', 'y', 'k']


	def multilineGeneralWithTitle(self, iLoop, numFigure, xSeq, meanCurve3val, title, label, printLegend=True, xlabel='', ylabel=''):
		which = str(self.numLines) + str(self.numCols) + str(numFigure)
		ax = self.fig.add_subplot(which)


		'''
		wspace = 0.2   # the amount of width reserved for blank space between subplots
		hspace = 0.2   # the amount of height reserved for white space between subplots
		'''
		plt.subplots_adjust(left=0.1, right=0.90, bottom=0.15, top=0.85, wspace = 0.4)

		ax.set_title(title, fontsize=self.fontsize)
		ax.set_xlabel(xlabel, alpha=0.5, fontsize=self.fontsize-1)
		ax.set_ylabel(ylabel, alpha=0.5, fontsize=self.fontsize-1)

		#ax.tick_params(axis='both', which='major', labelsize=self.fontsize-2)
		#ax.tick_params(axis='both', which='minor', labelsize=self.fontsize-2)
		plt.tick_params(labelsize=self.fontsize-2)

		try:
			color1 = self.colors[iLoop]
		except:
			color1 = 'blue'

		value = []
		sdvCurve = []
		for i in range(len(meanCurve3val)):
			x = meanCurve3val[i]
			value.append(x[0])
			sdvCurve.append(x[2])

		if (color1):
			ax.plot(xSeq, value, color=color1, label=label)
			# http://www.thetechrepo.com/main-articles/469-how-to-change-line-properties-in-matplotlib-python
			ax.plot(xSeq, sdvCurve, color=color1, linestyle='-')

		else:
			ax.plot(xSeq, value, label=label)
			ax.plot(xSeq, sdvCurve, linestyle='dotted')

		if printLegend:
			ax.legend(loc="center left", bbox_to_anchor=[self.left, self.top],
					   ncol=self.legColumns, shadow=True, title=self.legendTitle)
			ax.get_legend().get_title().set_color("blue")



	def multilineGeneral_TR_WithTitle(self, iLoop, numFigure, xSeq, entropyCurve3Val, title, label, printLegend=True, xlabel='', ylabel=''):
		which = str(self.numLines) + str(self.numCols) + str(numFigure)
		ax = self.fig.add_subplot(which)

		ax.set_title(title)
		ax.set_xlabel(xlabel, alpha=0.5, fontsize=self.fontsize)
		ax.set_ylabel(ylabel, alpha=0.5, fontsize=self.fontsize)

		try:
			color1 = self.colors[iLoop]
		except:
			color1 = None

		meanCurveMax = []
		meanCurve = []

		# Tsallis and Renyi have errors
		for i in range(len(entropyCurve3Val)):
			x = entropyCurve3Val[i]
			print('x', x)
			value = x[0]
			sdv   = x[2]

			print('value', value, 'sdv', sdv)

			meanCurveMax.append(sdv)
			meanCurve.append(value)

		if (color1):
			ax.plot(xSeq, meanCurveMax, color=color1, linestyle='dotted')
			ax.plot(xSeq, meanCurve, color=color1, label=label)

		else:
			ax.plot(xSeq, meanCurveMax, linestyle='dotted')
			ax.plot(xSeq, meanCurve, label=label)

		if printLegend:
			ax.legend(loc="center left", bbox_to_anchor=[self.left, self.top],
					   ncol=self.legColumns, shadow=True, title=self.legendTitle)
			ax.get_legend().get_title().set_color("blue")


	def multilineGeneralWithTitle_and_HorizError(self, iLoop, numFigure, xSeq, meanCurve, error, title, label, printLegend=True, xlabel='', ylabel=''):

		which = str(self.numLines) + str(self.numCols) + str(numFigure)
		ax = self.fig.add_subplot(which)

		ax.set_title(title)
		ax.set_xlabel(xlabel, alpha=0.5, fontsize=self.fontsize)
		ax.set_ylabel(ylabel, alpha=0.5, fontsize=self.fontsize)

		try:
			color1 = self.colors[iLoop]
		except:
			color1 = None

		if (color1):
			ax.plot(xSeq, meanCurve, color=color1, label=label)
		else:
			ax.plot(xSeq, meanCurve, label=label)

		if printLegend:
			ax.legend(loc="center left", bbox_to_anchor=[self.left, self.top],
					   ncol=self.legColumns, shadow=True, title=self.legendTitle)
			ax.get_legend().get_title().set_color("blue")


		for i in range(len(xSeq)):
			x = xSeq[i]
			y = meanCurve[i]
			size = error[i]

			square = plt.Rectangle((x - size / 2, y - size / 2), size, size, facecolor='gray', color=color1)
			ax.add_patch(square)

		horZ = []
		print ('xSeq', np.max(xSeq), '  error', np.max(error))
		num = round((np.max(xSeq) + np.max(error)) / 10., 4)
		val = -num
		for i in range(11):
			horZ.append(val)
			val += num
		plt.xticks(horZ)
		ax.grid(True)


	def multilineGeneralWithTitle_and_VertError(self, iLoop, numFigure, xSeq, incX, meanCurve, error, title, label, printLegend=True, xlabel='', ylabel='', goOrigem=False):

		which = str(self.numLines) + str(self.numCols) + str(numFigure)
		ax = self.fig.add_subplot(which)

		ax.set_title(title)
		ax.set_xlabel(xlabel, alpha=0.5, fontsize=self.fontsize)
		ax.set_ylabel(ylabel, alpha=0.5, fontsize=self.fontsize)

		try:
			color1 = self.colors[iLoop]
		except:
			color1 = None

		if (color1):
			ax.plot(xSeq, meanCurve, color=color1, label=label)
		else:
			ax.plot(xSeq, meanCurve, label=label)

		if printLegend:
			ax.legend(loc="center left", bbox_to_anchor=[self.left, self.top],
					   ncol=self.legColumns, shadow=True, title=self.legendTitle)
			ax.get_legend().get_title().set_color("blue")


		horZ = []

		numPos = np.max(xSeq)
		numNeg = np.min(xSeq)

		if goOrigem:
			if numNeg < .5 and numNeg >= 0:
				numNeg = 0
			if numPos > -.5 and numPos < .5:
				numPos = 0

		tot = numPos - numNeg
		widthX = tot / 1000. # 0.2 / (iLoop+1)
		maxi = np.max(meanCurve)
		withY = maxi / 1000.

		for i in range(len(xSeq)):
			x = xSeq[i]
			y = meanCurve[i]
			size = error[i]

			square = plt.Rectangle((x - widthX, y - size / 2.), 2.*widthX, size, facecolor='gray', color=color1)
			ax.add_patch(square)

			square = plt.Rectangle((x - (5 * widthX), y - size / 2.), (10.*widthX), withY, facecolor='gray', color=color1)
			ax.add_patch(square)

			square = plt.Rectangle((x - (5 * widthX), y + size / 2.), (10.*widthX), withY, facecolor='gray', color=color1)
			ax.add_patch(square)

		if incX:
			inc = incX
		else:
			inc = abs(round(tot / 20., 2))

		val = numNeg
		for i in range(int(round(tot / inc, 0) + 1)):
			horZ.append(val)
			val += inc

		plt.xticks(horZ)
		ax.grid(True)



class MultiLine:
	def __init__(self, numLines, numCols, left, top, legColumns, legendTitle):
		self.numLines = numLines
		self.numCols = numCols
		self.left = left
		self.top = top
		self.legColumns = legColumns

		self.legendTitle = legendTitle

		self.colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k',
					   'aqua', 'teal', 'lightslategray',
					   'maroon', 'chartreuse', 'fuchsia',
					   'navy', 'purple', 'silver', 'violet', 'white']
		'''
		b : blue
		g : green
		r : red
		c : cyan
		m : magenta
		y : yellow
		k : black
		w : white

		http://www.tedmontgomery.com/tutorial/colors.html

		'''

	''' iLoop goes from 0 to numLines -1 '''
	def multiline(self, iLoop, numFigure, meanCurve, vc, prot, title):

		ax = plt.subplot(1, self.numCols, numFigure, title=title)

		try:
			color1 = self.colors[iLoop]
		except:
			color1 = None

		''' http://matplotlib.sourceforge.net/api/pyplot_api.html '''
		label = prot.name.replace('DROSOPHILA ', '') + '(ind=' + str(prot.numIndiv) + ' vc=' + str(round(vc * 100, 3)) + '%)'

		if (color1):
			meanCurveAux = []
			if (vc > 0):
				for i in range(len(meanCurve)):
					meanCurveAux.append(meanCurve[i] * (1 + vc))
				ax.plot(prot.alfa, meanCurveAux, color=color1, linestyle='dotted')

			ax.plot(prot.alfa, meanCurve, color=color1, label=label)

			if (vc > 0):
				for i in range(len(meanCurve)):
					meanCurveAux[i] = meanCurve[i] * (1 - vc)
				ax.plot(prot.alfa, meanCurveAux, color=color1, linestyle='dotted')
		else:
			meanCurveAux = []
			for i in range(len(meanCurve)):
				meanCurveAux.append(meanCurve[i] * (1 + vc))
			ax.plot(prot.alfa, meanCurve, label=prot.name.replace('DROSOPHILA ', ''))
			ax.plot(prot.alfa, meanCurveAux, linestyle='dotted')
			for i in range(len(meanCurve)):
				meanCurveAux[i] = meanCurve[i] * (1 - vc)
			ax.plot(prot.alfa, meanCurveAux, linestyle='dotted')

		if (numFigure == 1):
			ax.legend(loc="center left", bbox_to_anchor=[self.left, self.top],
					   ncol=self.legColumns, shadow=True, title=self.legendTitle)
			ax.get_legend().get_title().set_color("blue")

		ax.set_backgroundcolor(facecolor='b', alpha=0.5)
		plt.xticks(prot.alfa)

		# print 'iLoop', iLoop, 'numFigure', numFigure
		if (iLoop == (self.numLines - 1) and (numFigure == self.numCols)):
			plt.draw()
			plt.show()


	def multilineGeneral(self, iLoop, numFigure, xSeq, meanCurve, vc, title, label, printLegend=True, xlable='', ylable=''):
		ax = plt.subplot(self.numLines, self.numCols, numFigure, title=title)
		plt.ylabel(ylable)
		plt.xlabel(xlable)

		try:
			color1 = self.colors[iLoop]
		except:
			color1 = None

		if (numFigure == 1):
			print ('color', color1, 'iLoop', iLoop, 'figure', numFigure, 'title', title)

		''' http://matplotlib.sourceforge.net/api/pyplot_api.html '''
		if (color1):
			meanCurveAux = []
			if (vc > 0):
				for i in range(len(meanCurve)):
					meanCurveAux.append(meanCurve[i] * (1 + vc))
				ax.plot(xSeq, meanCurveAux, color=color1, linestyle='dotted')

			# print 'mandei bala', xSeq, meanCurve
			ax.plot(xSeq, meanCurve, color=color1, label=label)

			if (vc > 0):
				for i in range(len(meanCurve)):
					meanCurveAux[i] = meanCurve[i] * (1 - vc)
				ax.plot(xSeq, meanCurveAux, color=color1, linestyle='dotted')
		else:
			meanCurveAux = []
			if (vc > 0):
				for i in range(len(meanCurve)):
					meanCurveAux.append(meanCurve[i] * (1 + vc))
					ax.plot(xSeq, meanCurveAux, linestyle='dotted')

			ax.plot(xSeq, meanCurve, label=label)

			if (vc > 0):
				for i in range(len(meanCurve)):
					meanCurveAux[i] = meanCurve[i] * (1 - vc)
				ax.plot(xSeq, meanCurveAux, linestyle='dotted')

		if printLegend:
			ax.legend(loc="center left", bbox_to_anchor=[self.left, self.top],
					   ncol=self.legColumns, shadow=True, title=self.legendTitle)
			ax.get_legend().get_title().set_color("blue")


	def multiLinePrint(self):
		plt.draw()
		plt.show()

class HeatMap:
	# tirei o save e criei o print
	def __init__(self, is3D, lenSeq, mnat:str, tickWidth=50, rand_method_var:str='shuffle',
			     dna_prot='DNA', isLog:bool=False, dpi:int=300):

		self.dpi = dpi

		self.dna_prot = dna_prot
		self.mnat = mnat
		self.isLog = isLog
		self.rand_method_var = rand_method_var

		''' http://stackoverflow.com/questions/22408237/named-colors-in-matplotlib '''
		self.colorList = ['darkblue', 'lightblue', 'g', 'purple', 'r']
		self.markerList = ['.','1','2','D','o']

		plt.close("all")
		plt.clf()

		if is3D:
			self.fig = plt.figure(1, dpi=dpi)
		else:
			self.fig = plt.figure(1, figsize=(8, 8), dpi=dpi)

		self.lines = 1
		self.columns = 1

		self.fontsize = 24- 2* int(round( (dpi-100)/ 100.,0))

		self.left = 0.10
		self.bottom = 0.10
		self.width = 0.80
		self.height = .75

		self.tick = [i*tickWidth for i in range(lenSeq / tickWidth)]


	def colors(self):
		pass
		'''
		cnames = {
		'aliceblue':			'#F0F8FF',
		'antiquewhite':		 '#FAEBD7',
		'aqua':				 '#00FFFF',
		'aquamarine':		   '#7FFFD4',
		'azure':				'#F0FFFF',
		'beige':				'#F5F5DC',
		'bisque':			   '#FFE4C4',
		'black':				'#000000',
		'blanchedalmond':	   '#FFEBCD',
		'blue':				 '#0000FF',
		'blueviolet':		   '#8A2BE2',
		'brown':				'#A52A2A',
		'burlywood':			'#DEB887',
		'cadetblue':			'#5F9EA0',
		'chartreuse':		   '#7FFF00',
		'chocolate':			'#D2691E',
		'coral':				'#FF7F50',
		'cornflowerblue':	   '#6495ED',
		'cornsilk':			 '#FFF8DC',
		'crimson':			  '#DC143C',
		'cyan':				 '#00FFFF',
		'darkblue':			 '#00008B',
		'darkcyan':			 '#008B8B',
		'darkgoldenrod':		'#B8860B',
		'darkgray':			 '#A9A9A9',
		'darkgreen':			'#006400',
		'darkkhaki':			'#BDB76B',
		'darkmagenta':		  '#8B008B',
		'darkolivegreen':	   '#556B2F',
		'darkorange':		   '#FF8C00',
		'darkorchid':		   '#9932CC',
		'darkred':			  '#8B0000',
		'darksalmon':		   '#E9967A',
		'darkseagreen':		 '#8FBC8F',
		'darkslateblue':		'#483D8B',
		'darkslategray':		'#2F4F4F',
		'darkturquoise':		'#00CED1',
		'darkviolet':		   '#9400D3',
		'deeppink':			 '#FF1493',
		'deepskyblue':		  '#00BFFF',
		'dimgray':			  '#696969',
		'dodgerblue':		   '#1E90FF',
		'firebrick':			'#B22222',
		'floralwhite':		  '#FFFAF0',
		'forestgreen':		  '#228B22',
		'fuchsia':			  '#FF00FF',
		'gainsboro':			'#DCDCDC',
		'ghostwhite':		   '#F8F8FF',
		'gold':				 '#FFD700',
		'goldenrod':			'#DAA520',
		'gray':				 '#808080',
		'green':				'#008000',
		'greenyellow':		  '#ADFF2F',
		'honeydew':			 '#F0FFF0',
		'hotpink':			  '#FF69B4',
		'indianred':			'#CD5C5C',
		'indigo':			   '#4B0082',
		'ivory':				'#FFFFF0',
		'khaki':				'#F0E68C',
		'lavender':			 '#E6E6FA',
		'lavenderblush':		'#FFF0F5',
		'lawngreen':			'#7CFC00',
		'lemonchiffon':		 '#FFFACD',
		'lightblue':			'#ADD8E6',
		'lightcoral':		   '#F08080',
		'lightcyan':			'#E0FFFF',
		'lightgoldenrodyellow': '#FAFAD2',
		'lightgreen':		   '#90EE90',
		'lightgray':			'#D3D3D3',
		'lightpink':			'#FFB6C1',
		'lightsalmon':		  '#FFA07A',
		'lightseagreen':		'#20B2AA',
		'lightskyblue':		 '#87CEFA',
		'lightslategray':	   '#778899',
		'lightsteelblue':	   '#B0C4DE',
		'lightyellow':		  '#FFFFE0',
		'lime':				 '#00FF00',
		'limegreen':			'#32CD32',
		'linen':				'#FAF0E6',
		'magenta':			  '#FF00FF',
		'maroon':			   '#800000',
		'mediumaquamarine':	 '#66CDAA',
		'mediumblue':		   '#0000CD',
		'mediumorchid':		 '#BA55D3',
		'mediumpurple':		 '#9370DB',
		'mediumseagreen':	   '#3CB371',
		'mediumslateblue':	  '#7B68EE',
		'mediumspringgreen':	'#00FA9A',
		'mediumturquoise':	  '#48D1CC',
		'mediumvioletred':	  '#C71585',
		'midnightblue':		 '#191970',
		'mintcream':			'#F5FFFA',
		'mistyrose':			'#FFE4E1',
		'moccasin':			 '#FFE4B5',
		'navajowhite':		  '#FFDEAD',
		'navy':				 '#000080',
		'oldlace':			  '#FDF5E6',
		'olive':				'#808000',
		'olivedrab':			'#6B8E23',
		'orange':			   '#FFA500',
		'orangered':			'#FF4500',
		'orchid':			   '#DA70D6',
		'palegoldenrod':		'#EEE8AA',
		'palegreen':			'#98FB98',
		'paleturquoise':		'#AFEEEE',
		'palevioletred':		'#DB7093',
		'papayawhip':		   '#FFEFD5',
		'peachpuff':			'#FFDAB9',
		'peru':				 '#CD853F',
		'pink':				 '#FFC0CB',
		'plum':				 '#DDA0DD',
		'powderblue':		   '#B0E0E6',
		'purple':			   '#800080',
		'red':				  '#FF0000',
		'rosybrown':			'#BC8F8F',
		'royalblue':			'#4169E1',
		'saddlebrown':		  '#8B4513',
		'salmon':			   '#FA8072',
		'sandybrown':		   '#FAA460',
		'seagreen':			 '#2E8B57',
		'seashell':			 '#FFF5EE',
		'sienna':			   '#A0522D',
		'silver':			   '#C0C0C0',
		'skyblue':			  '#87CEEB',
		'slateblue':			'#6A5ACD',
		'slategray':			'#708090',
		'snow':				 '#FFFAFA',
		'springgreen':		  '#00FF7F',
		'steelblue':			'#4682B4',
		'tan':				  '#D2B48C',
		'teal':				 '#008080',
		'thistle':			  '#D8BFD8',
		'tomato':			   '#FF6347',
		'turquoise':			'#40E0D0',
		'violet':			   '#EE82EE',
		'wheat':				'#F5DEB3',
		'white':				'#FFFFFF',
		'whitesmoke':		   '#F5F5F5',
		'yellow':			   '#FFFF00',
		'yellowgreen':		  '#9ACD32'}
		'''


	def plotMI(self, title, seq, limSup):
		bottom = 0.08

		ax1 = plt.subplot(1, 2, 1)
		plt.subplots_adjust(left=0.05, right=0.97, bottom=bottom, top=0.82, wspace=0.3)

		if limSup:
			# , cmap=plt.cm.YlGnBu   cmap=plt.cm.get_cmap('OrRd'),
			im = ax1.matshow(seq, aspect='auto', origin='lower', cmap=plt.cm.get_cmap(self.heatmap_color), vmin=0, vmax=limSup)
		else:
			im = ax1.matshow(seq, aspect='auto', origin='lower') # , cmap=plt.cm.YlGnBu

		ax1.set_xticks(self.tick)
		ax1.set_yticks(self.tick)

		#plt.title("")
		self.fig.text(.5, .92, title, ha='center', fontsize=self.fontsize, color="black")

		ax1.set_xticklabels(self.tick,minor=False,fontsize=self.fontsize-4)
		ax1.set_yticklabels(self.tick,minor=False,fontsize=self.fontsize-4)

		if self.dna_prot=='DNA':
			label = 'bp-nuc'
			plt.xlabel(label, fontsize=self.fontsize)
			plt.ylabel(label, fontsize=self.fontsize)
		else:
			label = 'aa'
			plt.xlabel(label, fontsize=self.fontsize)
			plt.ylabel(label, fontsize=self.fontsize)


		width = 0.02
		height =.74

		axcolor = self.fig.add_axes([.46, bottom, width, height])
		cbar = plt.colorbar(im, cax=axcolor)
		cbar.ax.tick_params(labelsize=self.fontsize-2)
		plt.axhline()


	def densityHeatmapBar(self, seq, limSup, bins:int=20, par_color:str="blue"):
		if self.mnat:
			unit = 'mnat'
		else:
			unit = 'nat'

		ylabel = 'frequency (log10(#)'

		if self.isLog:
			xlabel = 'log(<VMI>) (%s)'%(unit)
		else:
			xlabel = '<VMI> (%s)'%(unit)

		meanY = np.mean(seq)
		stdY = np.sqrt(np.var(seq))
		medianY = np.median(seq)

		if self.rand_method_var == "shuffle" or self.rand_method_var == "random":
			roundVal = 3
		else:
			if self.mnat:
				roundVal = 2
			else:
				roundVal = 1

		strfloat='.'+str(roundVal)
		tit = 'Frequency Distribution \n mean=%xxxf(%xxxf); median=%xxxf %s'
		tit = tit.replace('xxx', strfloat)
		titleG = tit %(meanY, stdY, medianY, unit)


		ax = plt.subplot(1, 2, 2)
		(n, _, _) = ax.hist(seq,bins=bins,color=par_color,log=True)
		# (n, bins, patches) = ax.hist(seq,bins=bins,color=par_color)

		yMax = np.max(n)
		yMax2 = yMax*1.05


		plt.plot([meanY, meanY], [0, yMax2], 'k-', lw=2, color='black')
		plt.plot([meanY+stdY, meanY+stdY], [0, yMax2], '--', lw=2, color='red')
		plt.plot([meanY-stdY, meanY-stdY], [0, yMax2], '--', lw=2, color='red')
		plt.plot([medianY, medianY], [0, yMax2], 'k-', lw=2, color='yellow')

		xPosAnn = meanY+stdY*1.05

		ax.annotate(r'$1\sigma$', xy=(xPosAnn, yMax*.75), color='red') #oi

		plt.ylim(0, yMax2)
		plt.xlim(0, limSup)

		plt.title(titleG, fontsize=self.fontsize)

		plt.xlabel(xlabel, fontsize=self.fontsize)
		plt.ylabel(ylabel, fontsize=self.fontsize)
		plt.tick_params(labelsize=self.fontsize-2)

		return

	def plotMI_3D(self, title, seqX, seqY, seqZ, L, limSup, unit):
		"""
		Plot 3D mutual information (MI).

		Parameters
		----------
		title : str
			Plot title.
		seqX, seqY, seqZ : list-like
			Lists of x, y, z coordinate arrays.
			Expected shape: seqX[i], seqY[i], seqZ[i] for each group i.
		L : int or float
			Maximum x/y axis limit.
		limSup : int or float
			Maximum z axis limit.
		unit : str
			Unit for z-axis label.
		"""

		self.fig.clf()
		ax = self.fig.add_subplot(1, 1, 1, projection='3d')

		# Axis ticks
		ticks = np.arange(0, L + 1, 50)
		ax.set_xticks(ticks)
		ax.set_yticks(ticks)

		ax.tick_params(axis="x", labelrotation=90)
		ax.tick_params(axis="y", labelrotation=90)

		# Number of groups to plot
		n_groups = min(len(seqX), len(seqY), len(seqZ))

		for i in range(n_groups):
			color = self.colorList[i % len(self.colorList)]
			marker = self.markerList[i % len(self.markerList)]

			ax.scatter(
				seqX[i],
				seqY[i],
				seqZ[i],
				c=color,
				marker=marker,
				label=f"Group {i + 1}",
			)

		ax.set_xlabel("nuc", fontsize=12)
		ax.set_ylabel("nuc", fontsize=12)
		ax.set_zlabel(f"<VMI> ({unit})", fontsize=12)

		ax.set_xlim3d(0, L)
		ax.set_ylim3d(0, L)
		ax.set_zlim3d(0, limSup)

		ax.set_title(title, fontsize=14)

		# Optional legend
		if n_groups > 1:
			ax.legend()

		self.fig.tight_layout()

		return ax


	def printBar(self):
		plt.show()

	def savePlot(self, desk, pictureName):
		plt.savefig(pictureName, format=desk.imgType, dpi=desk.dpi)


class PrintDendogram:
	def __init__(self, seq, names):

		"""
		linkage(y, method='single', metric='euclidean'):

		Performs hierarchical/agglomerative clustering on the
		condensed distance ma  number of original observations paired
		in the distance matrix. The behavior of this function is very
		similar to the MATLAB(TM) linkage function.
		"""

		Z = linkage(seq)

		d = dendrogram(Z)

		d['leaves'] = names
		plt.show()



class HistogramGraphic:
	def __init__(self, fontsize=10):
		self.fontsize = fontsize

	def gHist(self, numFig, num, seq, par_xlabel, par_ylabel, par_title, par_color="blue", bins=20):
		plt.subplot(1, numFig, num)
		plt.hist(seq, bins=bins, color=par_color)
		'''
		(n, bins, patches) = plt.hist(parY)

		print 'n', n
		print 'bins', bins
		print 'patches', patches
		'''
		plt.ylabel(par_ylabel, fontsize=self.fontsize-2)
		plt.xlabel(par_xlabel, fontsize=self.fontsize-2)
		plt.title(par_title, fontsize=self.fontsize)

	def printBar(self):
		plt.show()


class PlotGraphic:
	def __init__(self, num, par_xseq, par_ySeq, par_xmin, par_xmax, par_yMin, par_ymax, par_xlabel, par_ylabel, par_title):
		xseq = []
		ySeq = []
		xMin = []
		xMax = []
		yMin = []
		yMax = []
		xlabel = []
		ylabel = []
		title = []
		figure = []

		self.xseq = xseq
		self.ySeq = ySeq
		self.xMin = xMin
		self.xMax = xMax
		self.yMin = yMin
		self.yMax = yMax
		self.xlabel = xlabel
		self.ylabel = ylabel
		self.title = title

		self.figure = figure

		self.init(num, par_xseq, par_ySeq, par_xmin, par_xmax, par_yMin, par_ymax, par_xlabel, par_ylabel, par_title)

		# print "inicializing ...", num, " '" + self.title[num] + "'"

	def init(self, num, par_xseq, par_ySeq, par_xmin, par_xmax, par_yMin, par_ymax, par_xlabel, par_ylabel, par_title):
		self.xseq.append(par_xseq)
		self.ySeq.append(par_ySeq)
		self.xMin.append(par_xmin)
		self.xMax.append(par_xmax)
		self.yMin.append(par_yMin)
		self.yMax.append(par_ymax)
		self.xlabel.append(par_xlabel)
		self.ylabel.append(par_ylabel)
		self.title.append(par_title)

		# self.figure.append(plt.figure(num))
		# plt.subplot(1, 2, num)

		# print "inicializing ...", num, " '" + self.title[num] + "'"



	def buildPlot(self, num, lineColor, backColor):
		fig = self.figure[num]

		ax = fig.add_subplot(111)
		# the histogram of the data
		ax.hist(self.xseq[num], color=lineColor, facecolor=backColor)


		# l = ax.plot(bincenters, y, 'r--', linewidth=1)
		ax.plot(self.xseq[num], self.ySeq[num], 'r--', linewidth=1)

		ax.set_xlabel(self.xlabel[num])
		ax.set_ylabel(self.ylabel[num])
		ax.set_title(self.title[num])

		# print 'x length ', len(self.xseq[num])
		ax.set_xlim(self.xMin[num], self.xMax[num])
		# print 'np.max value ',  np.max(self.ySeq)
		ax.set_ylim(self.yMin[num], self.yMax[num])

		ax.grid(True)



	def printPlot(self):
		plt.show()


class GeneralGraphic:
	def __init__(self, desk):
		self.myPlot = Plot()
		self.fontsize = 8
		self.fig = plt.figure(1, dpi=desk.dpi)
		plt.subplots_adjust(left=0.1, right=0.95, bottom=0.15, top=0.90)

	def histogram_H0_Ha(self, desk, listH0, listHa, label_random):
		maxi = np.max(listH0)
		maxi2 = np.max(listHa)

		if maxi2 > maxi:
			maxi = maxi2

		# fig = plt.figure(1, dpi=desk.dpi)
		ax = self.fig.add_subplot("122")

		ax.hist(listH0, bins=20, color="blue",alpha=0.3, label='H0: same species', edgecolor = "blue")
		ax.hist(listHa, bins=20, color="red", alpha=0.5, label='Ha: diff. species', edgecolor = "red")

		plt.ylabel("frequency", fontsize=self.fontsize-2)
		plt.xlabel("HMI in %s"%(desk.unit), fontsize=self.fontsize-2)
		plt.xlim(0, maxi)
		plt.legend(loc='upper right', fontsize=self.fontsize-5)
		if label_random != "":
			label_random = " " + label_random
		plt.title("JSD[HMI(%s%s)]%s"%(desk.minmax,desk.str_correction,label_random), fontsize=self.fontsize-1)
		plt.tick_params(labelsize=self.fontsize-2)


	def ROC_Curve(self, desk, sensitivityList, specificityList, label_random):
		y = np.array(sensitivityList)
		x = 1.-np.array(specificityList)

		ax = self.fig.add_subplot("121")

		ax.scatter(x,y, edgecolor = "none")
		plt.ylabel("TPR = sensitivity", fontsize=self.fontsize-2)
		plt.xlabel("FPR = 1-specificity", fontsize=self.fontsize-2)
		plt.xlim(0, 1)
		plt.ylim(0, 1)

		if label_random != "":
			label_random = " " + label_random

		plt.title("ROC for HMI(%s%s)%s"%(desk.minmax,desk.str_correction,label_random), fontsize=self.fontsize-1)

		ax.tick_params(labelsize=self.fontsize-2)




"""
Draws a 3D barchart
:param labels: Array_like of bar labels
:param z_data: Array_like of bar heights (data coords)
:param title: Chart title
:param z_title: Z-axis title
:param n_row: Number of x-rows
:param width: Chart width (px)
:param height: Chart height (px)
:param thikness: Bar thikness (0; 1)
:param colorscale: Barchart colorscale
:param **kwargs: Passed to Mesh3d()
:return: 3D barchart figure
"""
def barchart3d(labels, z_data, title, z_title,
			   n_row=0, width=900, height=900, thikness=0.7, colorscale='Viridis',
			   **kwargs):


	if n_row < 1:
		n_row = math.ceil(math.sqrt(len(z_data)))
	thikness *= 0.5
	ann = []

	fig = go.Figure()

	for iz, z_max in enumerate(z_data):
		x_cnt, y_cnt = iz % n_row, iz // n_row
		x_min, y_min = x_cnt - thikness, y_cnt - thikness
		x_max, y_max = x_cnt + thikness, y_cnt + thikness

		fig.add_trace(go.Mesh3d(
			x=[x_min, x_min, x_max, x_max, x_min, x_min, x_max, x_max],
			y=[y_min, y_max, y_max, y_min, y_min, y_max, y_max, y_min],
			z=[0, 0, 0, 0, z_max, z_max, z_max, z_max],
			alphahull=0,
			intensity=[0, 0, 0, 0, z_max, z_max, z_max, z_max],
			coloraxis='coloraxis',
			hoverinfo='skip',
			**kwargs))

		ann.append(dict(
			showarrow=False,
			x=x_cnt, y=y_cnt, z=z_max,
			text=f'<b>#{iz+1}</b>',
			font=dict(color='white', size=11),
			bgcolor='rgba(0, 0, 0, 0.3)',
			xanchor='center', yanchor='middle',
			hovertext=f'{z_max} {labels[iz]}'))

	# mesh3d doesn't currently support showLegend param, so
	# add invisible scatter3d with names to show legend
	for i, label in enumerate(labels):
		fig.add_trace(go.Scatter3d(
			x=[None], y=[None], z=[None],
			opacity=0,
			name=f'#{i+1} {label}'))

	fig.update_layout(
		width=width, height=height,
		title=title, title_x=0.5,
		scene=dict(
			xaxis=dict(showticklabels=False, title=''),
			yaxis=dict(showticklabels=False, title=''),
			zaxis=dict(title=z_title),
			annotations=ann),
		coloraxis=dict(
			colorscale=colorscale,
			colorbar=dict(
				title=dict(
					text=z_title,
					side='right'),
				xanchor='right', x=1.0,
				xpad=0,
				ticks='inside')),
		legend=dict(
			yanchor='top', y=1.0,
			xanchor='left', x=0.0,
			bgcolor='rgba(0, 0, 0, 0)',
			itemclick=False,
			itemdoubleclick=False),
		showlegend=True)
	return fig

'''
title = "Efeitos Adversos"
fields = ['AE_wo', 'AE_any', 'AE_many']
names  = ['Sem EA', 'Alguns EA', 'Muitos EA']
colors = ['cyan', 'orange', 'red']


title = "Efeitos Adversos Fracos"
fields = ['AE_systemic', 'local_pain']
names  = ['EA sistêmicos', 'Dor local']
colors = ['darkcyan', 'lightcoral'] # dimgray, dimgrey, 'orangered' 'lightpink'

title = "Efeitos Adversos Médios e Fortes"
fields = ['headache', 'fatigue', 'fever','AE_3degree']
names  = ['Dor de cabeça', 'Fatiga', 'Febre', 'EA fortes']
colors = ['darkcyan', 'rosybrown', 'coral', 'red'] # dimgray, dimgrey, 'orangered' 'lightpink'
from plotly.offline import init_notebook_mode, iplot
init_notebook_mode()
'''
def barplot_Efeitos_Adversos(dfe, title, root_result, fields, names, colors, fontsize=14, fontcolor='black', width=1000, height=1000, template='plotly_white'):
	traces = {}; count= 0
	for col in fields:
		vals = [None if pd.isnull(val) else val*100  for val in dfe[col]]
		traces[col] = go.Bar(x=dfe.company, y=vals, name=names[count],  marker_color=colors[count])
		count += 1

	data = [traces[key] for key in traces.keys()]

	layout = go.Layout(
		barmode = 'group',
		autosize=False,
		title=title,
		width=width,
		height=height,
		template=template,
		margin=dict( l=80, r=80, b=80, t=80, pad=4),
		font=dict(
			family="Arial, bold, monospace",
			size=fontsize,
			color=fontcolor
		),
		xaxis_title="companie studies",
		yaxis_title="% Efeitos Adversos",
		paper_bgcolor="whitesmoke",
		plot_bgcolor= "whitesmoke", # lightgrey ivory gainsboro whitesmoke lightsteelblue 'lightcyan' 'azure', white, lightgrey, snow ivory beige powderblue
		showlegend  = True
	)

	fig = go.Figure(data=data, layout=layout)

	title2 = title_replace(title) + ".html"
	filefig = os.path.join(root_result, title2)
	fig.write_html(filefig)

	filefig = filefig.replace(".html", ".png")
	fig.write_image(filefig)
	print("file(s) saved '%s'"%(filefig))

	return fig


def rgb_to_hex(r, g, b):
	return '#{:02x}{:02x}{:02x}'.format(int(r), int(g), int(b))

def calc_color(lfc, max_red = 4, min_blue = -4):

	perc = (lfc+max_red)/(max_red-min_blue)
	if perc > 1:
		perc = 1
	elif perc < 0:
		perc = 0

	pct_diff = 1.0 - perc
	red_color  = min(255, pct_diff*2 * 255)
	blue_color = min(255, perc*2 * 255)

	pcg_green = perc if perc <= 0.5 else 1-perc
	green_color = 255 * 2*pcg_green

	#print(perc, red_color, green_color, blue_color)
	return rgb_to_hex(red_color, green_color, blue_color)

