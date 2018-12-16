# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 16:28:03 2018

@author: DaVo
"""

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib
import numpy as np

from plotly.offline import download_plotlyjs, init_notebook_mode, iplot
init_notebook_mode(True)
import plotly.graph_objs as go
from IPython.display import Image


def showFigure(x):
	"""
	:param x: the name of the png figure to be plotted
	:type x: string    
	Plots the passed figure in the notebook.
	"""
		
			
	img = mpimg.imread(x)	
	plt.figure(figsize = (8, 6.4))
	imgplot = plt.imshow(img, aspect='auto')
