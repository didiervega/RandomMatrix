# RandomMatrix

# RandomMatrix: How to use


```python
%matplotlib inline
#from mlab.releases import latest_release as matlab
import matplotlib
import numpy as np
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

from plotly.offline import download_plotlyjs, init_notebook_mode, iplot
init_notebook_mode(True)
import plotly.graph_objs as go
from IPython.display import Image

%reload_ext autoreload
%autoreload 1
%aimport util

```


<script type="text/javascript">window.PlotlyConfig = {MathJaxConfig: 'local'};</script><script type="text/javascript">if (window.MathJax) {MathJax.Hub.Config({SVG: {font: "STIX-Web"}});}</script><script>requirejs.config({paths: { 'plotly': ['https://cdn.plot.ly/plotly-latest.min']},});if(!window._Plotly) {require(['plotly'],function(plotly) {window._Plotly=plotly;});}</script>


This is a notebook example of how to use the RandomMatrix code for multifractality in complex networks. The code is developed in MATLAB, and some toolboxes are necessary to run the commands. 



This is a Python interface for illustrative purpose of using and calling the commands
You run the MATLAB functions in bash/shell prompt as follow:


    matlab -nodisplay -nodesktop -nosplash -r "try; yourMATLABFunction(your parameters); catch; end; quit" 

Note: `matlab` is the environment variable in your system that calls the MATLAB program. Here, because of the Python interface, we prefix the command `%system`. But, do not need it if you run directly in the shell.



### Functions

The `getEtaUGraph` function shows the relative fluctations of the participation numbers ($\eta$) as a function of $\mu$ given the parameters

[PointsVar] = getEtaUGraph(Q,B,U,Ns,a)



```python
# Q: is for a specific q value
# B: the specific band
# U: the set of evaluation points for the participation numbers (x-axis)
# Ns: this is the vector of expoente network size, in the form 2^Ns[i]
# A: this is the sparcity parameter of the Power-Law Banded Random Matrix (PBRM) model.
```

##### Example of how to call the getEtaUGraph function


```python
%system matlab -nodisplay -nodesktop -nosplash -r "try; RandomMatrix.getEtaUGraph([2],[1],[0.6, 0.9, 1.0, 1.1],[6,7,8,9],1.0); catch; end; quit"

#just for visualizing the image here
Image(url = 'sBaFigure.png')


```




<img src="sBaFigure.png"/>



It presents the curves of $\eta$ vs $\mu$ for the dPBRM model with sparsity values $\alpha$ = 1 (no sparsity). 
The critical point is when the lines intercept.  
Remember, the `%system` command is only for the purpose of running here in this notebook.


The next function calculates the fractal dimensions of our studied model (dPBRM) and prints the figure, i.e, the signature of the multifractality of eigenfunctions of our network model.



[Data, eData] = getDqAlphaGraphMAT(Q,B,u,Ns,A)


```python
# Q: is a vector of all the points (y-axis Figure XX of the paper)
# B: this is a vector for all the band values to be calculated
# u: it is the critical point u_c
# Ns: this is the vector of expoente network size, in the form 2^Ns[i]
# A: this is the sparcity parameter of the Power-Law Banded Random Matrix (PBRM) model.
```

##### Example of how to call the getDqAlphaGraphMAT function


```python
%system matlab -nodisplay -nodesktop -nosplash -r "try; RandomMatrix.getDqAlphaGraphMAT([0.5,1,1.5,2,3,4,5],[1],1,[6,7,8,9],[0.3,1.0]); catch; end; quit"

#Remember, the `%system` command is only for the purpose of running here in this notebook.

#just for visualizing the image here
Image(url = 'DqVSq.png')


```




<img src="DqVSq.png"/>



It creates a figure that shows the multifractal dimensions $D_q$ vs $q$ for our model. (It will take a while)

![DqVSq.png](attachment:DqVSq.png)

