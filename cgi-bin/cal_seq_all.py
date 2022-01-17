#!/home/JFKlab/anaconda3/bin/python
##!/usr/bin/python

# Import sys module to output a graph in a standard output
import sys

# Import image downloading module
# import urllib
# urllib.urlretriver(url, image.png)

# Import modules for CGI handling
import cgitb
cgitb.enable()
import cgi

# Import plotting modules
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np; np.random.seed(1)
import matplotlib as mpl
import matplotlib.colors as colors
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Import image embedding module
from io import StringIO, BytesIO
import base64

# Create instance of FieldStorage
form = cgi.FieldStorage()

# Get data from fields
sequence = form.getvalue('sequence')
# sequence = form.getvalue('sequence').replace('\r','').replace('\n','').upper()

if sequence == None:
	print("Content-type:text/html\n")
	print("<html>")
	print("<head>")
	print("<title>Error</title>")
	print("</head>")
	print("<body>")
	print("<h2>Please input your sequence!</h2>")
	exit()
else:
	sequence = form.getvalue('sequence').replace(' ', '').replace('\r','').replace('\n','').upper()
	pass

# Check if the sequence is valid
aalist = "GALMFWKQESPVICYHRNDT"

for aa in sequence:
	if aa not in aalist:
		print("Content-type:text/html\n")
		print("<html>")
		print("<head>")
		print("<title>Error</title>")
		print("</head>")
		print("<body>")
		print("<h2>The input sequence contain invalid amino acid character out of this list</h2>")
		print("<h2>GALMFWKQESPVICYHRNDT</h2>")
		exit()
	else:
		pass

if len(sequence)<140:
	print("Content-type:text/html\n")
	print("<html>")
	print("<head>")
	print("<title>Warning: sequence is too short!</title>")
	print("</head>")
	print("<body>")
	print("<h2>Your sequence must be no shorter than 140!!!</title></h2>")
	print("</body>")
	print("</html>")

else:
	# Import pscore table calculation algorithm
	import pscore_table
		
	# Calculate Pscore
	table = pscore_table.table(sequence)
	overallscore = table.split("<br>")[0].replace('\t', ' ')
		
	# Import pscore plot calculation algorithm
	import pscore_plot

	# Calculate pscore for each seq_num and return two columns of data for plotting
	seq_n, pscore = pscore_plot.plot(sequence)
	x = np.array(seq_n)
	y = np.array(pscore)

	# plot data (sequence pscore heatmap + line plot)
	fig = plt.figure(figsize=[6.4, 4.8], dpi=600)
	
	# font definition
	# Set the font dictionaries (for plot title and axis titles)
	title_font = {'fontname':'Arial', 'size':'14', 'color':'black', 'weight':'normal',
        	      'verticalalignment':'bottom'} # Bottom vertical alignment for more space
	axis_font = {'fontname':'Arial', 'size':'12'}
	
	# cmap definition
	cmap = LinearSegmentedColormap.from_list('mycmap', ['#C2A5CF', '#E7D4E8', '#FFFFFF',  
						'#D9F0D3', '#A6DBA0', '#5AAE61','#1B7837'], N=500)
	cmap.set_over('#16632D')
	cmap.set_under('#9970AB')

	# subplot 1: heatmap of protein pscores over full sequence
	ax1 = plt.subplot(211)

	im1 = ax1.imshow(y[np.newaxis,:], cmap=cmap, 
					 vmin=-2, vmax=4, aspect='auto', interpolation="nearest", origin='lower'
					 ) # heatmap

	axins1 = inset_axes(ax1, width="3%", height="100%", loc="lower left",
						bbox_to_anchor=(1.02, 0., 1, 1),
						bbox_transform=ax1.transAxes,
						borderpad=0) # colorbar location setup
	ax1.set_title("protein " + overallscore, **title_font)
	ax1.autoscale_view()
	# ax1.set_xlim(int(min(x))-1, int(max(x)))
	# ax1.set_xlim(min(x), max(x))
	ax1.set_yticks([])
	ax1.set_yticklabels([], **axis_font)
	ax1.set_xticks([])
	ax1.set_xticklabels([], **axis_font)
	ax1.tick_params(axis = 'both', which = 'major')

	# subplot-2: line plot of protein pscores over full sequence
	ax2 = plt.subplot(212)
	ax2.set_xlim(int(min(x))-1, int(max(x)))
	# ax2.set_xlim(min(x), max(x))
	ax2.plot(x,y)
	ax2.axhline(y=0, color="grey", linestyle='-.', linewidth=1) # PScore=0 -> PDB average value
	ax2.axhline(y=4, color="red", linestyle='-.', alpha=0.5, linewidth=1) # PScore=4 -> confidence threshold
	xticks = np.arange(int(min(x))-1, int(max(x)), 10) # set up x-axis tick labels
	ax2.autoscale_view()
	# ax2.set_xticks(xticks)
	# ax2.set_xticklabels(xticks)
	ax2.tick_params(direction='in')
	ax2.set_xlabel('protein sequence', **axis_font)
	ax2.set_ylabel('PScore')
	ax2.text(int(max(x))+26, 0, 'PDB\naverage', horizontalalignment='center')
	ax2.text(int(max(x))+32, 4, 'Confidence\nthreshold', horizontalalignment='center')
#	ax2.text(int(min(x))+16, 2, "protein " + overallscore, horizontalalignment='left', bbox=dict(facecolor='lightgrey', alpha=0.7, color='lightgrey'))
	# plt.subplots_adjust(bottom=0.1, right=0.9, top=2)
	plt.subplots_adjust(hspace=.0, left=0.08, right=0.85, bottom=0.1, top=0.9)
	cbar = fig.colorbar(im1, extend='both', cax=axins1, )
	
	# print("Content-type: image/png\n")
	# plt.savefig(sys.stdout, format='png')
	# sio = StringIO()
	sio = BytesIO()

	plt.savefig(sio, format="jpeg")
	print("Content-type: text/html\n")
	
	print("<html>")
	print("<head>")
	print("<title>PScore table and plot</title>")
	print("</head>")
	print("<body>")
	print("<p>PScore plot:</p>")
	print("""<img src="data:image/jpeg;base64,{}"/ alt="pscore plot" width="640" height="480">""".format(base64.encodebytes(sio.getvalue()).decode()))	
	# print("""<img src="data:image/png;base64,%s" alt="pscore plot" width="500" height="600" />""") % sio.getvalue().encode("base64").strip()
	print("<br>")
	print("<p>%s</p>" % (table))
	print("</body></html>")
