import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from itertools import chain
from collections.abc import Iterable
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import colors
import math
import argparse

parser = argparse.ArgumentParser(description='Test for argparse')
parser.add_argument('--rnasamba', '-a', help='input rnasamba result',required=True)
parser.add_argument('--cpat', '-b', help='inpu cpat result',required=True)
parser.add_argument('--lncfinder', '-c', help='input lncfinder result',required=True)
parser.add_argument('--venn', '-d', help='output veen file',default="./veen.pdf")
parser.add_argument('--overlap', '-e', help='output overlap file',default="./overlap.txt")
args = parser.parse_args()

default_colors = [
	# r, g, b, a
	[92, 192, 98, 0.5],
	[90, 155, 212, 0.5],
	[246, 236, 86, 0.6],
	[241, 90, 96, 0.4],
	[255, 117, 0, 0.3],
	[82, 82, 190, 0.2],
]
default_colors = [
	[i[0] / 255.0, i[1] / 255.0, i[2] / 255.0, i[3]]
	for i in default_colors
]

def draw_ellipse(fig, ax, x, y, w, h, a, fillcolor):
	e = patches.Ellipse(
		xy=(x, y),
		width=w,
		height=h,
		angle=a,
		color=fillcolor)
	ax.add_patch(e)

def draw_triangle(fig, ax, x1, y1, x2, y2, x3, y3, fillcolor):
	xy = [
		(x1, y1),
		(x2, y2),
		(x3, y3),
	]
	polygon = patches.Polygon(
		xy=xy,
		closed=True,
		color=fillcolor)
	ax.add_patch(polygon)

def draw_text(fig, ax, x, y, text, color=[0, 0, 0, 1], fontsize=14, ha="center", va="center"):
	ax.text(
		x, y, text,
		horizontalalignment=ha,
		verticalalignment=va,
		fontsize=fontsize,
		color="black")

def draw_annotate(fig, ax, x, y, textx, texty, text, color=[0, 0, 0, 1], arrowcolor=[0, 0, 0, 0.3]):
	plt.annotate(
		text,
		xy=(x, y),
		xytext=(textx, texty),
		arrowprops=dict(color=arrowcolor, shrink=0, width=0.5, headwidth=8),
		fontsize=14,
		color=color,
		xycoords="data",
		textcoords="data",
		horizontalalignment='center',
		verticalalignment='center'
	)

def get_labels(data, fill=["number"]):
	"""
	get a dict of labels for groups in data
	@type data: list[Iterable]
	@rtype: dict[str, str]
	input
	  data: data to get label for
	  fill: ["number"|"logic"|"percent"]
	return
	  labels: a dict of labels for different sets
	example:
	In [12]: get_labels([range(10), range(5,15), range(3,8)], fill=["number"])
	Out[12]:
	{'001': '0',
	 '010': '5',
	 '011': '0',
	 '100': '3',
	 '101': '2',
	 '110': '2',
	 '111': '3'}
	"""
	N = len(data)
	sets_data = [set(data[i]) for i in range(N)]  # sets for separate groups
	s_all = set(chain(*data))					 # union of all sets
	# bin(3) --> '0b11', so bin(3).split('0b')[-1] will remove "0b"
	set_collections = {}
	for n in range(1, 2**N):
		key = bin(n).split('0b')[-1].zfill(N)
		value = s_all
		sets_for_intersection = [sets_data[i] for i in range(N) if  key[i] == '1']
		sets_for_difference = [sets_data[i] for i in range(N) if  key[i] == '0']
		for s in sets_for_intersection:
			value = value & s
		for s in sets_for_difference:
			value = value - s
		set_collections[key] = value
	labels = {k: "" for k in set_collections}
	if "logic" in fill:
		for k in set_collections:
			labels[k] = k + ": "
	if "number" in fill:
		for k in set_collections:
			labels[k] += str(len(set_collections[k]))
	if "percent" in fill:
		data_size = len(s_all)
		for k in set_collections:
			labels[k] += "(%.1f%%)" % (100.0 * len(set_collections[k]) / data_size)
	return labels

def venn3(labels, names=['A', 'B', 'C'], **options):
	"""
	plots a 3-set Venn diagram
	@type labels: dict[str, str]
	@type names: list[str]
	@rtype: (Figure, AxesSubplot)
	input
	  labels: a label dict where keys are identified via binary codes ('001', '010', '100', ...),
			  hence a valid set could look like: {'001': 'text 1', '010': 'text 2', '100': 'text 3', ...}.
			  unmentioned codes are considered as ''.
	  names:  group names
	  more:   colors, figsize, dpi, fontsize
	return
	  pyplot Figure and AxesSubplot object
	"""
	colors = options.get('colors', [default_colors[i] for i in range(3)])
	figsize = options.get('figsize', (9, 9))
	dpi = options.get('dpi', 96)
	fontsize = options.get('fontsize', 14)
	fig = plt.figure(0, figsize=figsize, dpi=dpi)
	ax = fig.add_subplot(111,aspect='equal')
	ax.set_axis_off()
	ax.set_ylim(bottom=0.0, top=1.0)
	ax.set_xlim(left=0.0, right=1.0)
	# body
	draw_ellipse(fig, ax, 0.333, 0.633, 0.5, 0.5, 0.0, colors[0])
	draw_ellipse(fig, ax, 0.666, 0.633, 0.5, 0.5, 0.0, colors[1])
	draw_ellipse(fig, ax, 0.500, 0.310, 0.5, 0.5, 0.0, colors[2])
	draw_text(fig, ax, 0.50, 0.27, labels.get('001', ''), fontsize=fontsize)
	draw_text(fig, ax, 0.73, 0.65, labels.get('010', ''), fontsize=fontsize)
	draw_text(fig, ax, 0.61, 0.46, labels.get('011', ''), fontsize=fontsize)
	draw_text(fig, ax, 0.27, 0.65, labels.get('100', ''), fontsize=fontsize)
	draw_text(fig, ax, 0.39, 0.46, labels.get('101', ''), fontsize=fontsize)
	draw_text(fig, ax, 0.50, 0.65, labels.get('110', ''), fontsize=fontsize)
	draw_text(fig, ax, 0.50, 0.51, labels.get('111', ''), fontsize=fontsize)
	# legend
	draw_text(fig, ax, 0.15, 0.87, names[0], colors[0], fontsize=fontsize, ha="right", va="bottom")
	draw_text(fig, ax, 0.85, 0.87, names[1], colors[1], fontsize=fontsize, ha="left", va="bottom")
	draw_text(fig, ax, 0.50, 0.02, names[2], colors[2], fontsize=fontsize, va="top")
	leg = ax.legend(names, loc='center left', bbox_to_anchor=(1.0, 0.5), fancybox=True)
	leg.get_frame().set_alpha(0.5)
	return fig, ax

def venn_plot(rnasamba,cpat,lncfinder,venn,overlap):
	#读取数据
	rnasamba= pd.read_table(rnasamba,header=None)
	rename=["id"]
	rnasamba.columns=rename
	rnasamba=set(rnasamba['id'])
	CPAT= pd.read_table(cpat,header=None)
	rename2=["id"]
	CPAT.columns=rename2
	CPAT=set(CPAT['id'])
	lncfinder= pd.read_table(lncfinder,header=None)
	lncfinder.columns=rename2
	lncfinder=set(lncfinder['id'])
	veen_list=[rnasamba,CPAT,lncfinder]
	labels =get_labels(veen_list, fill=['number'])
	fig, ax = venn3(labels, names=["lncScan_svm","CPAT","lncScan_XGB"],dpi=96)
	fig.savefig(venn,bbox_inches='tight',dpi=96)
	s1=rnasamba.intersection(CPAT)
	overlap_file=s1.intersection(lncfinder)
	overlap_file1=pd.Series(list(overlap_file))
	overlap_file2=pd.DataFrame(overlap_file1)
	overlap_file2.to_csv(overlap, sep='\t', index=False,header=False)

if __name__ == '__main__':
	try:
		venn_plot(args.rnasamba, args.cpat, args.lncfinder,args.venn,args.overlap)
	except Exception as e:
		print(e)
