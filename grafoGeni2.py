#libries
import pickle
import numpy as np
import pandas as pd
import networkx as nx
from bisect import bisect_left
from bokeh.io import show, curdoc
from bokeh.models import Plot, Range1d, MultiLine, Circle, HoverTool, TapTool, BoxSelectTool
from bokeh.models.graphs import from_networkx, NodesAndLinkedEdges, EdgesAndLinkedNodes
from bokeh.palettes import Spectral4
#from bokeh.plotting import figure, gridplot
from bokeh.models import ColumnDataSource, LabelSet
from bokeh.layouts import widgetbox, column, row
from bokeh.models.widgets import Button, Slider, DataTable, TableColumn, TextInput

#variables
sorte=None
prev_chart=None
prev_chart2=None
new_network=None	
total_network=None
threshold=0.35
inputGenes=[]
gene_indice={}
pickleScore={}
done=[]
ed=None
geni=None
G=None
maxs=None
eligible_edges=None
scoreTot=None
geniuno=[]
genidue=[]
pvalues=[]
scoreees=[]
posss=[]
tablePval=None
pvalcalculated=False

print ('Lets begin')
with open('exe2_005filename.pickle', 'rb') as handle: #apri il file con i dati di geni con le frequenze
	print ('start')
	geni = pickle.load(handle)
tmp=[]	

#enter data
gene_input = TextInput(value="", title="Enter comma-separated genes: (eg. vit_12s0055g01000,vit_12s0055g01020)")
edges_input = TextInput(value="", title="Number of edges:")
#slider_edges=None
tr_input=TextInput(title="Enter threshold (we suggest 0.35):")
buttonM=Button(label='-', button_type='success')
buttonP=Button(label='+', button_type='success')

#get score data from file (random coples score)
def getScoreTot():
	global scoreTot
	with open('score10000Tot.pickle', 'rb') as handle:
		print ('scoreTot')
		scoreTot = pickle.load(handle)

#enter genes of the input local genes network
def input_genes(attrname, old, new): 
	global inputGenes
	global G
	global prev_chart
	global pickleScore
	inputGenes = gene_input.value.split(',')
	for gene in inputGenes: 
		done.append(gene) 
		pickleScore[gene]={} 
		for genedue in geni[gene]: 
			if gene!=genedue and genedue not in done: 
				a=geni[gene][genedue] 
				if genedue in geni.keys() and gene in geni[genedue].keys():
					b=geni[genedue][gene]
				else:
					b=0
			pickleScore[gene][genedue]={'weight':(a+b)/2} 

	G=nx.from_dict_of_dicts(pickleScore) 
	curdoc().remove_root(prev_chart) #remove input gene text box
	prev_chart=widgetbox(tr_input) 
	curdoc().add_root(prev_chart)#add input threshold text box

def calcPvalue():
	global pickleScore
	global scoreTot
	global geniuno
	global genidue
	global pvalues
	global tablePval
	global pvalcalculated
	print('getScoreTot')
	getScoreTot() 
	print('calcPvalue')
	scoreTot.sort()
	scoreDict={}
	tmp=0
	tmpp=[]
	for i in scoreTot: 
		scoreDict[round(i, 4)]=tmp 
		tmpp.append(round(i, 4))
		tmp=tmp+1
	tmpp.sort()
	a=np.array(tmpp) 
	print('score')
	print(type(scoreTot))
	print(type(scoreDict))

	for gene in pickleScore: 
		for genedue in pickleScore[gene]:
			geniuno.append(gene)
			genidue.append(genedue)
			score=pickleScore[gene][genedue]['weight']
			
			if score==1.0:
				print('aaaaaaaaaaaaaaaaaaaaaa')
				pvalues.append(0.0)
				scoreees.append(score)
			else:
				score=round(score, 4)
				if score in scoreDict:
					print('z')
					indice=scoreDict[score] 
				else:
					print('v') 
					indice=np.argmin(np.abs(a-score)) 
				percentuale=float(float(indice)/49995000)
				pvalue=round(1.0-percentuale, 16) #calculate pvalue
				pvalues.append(pvalue) 
				scoreees.append(score) 
	#generate the data frame (table)
	datafr=dict(geneUno = geniuno, geneDue = genidue, scoRe = scoreees, pvalue = pvalues) #creo la tabella dei pvalue
	sourcePval = ColumnDataSource(datafr)
	columnsPval = [
				TableColumn(field="geneUno", title="Gene 1"),
				TableColumn(field="geneDue", title="Gene 2"),
				TableColumn(field="scoRe", title="Score"),
				TableColumn(field="pvalue", title="Pvalue")
				]
	tablePval = DataTable(source=sourcePval, columns=columnsPval, width=700, height=280)
	pvalcalculated=True

#input the threshold, number of edges and buttons
def input_tr(attrname, old, new):
	global geni
	global gene_indice
	global done
	global eligible_edges
	global edges_input
	global prev_chart
	global ed
	print('tr')
	threshold=float(tr_input.value)
	print (threshold)
	eligible_edges = [(from_node,to_node, edge_attributes) for from_node,to_node,edge_attributes in G.edges(data=True) if edge_attributes['weight'] > threshold]
	ed=row(widgetbox(edges_input), buttonM, buttonP)

	new_root=column(widgetbox(tr_input), ed) 
	curdoc().remove_root(prev_chart)
	prev_chart=new_root
	curdoc().add_root(prev_chart) 
	print('in')
	buttonM.on_click(minus)
	buttonP.on_click(plus)
	edges_input.on_change("value", max_edges)

#min button
def minus():
	global edges_input
	#slider_edges.on_change("value", max_edges)#dividi la max_edges per poter richiamare la seconda parte
	if int(edges_input.value)>0:
		edges_input.value=str(int(edges_input.value)-1)

#plus button
def plus():
	global edges_input
	edges_input.value=str(int(edges_input.value)+1)

#generate graph
def max_edges(attrname, old, new):
	global prev_chart
	global new_network
	global tmp
	global eligible_edges
	global tablePval
	global edges_input
	global sorte
	global pvalcalculated
	new_network = nx.DiGraph()
	maxs = edges_input.value
	print('edge')
	total_network= nx.DiGraph()
	sorte=sorted(eligible_edges, reverse=True, key=lambda edge:edge[2]['weight'])
	z=0
	tmp=[]
	for from_node,to_node,edge_attributes in sorte:
			if z<int(maxs):
				tmp.append((from_node, to_node))
				z=z+1
	print(maxs)
	new_network.add_edges_from(tmp)
	pos=nx.spring_layout(new_network)
	total_network.add_edges_from(eligible_edges)
	posTot=nx.spring_layout(total_network)
	labels_dict={}
	i=1
	for j in pos.keys():
		labels_dict[j]=i
		i=i+1
	labels_dict_tot={}
	i=1
	for j in posTot.keys():
		labels_dict_tot[j]=i
		i=i+1
	df=[]
	for i in posTot.keys():
		df.append({'x': labels_dict_tot[i], 'y': i})
	#----
	plot = Plot(plot_width=600, plot_height=500,
	        x_range=Range1d(-1.1,1.1), y_range=Range1d(-1.1,1.1))
	plot.title.text = "Gene network"
	plot.add_tools(HoverTool(tooltips=None), TapTool(), BoxSelectTool())
	graph_renderer = from_networkx(new_network, nx.circular_layout, scale=1, center=(0,0))
	graph_renderer.node_renderer.glyph = Circle(size=15, fill_color=Spectral4[0])
	graph_renderer.node_renderer.selection_glyph = Circle(size=15, fill_color=Spectral4[2])
	graph_renderer.node_renderer.hover_glyph = Circle(size=15, fill_color=Spectral4[1])
	graph_renderer.node_renderer.glyph.properties_with_values()
	graph_renderer.edge_renderer.glyph = MultiLine(line_color="#CCCCCC", line_alpha=0.8, line_width=5)
	graph_renderer.edge_renderer.selection_glyph = MultiLine(line_color=Spectral4[2], line_width=5)
	graph_renderer.edge_renderer.hover_glyph = MultiLine(line_color=Spectral4[1], line_width=5)
	graph_renderer.selection_policy = NodesAndLinkedEdges()
	graph_renderer.inspection_policy = EdgesAndLinkedNodes()
	x, y = zip(*graph_renderer.layout_provider.graph_layout.values())
	plot.renderers.append(graph_renderer)
	sourceLab = ColumnDataSource({'x': x, 'y': y,
	                           'geneL': [labels_dict[i] for i in pos.keys()]})
	labels = LabelSet(x='x', y='y', text='geneL', source=sourceLab,
	                  background_fill_color='white')
	sourceDf=ColumnDataSource(data=pd.DataFrame(df))
	plot.renderers.append(labels)

	columns = [
	        TableColumn(field="x", title="Label"),
	        TableColumn(field="y", title="Gene")
	    ]
	data_table = DataTable(source=sourceDf, columns=columns, width=400, height=280)

	curdoc().remove_root(prev_chart)
	if pvalcalculated==False:
		calcPvalue()
	new_chart=column(ed, row(plot, widgetbox(data_table)), widgetbox(tablePval))
	curdoc().add_root(new_chart)
	prev_chart=new_chart

gene_input.on_change("value", input_genes)
tr_input.on_change("value", input_tr)
prev_chart=widgetbox(gene_input)
curdoc().add_root(prev_chart)
curdoc().title = "grafoGeni"