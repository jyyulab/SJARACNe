##############################################################
## Python script to generate enhanced consensus network 	##
## (ncol) from a weak network (3col) by adding the gene 	##
## symbols and other attributes like: pearson, spearman,	##
## linear regression coefficient and p-value of the pairs	##
## in the network. This script will produce NDex and 		##
## Cytoscape visualization format and .xlsx format, as well.##
## This script submits the job to the computing clusters    ##
## at St. Jude Children's Research Hospital.                ##
## Author: Alireza Khatamian				    			##
## Email: akhatami@stjude.org				    			##
##############################################################
import os, sys, subprocess, sys, csv, re, igraph
from scipy import stats
import xlsxwriter

if len(sys.argv) < 3:
	print "Insufficient parameters!"
	print "Usage: python getenhancedconsensusnetwork.py <exp_file> <3col_network> <output_dir>"
	print "Example: python getenhancedconsensusnetwork.py input/infile.exp ***_3col.txt input/***/***.final"
	exit()

marker_set = {}
with open(sys.argv[1]) as exp:
	reader = csv.reader(exp, delimiter='\t')
	columns = reader.next()
	reader = csv.DictReader(exp, columns, delimiter='\t')
	for row in reader:
		samples = []
		for i in range(2, len(row)):
			samples.append(float(row[columns[i]]))
		marker_set[row[columns[0]]] = dict(symbol=row[columns[1]], samples=samples)

gene_network = {}
with open(sys.argv[2]) as net:
	reader = csv.reader(net, delimiter='\t')
	columns = reader.next()
	reader = csv.DictReader(net, columns, delimiter='\t')
	for row in reader:
		gene_network[row[columns[0]] + '--' + row[columns[1]]] = float(row[columns[2]])

for src_gene in marker_set.keys():
	for tar_gene in marker_set.keys():
		if src_gene != tar_gene:
			if src_gene + '--' + tar_gene in gene_network:
				src_marker = marker_set[src_gene]
				tar_marker = marker_set[tar_gene]
				src_symbol = src_marker['symbol']
				tar_symbol = tar_marker['symbol']
				src_samples = src_marker['samples']
				tar_samples = tar_marker['samples']
				rho, intercept, r, p, stderr = stats.linregress(src_samples, tar_samples)
				scc, sp = stats.spearmanr(src_samples, tar_samples)
				pcc, pp = stats.pearsonr(src_samples, tar_samples)

				mi = gene_network[src_gene + '--' + tar_gene]

				gene_network[src_gene + '--' + tar_gene] = dict(srcSymbol=src_symbol, tarSymbol=tar_symbol, mi=mi, rho=rho, intercept=intercept, r=r, p=p, stderr=stderr, scc=scc, sp=sp, pcc=pcc, pp=pp)

input_net_path_tokens = sys.argv[2].split('/')
input_net_name = input_net_path_tokens[len(input_net_path_tokens)-1]
input_net_name_tokens = re.split('_|\.', input_net_name)
out_name = '_'.join(input_net_name_tokens[0:2]) + '_ncol_' + '_'.join(input_net_name_tokens[3:len(input_net_name_tokens)-1])

out_subnet = None
subnet_list = []
if len(sys.argv) == 5:
	out_subnet = open(sys.argv[3] + 'sub_' + out_name + '.txt', 'w')
	with open(sys.argv[4]) as subnet_file:
		for _id in subnet_file:
			subnet_list.append(_id.split('\n')[0].strip())

out = open(sys.argv[3] + out_name + '.txt', 'w')
out_xlsx = xlsxwriter.Workbook(sys.argv[3] + out_name + '.xlsx')
out_xlsx_sheet = out_xlsx.add_worksheet()
header = ('source','target','source.symbol','target.symbol','MI','pearson','spearman','rho','p-value')
out.write('\t'.join(header) + '\n')
if out_subnet != None:
	out_subnet.write('\t'.join(header) + '\n')
row_index = 0
out_xlsx_sheet.write_row(row_index, 0, header)
row_index += 1
for edge in gene_network.keys():
	info = gene_network[edge]
	tokens = edge.split('--')
	src = tokens[0]
	tar = tokens[1]
	data = (src, tar, info['srcSymbol'], info['tarSymbol'], '{0:.4f}'.format(info['mi']), '{0:.4f}'.format(info['pcc']), '{0:.4f}'.format(info['scc']), '{0:.4f}'.format(info['rho']), '{0:.4f}'.format(info['p']))
	out.write('\t'.join(data) + '\n')
	out_xlsx_sheet.write_row(row_index, 0, data)
	if out_subnet != None:
		if info['srcSymbol'].strip() in subnet_list or info['tarSymbol'].strip() in subnet_list:
			out_subnet.write('\t'.join(data) + '\n')
	row_index += 1

out_xlsx.close()
out.close()
if out_subnet != None:
	out_subnet.close()

graph_name = '_'.join(input_net_name_tokens[0:2]) + '_' + '_'.join(input_net_name_tokens[3:len(input_net_name_tokens)-1])
graph = igraph.Graph.Read_Ncol(sys.argv[2])
graph.write_graphml(sys.argv[3] + graph_name + '.graphml')

'''
max_degree =  graph.maxdegree()
weights = graph.es['weight']
min_weight = min(weights)
max_weight = max(weights)
vertex_color = []
edge_color = []
for v in graph.vs(_degree=max_degree):
	print v

lgl_layout = igraph.Graph.layout_lgl(graph, maxiter=150, maxdelta=4000, area=4000*4000, coolexp=1.5, repulserad=4000*4000, root=graph.vs(_degree=max_degree)[0])
pal = igraph.RainbowPalette(n=256, s=0.8, v=0.8, alpha=0.8)
edge_color = [pal.get(int((e['weight'] - min_weight) / (max_weight - min_weight) * 255)) for e in graph.es]
vertex_color = [pal.get(int(v.degree() / max_degree * 255)) for v in graph.vs]
visual_style = {}
visual_style['vertex_size'] = [v.degree()/float(max_degree) * 10 for v in graph.vs]
visual_style['vertex_color'] = vertex_color
visual_style['edge_width'] = 1
visual_style['edge_color'] = edge_color
visual_style['layout'] = lgl_layout
visual_style['bbox'] = (4000,4000)
visual_style['margin'] = 20
visual_style['edge_arrow_width'] = 0.001
visual_style['edge_arrow_size'] = 0.001
igraph.plot(graph, sys.argv[3] + graph_name + '_lgl_1.png', **visual_style)


gfr_layout = igraph.Graph.layout_grid_fruchterman_reingold(graph, maxiter=150, maxdelta=len(graph.vs), area=len(graph.vs)*len(graph.vs)*len(graph.vs), coolexp=1.5, repulserad=len(graph.vs)*len(graph.vs)*len(graph.vs))

pal = igraph.RainbowPalette(n=256, s=0.8, v=0.8, alpha=0.8)
edge_color = [pal.get(int((e['weight'] - min_weight) / (max_weight - min_weight) * 255)) for e in graph.es]
vertex_color = [pal.get(int(v.degree() / max_degree * 255)) for v in graph.vs]
visual_style = {}
visual_style['vertex_size'] = [v.degree()/float(max_degree) * 100 for v in graph.vs]
visual_style['vertex_color'] = vertex_color
visual_style['edge_width'] = 1
visual_style['edge_color'] = edge_color
visual_style['layout'] = gfr_layout
visual_style['bbox'] = (4000,4000)
visual_style['margin'] = 20
visual_style['edge_arrow_width'] = 0.001
visual_style['edge_arrow_size'] = 0.001
igraph.plot(graph, sys.argv[3] + graph_name + '_gfr_4.png', **visual_style)
'''
