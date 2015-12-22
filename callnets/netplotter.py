#!/usr/bin/env python

# code to create latex/tikz social network plots
# written by Dan Stowell 2015

import numpy as np
import subprocess, os


def runcmd(cmd, cwd=None, postoutput=True):
	print " ".join(cmd)
	if postoutput:
		stdout = None
	else:
		stdout = os.tmpfile()
	retcode = subprocess.call(cmd, cwd=cwd, stdout=stdout)
	if retcode != 0:
		raise RuntimeError("subprocess exit code was %i" % retcode)

def mat2str(mat):
	string = '['
	for row in mat:
		string += '[%s] ' % ["%.3g" % datum  for datum in row]
	string += ']'
	return string

def compose_tex_network_graph(title, matrixdata, K, normaliserval=None, nodelbls=None):
	"Returns a tex fragment for a network graph, to be embedded into a full latex file"
	rawmatrixdata = matrixdata

	if nodelbls==None:
		nodelbls = map(str, range(1, K+1))
	nodestyles = ["main node"] * K
	if K==4:
		# standard case
		trueK = K
	elif K<4:
		# upgrade to imitate K=4
		newmatrixdata = np.zeros((4,4))
		newmatrixdata[0:len(matrixdata), 0:len(matrixdata[0])] = matrixdata
		matrixdata = newmatrixdata
		trueK = K
		K = 4
		while len(nodelbls)<K:
			nodelbls.append("")
			nodestyles.append("draw=none")
	else:
		raise RuntimeError("compose_tex_network_graph() designed for K no more than 4 at the moment")
	print nodelbls
	matrixdata = np.array(matrixdata)
	if normaliserval==None:
		normaliserval = np.median(matrixdata[np.nonzero(np.maximum(0, matrixdata))])
	if np.isfinite(normaliserval):
		matrixdata = matrixdata / normaliserval
	loopers = ["loop above", "loop left", "loop below", "loop right"]
	arrows = []
	thicknesses = np.zeros((K,K))
	excitatory  = np.zeros((K,K), dtype=int)
	for frm in range(K):
		for too in range(K):
			thickness = np.abs(matrixdata[frm][too])
			excitatory[frm,too]  = matrixdata[frm][too] > 0
			if not np.isfinite(thickness):
				print("WARNING in compose_tex_network_graph(): matrixdata contains invalid data: %s" % str(matrixdata))
				thickness = 0
			if thickness > 0:
				thickness = pow(thickness, 1.5)
			thicknesses[frm,too] = thickness

	# new rescaling approach, hopefully applicable enough
	thicknesses = 7 * trueK * thicknesses / np.sum(thicknesses * excitatory)
	# since the inhibitors are not included in the rescaling, we hard-limit them to avoid plot craziness
	thicknesses = np.clip(thicknesses, 0, np.max(thicknesses * excitatory))

	for frm in range(K):
		for too in range(K):
			if thicknesses[frm,too] != 0:
				arrowcode = ['-|', '->'][excitatory[frm,too]]
				bendstring = ['bend right', loopers[frm]][frm==too]
				arrows.append(["  \\draw[%s,line width=%f pt, %s] (%i) to (%i);" % (arrowcode, thicknesses[frm,too], bendstring, frm+1, too+1)])

	texstring = """\\begin{minipage}[b]{3in}
\\begin{center}
%s

%% %s

%% %s

\\begin{tikzpicture}[->,>=stealth',shorten >=1pt,auto,node distance=3cm,
  thick,main node/.style={circle,fill=blue!20,draw,font=\\sffamily\\Large\\bfseries}]

  \\node[%s] (1) {%s};
  \\node[%s] (2) [below left of=1] {%s};
  \\node[%s] (3) [below right of=2] {%s};
  \\node[%s] (4) [below right of=1] {%s};

%s

 \\end{tikzpicture}
\\end{center}
\\end{minipage}""" % (title.replace('_', '\\_'), mat2str(rawmatrixdata), mat2str(matrixdata), nodestyles[0], nodelbls[0], nodestyles[1], nodelbls[1], nodestyles[2], nodelbls[2], nodestyles[3], nodelbls[3], "\n".join(["\n".join(row) for row in arrows]))

	return texstring

def plot_network_graphs(outstem, graphfragments, numcols=3):
	"Renders a set of network graphs into a PDF using tikz+latex"

	texstring = """\\documentclass{article}
\\usepackage[top=0.1in, bottom=0.1in, left=0.3in, right=0.3in, paperwidth=14in, paperheight=9in]{geometry}
\\usepackage{tikz}
\\usetikzlibrary{arrows}
\\begin{document}
\\begin{center}
"""

	counter = 0
	for frag in graphfragments:
		if frag!=None:
			texstring += frag
		counter += 1
		if (counter >= numcols) or (frag==None):
			texstring += "\n\n"
			counter = 0

	texstring += """

\\end{center}
\\end{document}
"""

	#print texstring
	with open("%s.tex" % outstem, 'wb') as outfp:
		outfp.write(texstring)
	runcmd(["pdflatex", "%s.tex" % os.path.basename(outstem)], cwd=os.path.dirname(outstem), postoutput=False)


#######################################################
if __name__=='__main__':
	frags = [
		compose_tex_network_graph("Example social network plot K=4", [[4,3,2,1], [1,1,1,1], [3,2,4,1], [3,2,4,1]], K=4),
		compose_tex_network_graph("Example social network plot K=3", [[4,3,1], [1,1,1], [3,2,4]], K=3),
	]
	plot_network_graphs('pdf/netplot_example', frags, numcols=3)

