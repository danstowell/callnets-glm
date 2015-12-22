#!/usr/bin/env python

import numpy as np
import csv, os
from scipy.stats import pearsonr

import matplotlib
#matplotlib.use('PDF') # http://www.astrobetter.com/plotting-to-a-file-in-python/
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages

from userconfig import *

##################################
# config

csvindir = os.path.expanduser('../code_GLM/outcsv')
plotdir = 'pdf'

# NB the order here is important: sessiongroups_within is day 0 day 1 just as is sessions_between. this is used when re-permuting the data
sessiongroups_within = [
['session2a1', 'session2a2', 'session2b1', 'session2b2'],
['session3a1', 'session3a2', 'session3b1', 'session3b2'],
]

sessions_between = ['session2full', 'session3full']

permutations = [
	('Individual', [0,1,2,3]),
	('Location',     [3,2,1,0]),
	('Null 1',   [2,3,0,1]),
	('Null 2',   [1,0,3,2]),
]

K = 4

colours = ['r', 'b', 'y', 'c']

##################################
# load each of the session CSVs
data = {}

def loadsomedata(sessname):
	if sessname in data:
		return
	with open('%s/zf4f_glm_stats_%s%s_2d.csv' % (csvindir, sessname, nlname), 'rb') as infp:
		rdr = csv.DictReader(infp)
		data_peakvals = np.zeros((K, K))
		for row in rdr:
			data_peakvals[int(row['frm'])-1, int(row['too'])-1] = float(row['peakval'])
		data[sessname] = data_peakvals

for agroup in sessiongroups_within:
	for sessname in agroup:
		loadsomedata(sessname)

for sessname in sessions_between:
	loadsomedata(sessname)

##################################
# utility func

def preprocess_strengthmatrix(matrix, whichsessnum, permmap=None):
	"permutation and other minor preproc of a matrix of strengths. REMOVES the self-self entries."
	# grab a version of the data after applying the permutation ONLY to the second session
	# NOTE: we skip the self-self
	if permmap!=None and whichsessnum==1:
		return np.array([[matrix[permmap[frm],permmap[too]] for too in range(K) if frm!=too] for frm in range(K)])
	else:
		return np.array([[matrix[        frm ,        too ] for too in range(K) if frm!=too] for frm in range(K)])

##################################
# analysis between the two full-sessions, conducted for each of the different permutations
results_perm = {}
results_perm_pval = {}
for permname, permmap in permutations:
	peaklists = {k:[] for k in sessions_between}
	for whichsessnum, whichsess in enumerate(sessions_between):
		sessx = preprocess_strengthmatrix(data[whichsess], whichsessnum, permmap)
		peaklists[whichsess].extend(sessx.flatten())

	correl, correl_pval = pearsonr((peaklists[sessions_between[0]]), (peaklists[sessions_between[1]]))
	print("Permutation %s: correlation coefficient is %g, standard p-value is %g" % (permname.ljust(9), correl, correl_pval))
	results_perm[permname] = correl
	results_perm_pval[permname] = correl_pval


plt.rcParams.update({'font.size': 8})
pdf = PdfPages('pdf/plot_analyseadj_between.pdf')
params = {
   'axes.labelsize': 8,
   'text.fontsize': 8,
   'legend.fontsize': 10,
   'xtick.labelsize': 10,
   'ytick.labelsize': 10,
   'text.usetex': False,
   'figure.figsize': [4.5, 4.5]
}
plt.rcParams.update(params)
fig = plt.figure(frameon=False)
plt.axes(frameon=False)
plotx, ploty = zip(*sorted(results_perm.items()))
plotpval = [results_perm_pval[x] for x in plotx]
plt.plot(ploty, 'x')
plt.xticks(range(len(plotx)), plotx)
# thanks https://github.com/jbmouret/matplotlib_for_papers#stars-statistical-significance
def stars(p):
	if p < 0.0001:
		return "****"
	elif (p < 0.001):
		return "***"
	elif (p < 0.01):
		return "**"
	elif (p < 0.05):
		return "*"
	else:
		return ""

for whichx, xval in enumerate(plotx):
	xstars = stars(plotpval[whichx])
	if xstars != "":
		#plt.annotate("", xy=(whichx-0.3, ploty[whichx]), xycoords='data',
		#	   xytext=(whichx+0.3, ploty[whichx]), textcoords='data',
		#	   arrowprops=dict(arrowstyle="-", ec='#aaaaaa',
		#		           connectionstyle="bar,fraction=0.2"))
		plt.text(whichx, ploty[whichx]+0.05, xstars,
		       horizontalalignment='center',
		       verticalalignment='center')

plt.xlim(-0.5, len(plotx)-0.5)
plt.ylim(ymax=1)
plt.axhline(0, color=[0.5]*3)
plt.axvline(plt.xlim()[0], color='k')
plt.xlabel("Data permutation")
plt.ylabel("Pearson correlation")
#plt.title('Predictability of inter-individual peak strengths,\nfrom one day to the next')
pdf.savefig()
plt.close()
pdf.close()





##################################
# analysis of 15-minute chunks, which does NOT relate to the permutation question
plotx = np.zeros((len(sessiongroups_within), K*K-K, len(sessiongroups_within[0])-1))
ploty = np.zeros((len(sessiongroups_within), K*K-K, len(sessiongroups_within[0])-1))
corrx = []
corry = []
# for each pair of adjacent segments
for whichgrp, agroup in enumerate(sessiongroups_within):
	for offset in range(len(agroup)-1):

		sessx = preprocess_strengthmatrix(data[agroup[offset  ]], whichgrp)
		sessy = preprocess_strengthmatrix(data[agroup[offset+1]], whichgrp)

		# list the adjacent segments' data as xy pairs
		corrx.extend(sessx.flatten())
		corry.extend(sessy.flatten())

		for too in range(K-1):
			for frm in range(K):
				plotx[whichgrp, frm+too*K, offset] = sessx[frm, too]
				ploty[whichgrp, frm+too*K, offset] = sessy[frm, too]

# find the pcorr of the collection
correl_15min, correl_15min_pval = pearsonr(corrx, corry)
print("Adjacent 15-minute segments: correlation coefficient is %g, standard p-value is %g" % (correl_15min, correl_15min_pval))

plt.rcParams.update({'font.size': 12})
pdf = PdfPages('pdf/plot_analyseadj_within.pdf')
fig = plt.figure(frameon=False)
plt.plot([np.min(ploty), np.max(ploty)], [np.min(ploty), np.max(ploty)], hold=True, linestyle='--', color=[0.7]*3)
print("shape of ploty: %s" % str(np.shape(ploty)))
for whichgrp, agroup in enumerate(sessiongroups_within):
	for too in range(K-1):
		for frm in range(K):
			frmtoo = frm+too*K
			colour = colours[frm]
			if too==0 and whichgrp==0:
				lbl = 'from %s' % (frm + 1)
			else:
				lbl = None
			plt.plot(plotx[whichgrp, frmtoo, :], ploty[whichgrp, frmtoo, :], hold=True, linestyle='-', marker='.', color=colour, label=lbl)
			plt.xlabel('Peak strengths in segment n')
			plt.ylabel('Peak strengths in segment n+1')
			plt.title('Predictability of inter-individual peak strengths,\nfrom one 15-minute segment to the next')
			plt.legend()

plt.text(np.max(ploty), np.min(ploty), "R=%.2g %s" % (correl_15min, stars(correl_15min_pval)), ha='right')

pdf.savefig()
plt.close()
pdf.close()


###############################################################################
# a sequential plot of the strengths, to show visually the continuity

chunkposses = []
for whichgrp, agroup in enumerate(sessiongroups_within):
	if whichgrp==0:
		chunkoffset = 0
	else:
		chunkoffset = chunkposses[-1][1]
	chunkposses.append((chunkoffset, chunkoffset + len(agroup)))

	sessx = np.array([preprocess_strengthmatrix(data[dataname], whichgrp) for dataname in agroup])
	print("shape of sessx is %s" % str(sessx.shape))

	if whichgrp==0:
		ploty = sessx
	else:
		ploty = np.vstack((ploty, sessx))

	print("shape of ploty is %s" % str(ploty.shape))

pdf = PdfPages('pdf/plot_analyseadj_15seq.pdf')
fig = plt.figure(frameon=False)
for whichchunk, (chunkstart, chunkend) in enumerate(chunkposses):
	for frm in range(K):
		colour = colours[frm]
		for tooish in range(K-1):
			if whichchunk==0 and tooish==0:
				lbl = 'from %s' % (frm + 1)
			else:
				lbl = None
			plotdata = ploty[chunkstart:chunkend, frm, tooish]
			plt.plot(range(chunkstart, chunkend), plotdata, hold=True, linestyle='-', marker='.', color=colour, label=lbl)

			# we also draw "connectors" between sessions
			if whichchunk != 0:
				plotdata = ploty[chunkstart-1:chunkstart+1, frm, tooish]
				plt.plot(range(chunkstart-1, chunkstart+1), plotdata, hold=True, linestyle=':', marker=None, color=colour, label=None)

plt.xlim(chunkposses[0][0]-0.5, chunkposses[-1][1]-0.5)
plt.xticks(range(ploty.shape[0]), [""] * ploty.shape[0])
plt.yticks(range(4))
plt.xlabel("Index of 15-minute session segment (sequential, across two days)")
plt.ylabel("Inter-individual peak strength")
plt.legend(loc='upper left', frameon=False)
plt.gca().set_frame_on(False)
plt.axhline(plt.ylim()[0], color='k', lw=1.5)
plt.axvline(plt.xlim()[0], color='k', lw=1.5)

print("presumed 3b1 data:")
print(ploty[6, :, :])

pdf.savefig()
plt.close()
pdf.close()

