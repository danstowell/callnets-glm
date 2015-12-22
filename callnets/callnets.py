#!/usr/bin/env python

# callnets: script to apply xcorr or code_GLM to animal call timing data
# (c) 2015 Dan Stowell

# NOTE: DATA FORMAT USED HERE IS (time, dur, id). consistent with the zcompiled etc

import numpy as np

import csv, copy, os, tempfile, shutil, re

import subprocess32 as subprocess

from operator import itemgetter

import scipy.stats

import matplotlib
#matplotlib.use('PDF') # http://www.astrobetter.com/plotting-to-a-file-in-python/
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams.update({'font.size': 6})

from netplotter import *
from userconfig import *

#############################################################################
# user config:

glmcodedir = os.path.expanduser('../code_GLM')  # where to find the Octave/Matlab code to run
csvstatssourcedir = os.path.expanduser('../code_GLM/outcsv')  # where to get CSV files to analyse, that have come from GLM analysis
rawdatasourcedir = os.path.expanduser('../data')  # where to get CSV files to analyse, the raw data files

# these options are not used by the GLM, but in some configurations of the cross-correlation
bindur = 0.005          # the size (in seconds) of the bins used for quantising/binarising data
maxlag_default = 1      # the size (in seconds) of the longest lag to consider/plot. don't make this too long or it'll take forever to run.
defaultcalldur = 0.1    # the duration of calls synthesised by our basic synthetic test cases

#############################################################################
# derived user config

defaultcalllen = int(defaultcalldur / bindur)
if defaultcalllen < 1:
	defaultcalllen = 1
	defaultcalldur = defaultcalllen * bindur
	print("Note: default call dur is smaller than 1 bin, so automatically adjusted to %g" % defaultcalldur)

#############################################################################

def generate_purerandom(k=3, rate=0.166667, dursecs=60):
	"Generate pure independent stationary poisson events - a null model as the simplest test case"
	allrate = rate * k
	# independent poisson events: generate independent exponential gap times, and independent identities, until the duration is complete
	# first we generate approximately the right amount
	gaps = np.random.exponential(1./allrate, dursecs * allrate * 1.1)
	times = np.cumsum(gaps)
	# then if we're still short we append
	while times[-1] < dursecs:
		times = np.hstack((times, times[-1]+np.random.exponential(1./allrate)))
	# now truncate to the correct place
	times = times[times < dursecs]
	# at last we can add IDs
	ids = np.random.randint(k, size=len(times))
	tstamps = [(times[i], defaultcalldur, ids[i]) for i in range(len(times))]
	return tstamps


def xcorr_pairwise(tstamps, k=3, maxlag=None, smoothdur=0):
	"""For each pair of IDs, turn the timestamp lists into quantised binary, and apply cross-correlation.
	Set smoothdur to 0 if you don't want to use smoothing, else to the smoothing scale desired, in secs.
		(Zero smoothing is sensible for durational data. For point data you do want some smoothing.)"""
	if maxlag==None:
		maxlag = maxlag_default
	# quantise the time data
	tstamps = [(int(np.round(time / bindur)), int(np.round(dur / bindur)), id) for time, dur, id in tstamps]
	maxlag  = int(np.round(maxlag / bindur))
	if smoothdur:
		smoothdur  = int(np.round(smoothdur / bindur))
		# create smoothing window
		# I'm using double the smoothdur in the window construction, so that it's approximately the std
		window = np.hanning(smoothdur * 4 + 1)
		if len(window)==0:
			window = np.array([1])
		window /= np.sum(window)
		print(" Smoothing window is of length %i" % len(window))

	# let's move the data so the first tstamp is at zero
	mintime = tstamps[0][0]
	tstamps = [(time-mintime, dur, id) for time, dur, id in tstamps]

	# create data structure
	# foreach individ, make a binary signal
	sigs = []
	for who in range(k):
		sig = np.zeros(tstamps[-1][0]+1) # zeros of full length used
		for time, dur, id in tstamps:
			if id==who:
				sig[time:time+dur] += 1
		if smoothdur:
			sig = np.convolve(sig, window, mode='same')
		sig -= np.mean(sig) # centring is useful because null is then at zero (since our raw data is nonneg this is not necc obv)
		l2norm = np.sqrt(np.sum(sig * sig))
		if l2norm != 0:
			sig /= l2norm # L2 normalisation - always gets our self-self-zero-lag peaks to 1
		sigs.append(sig)
	# foreach pair, do the xcorr
	print(" xcorring")
	xcorrs = {}
	for who1 in range(k):
		for who2 in range(who1, k):
			# here we do our own xcorr function, because np's method doesn't allow us to efficiently process long signals with short maxlag
			zeros = np.zeros(maxlag)

			sig1 = np.hstack((zeros, sigs[who1], zeros))
			sig2 = np.hstack((zeros, sigs[who2], zeros))

			truelen = len(sigs[0])
			xcorr = np.array([np.sum(sig2[offset:truelen+offset] * sig1[maxlag:truelen+maxlag]) for offset in range(maxlag*2+1)])

			if False:
				print(" Verifying our xcorr method gets the same result as numpy")

				# apply xcorr, truncate, and store result. the length of the array should always be odd, with the middle being zero lag.
				xcorr2 = np.correlate(sigs[who2], sigs[who1], mode='full')
				midpoint = (len(xcorr2)-1)/2
				assert maxlag < midpoint  # for sanity checking, not doing unsafe truncation. it also means we refuse to run on very short data.
				xcorr2 = xcorr2[midpoint-maxlag:midpoint+maxlag+1]

				diffmean = np.mean(xcorr - xcorr2)

				if diffmean > 1e-10:
					raise RuntimeError("diff mean: %g" % diffmean)

			xcorrs[(who1, who2)] = xcorr
	print(" done")

	return xcorrs

def plot_xcorrs(xcorrs, k, fig=None, plotmaxval=0.5, labels=None):
	"Nice standard xcorr multi-pane plot"
	if fig==None:
		fig = plt.figure(frameon=False)
	if labels==None:
		labels = map(str, range(1, k+1))
	xcorrlen = len(xcorrs[(0,0)])
	trumax = (xcorrlen-1)/2
	for who1 in range(k):
		for who2 in range(who1, k):
			plt.subplot(k,k,who2*k + who1 + 1)
			plt.plot(xcorrs[(who1, who2)])
			plt.ylim(-0.1, plotmaxval)
			ticks = plt.xticks()[0]
			plt.xticks(np.linspace(ticks[0], ticks[-1], 5), map(int, np.linspace(-trumax*bindur*1000, trumax*bindur*1000, 5)))
			if who1==0:
				plt.ylabel("to %s" % labels[who2])
			if who2==who1:
				plt.title("from %s" % labels[who1])
			if who2==k-1:
				plt.xlabel("lag (ms)")
	return fig

def xcorr_zap_selfpeaks(xcorrs, k):
	"Pre-process xcorrs data to remove the self-self peak: everything within the central two minima are set to the minima values"
	xcorrs = copy.deepcopy(xcorrs)
	centrepos = (len(xcorrs[(0,0)])-1)/2
	for who in range(k):
		anxcorrs = xcorrs[(who,who)]
		minpos_fwd = centrepos
		while minpos_fwd < len(anxcorrs)-1 and anxcorrs[minpos_fwd] >= anxcorrs[minpos_fwd+1]:
			minpos_fwd += 1
		minpos_bwd = centrepos
		while minpos_bwd > 0 and anxcorrs[minpos_bwd] >= anxcorrs[minpos_bwd-1]:
			minpos_bwd -= 1
		anxcorrs[centrepos:minpos_fwd] = anxcorrs[minpos_fwd]
		anxcorrs[minpos_bwd:centrepos] = anxcorrs[minpos_bwd]
	return xcorrs

def histogram_median(histo):
	"Find the median, given an equal-width-bin histogram. The return value is in bin indices, so you may want to convert to your preferred units."
	cdf = np.cumsum(np.abs(np.array(histo, dtype=float)))
	if cdf[-1]==0.:
		print("WARNING: histogram_median() was given an empty histogram")
		return 0
	cdf /= cdf[-1]
	# locate the position, and linear interpolate between adjacent bins
	i_hi = (cdf >= 0.5).nonzero()[0][0]
	i_lo = i_hi - 1
	c_hi = cdf[i_hi]
	c_lo = cdf[i_lo]
	weight_lo = c_hi - 0.5
	weight_hi = 0.5  - c_lo
	divisor = weight_hi + weight_lo
	if divisor != 0:
		i_interp = ((weight_lo * i_lo) + (weight_hi * i_hi)) / divisor
	else:
		i_interp = i_hi
	return i_interp

def find_xcorr_peaks(xcorrs, k):
	"Finds the peak positions and values from xcorr data. NOTE: we do extra processing to avoid the self-peak at zero in the self-self xcorrs."
	xcorrs_raw = xcorrs
	xcorrs = xcorr_zap_selfpeaks(xcorrs, k)  # for peak analysis, we zap the selfpeak since it's uninformative
	peakvals = np.zeros((k,k))
	peakposs = np.zeros((k,k))
	mednposs = np.zeros((k,k))
	centrepos = (len(xcorrs[(0,0)])-1)/2
	for who1 in range(k):
		for who2 in range(who1, k):
			# fwd
			half     = xcorrs[    (who1, who2)][centrepos+1:]
			half_raw = xcorrs_raw[(who1, who2)][centrepos+1:]
			apeakpos = np.argmax(half)
			peakposs[who1, who2] = apeakpos
			peakvals[who1, who2] = half[apeakpos]
			mednposs[who1, who2] = histogram_median(half_raw)
			# bwd
			half     = xcorrs[   (who1, who2)][centrepos-1::-1]
			half_raw = xcorrs_raw[(who1, who2)][centrepos-1::-1]
			apeakpos = np.argmax(half)
			peakposs[who2, who1] = apeakpos
			peakvals[who2, who1] = half[apeakpos]
			mednposs[who2, who1] = histogram_median(half_raw)

	peakposs *= bindur
	mednposs *= bindur

	return {'peakposs': peakposs, 'peakvals': peakvals, 'mednposs': mednposs}

def find_markov_tt(tstamps, k):
	"Simply retrieve the first-order Markov model transition table given some timestamp data"
	tt = np.zeros((k,k))
	for pos in range(len(tstamps)-1):
		tt[tstamps[pos][2]][tstamps[pos+1][2]] += 1
	for who in range(k):
		normer = np.sum(tt[who, :])
		if normer != 0:
			tt[who, :] /= normer
	return tt

def apply_glm_analysis(tstamps, k, runlabel='', workingfolder=None, regln=-1, resimuldur=0):
	"Invokes the matlab/octave analysis of data, via writing it out as CSV and running a subprocess"
	if workingfolder==None:
		# make temp folder
		tmpdir = tempfile.mkdtemp(prefix='callnetstooct')
	else:
		tmpdir = os.path.realpath(workingfolder)
	print("apply_glm_analysis(): using tmpdir %s" % tmpdir)

	tstamps.sort() # ensure time ordered
	firsttime = tstamps[0][0] # we'll subtract it, just because it makes the CSV data a bit clearer, and avoids some potential issues with v large timestamps

	# write data out as csv
	csvpath = "%s/indata.csv" % tmpdir;
	with open(csvpath, 'wb') as csvfp:
		for tstamp in tstamps:
			csvfp.write("%g,%g,%i\n" % (tstamp[0]-firsttime, tstamp[1], tstamp[2]))

	# invoke the analysis code
	cmds = ['octave', '-q', '--eval', "cd '%s'; pwd; setpaths; dofit_fromcsv_GLM_zf4f('%s', '%s', %s, %s, 0, 999999999999, %g, '%s', '%s', %i);"
			% (glmcodedir, csvpath, runlabel, k, str(range(1,k+1)), regln, [tmpdir, ''][workingfolder==None], tmpdir, resimuldur)]
	print cmds
	subprocess.check_call(cmds)

	numericaliser = {'runname':str, 'individ':int, 'frm':int, 'too':int, 'numcalls':int}
	numericalise = lambda row: {k:numericaliser.get(k, float)(v) for k, v in row.items()}

	# read the csvs back in that it's written
	with open('%s/zf4f_glm_stats__0d.csv' % tmpdir, 'rb') as infp:
		rdr = csv.DictReader(infp)
		for row in rdr:
			stats0d = numericalise(row)

	stats1d = [None for _ in range(k)]
	with open('%s/zf4f_glm_stats__1d.csv' % tmpdir, 'rb') as infp:
		rdr = csv.DictReader(infp)
		for row in rdr:
			row = numericalise(row)
			stats1d[row['individ']-1] = row

	stats2d = [[None for _ in range(k)] for _ in range(k)]
	with open('%s/zf4f_glm_stats__2d.csv' % tmpdir, 'rb') as infp:
		rdr = csv.DictReader(infp)
		for row in rdr:
			row = numericalise(row)
			stats2d[row['frm']-1][row['too']-1] = row

	kernels_disc = [[None for _ in range(k)] for _ in range(k)]
	with open('%s/zf4f_glm_stats__kernels_discret.csv' % tmpdir, 'rb') as infp:
		rdr = csv.reader(infp)
		for row in rdr:
			frm = int(row[1])
			too = int(row[2])
			kernels_disc[frm-1][too-1] = map(float, row[3:])

	with open('%s/zf4f_glm_stats__kernels_timeaxis.csv' % tmpdir, 'rb') as infp:
		rdr = csv.reader(infp)
		for row in rdr:
			timeaxis = map(float, row)

	# tidy up after ourselves; delete the temp stuff, but only if we haven't been told a special place to put things
	if workingfolder==None:
		shutil.rmtree(tmpdir)

	return (stats0d, stats1d, stats2d, kernels_disc, timeaxis)


def synthesise_abc(maxdursecs=6000):
	"Synthesises events dominated by an A->B->C influence pattern"

	#maxdursecs /= 10
	maxdursecs *= 16   # NB this gives me 2^16 data points max, which takes a while but shows a clear ABC convergence plot

	###########
	# params
	lambda_a = 0.1  # A events per second
	# here we define the two markov-renewal-type log-gaussian transitions A->B and B->C
	trans_probys = [0.9, 0.9]
	loggapmeans = np.log(np.array([1, 1]) * 0.75)
	loggapstdvs = [0.5, 0.5]


	events = [[0], [], []]


	# simulate A events by exponential gaps
	scale_a = 1. / lambda_a
	while events[0][-1] < maxdursecs:
		gapsize = np.random.exponential(scale_a)
		events[0].append(events[0][-1] + gapsize)


	# simulate follow-on events: foreach A, choose if B will occur, and then choose how far on
	for (frmwhich, toowhich) in [(0,1), (1,2)]:
		for frm in events[frmwhich]:
			if np.random.uniform() < trans_probys[frmwhich]:
				gapsize = np.exp(np.random.normal(loggapmeans[frmwhich], loggapstdvs[frmwhich]))
				events[toowhich].append(frm + gapsize)

	# now compile the tstamps all together
	tstamps = []
	for which, sublist in enumerate(events):
		tstamps.extend([(timepos, defaultcalldur, which) for timepos in sublist if timepos < maxdursecs])
	tstamps.sort()

	return tstamps

def datasize_convergence_analysis(tstamps, k, pdfoutstem, method='xcorr', ploteachstep=False, gtdata=None):
	toppow = int(np.floor(np.log2(len(tstamps))))
	print(" toplen is %i" % (2 ** toppow))
	results = {}
	usepows = np.arange(6, toppow + 0.0001, 1)

	if ploteachstep:
		pdf = PdfPages('%s_%s_eachstep.pdf' % (pdfoutstem, method))
	for whichcurpow, curpow in enumerate(usepows):
		curlen = int(2 ** curpow)
		curstart = (len(tstamps) - curlen) / 2
		curtstamps = tstamps[curstart:curstart+curlen]
		print(" %i/%i     (%i events)" % (whichcurpow+1, len(usepows), len(curtstamps)))

		if method=='xcorr':
			xcorrs = xcorr_pairwise(curtstamps, k=k, smoothdur=0)
			if whichcurpow==0:
				imgresult = np.zeros((len(xcorrs[(0,0)]), len(usepows), 3))
			imgresult[:, whichcurpow, 0] = xcorrs[(0,1)]
			imgresult[:, whichcurpow, 1] = xcorrs[(1,2)]
			imgresult[:, whichcurpow, 2] = xcorrs[(0,2)]
			if ploteachstep:
				fig = plt.figure(frameon=False)
				plot_xcorrs(xcorrs, k, fig)
				pdf.savefig()
				plt.close()
			results[curpow] = find_xcorr_peaks(xcorrs, k)
		elif method[:3]=='glm':
			if method[4:]=='reg':
				regln = -1
			elif method[4:]=='unreg':
				regln = 0
			else:
				raise ValueError("glm_ should be reg or unreg")
			(stats0d, stats1d, stats2d, kernels_disc, timeaxis) = apply_glm_analysis(curtstamps, k, regln=regln, workingfolder='working')
			timeaxis = np.array(timeaxis)
			if timeaxis[-1] > maxlag_default:
				timeaxis_cutoff = (timeaxis > maxlag_default).nonzero()[0][0]  # this allows us to limit the imgresult to the same timescale as in xcorr
			else:
				timeaxis_cutoff = len(timeaxis)
			# handle in the same way as the xcorr stuff above
			peakvals = np.array([[item['peakval'] for item in row] for row in stats2d])
			peakposs = np.array([[item['peakpos'] for item in row] for row in stats2d])
			mednposs = np.array([[histogram_median(kernels_disc[frm][too]) for too in range(3)] for frm in range(3)])

			# to compile the kernels down to 3-channel RGB, we concatenate the kernels of frmtoo and toofrm
			if whichcurpow==0:
				imgresult = np.zeros((timeaxis_cutoff * 2 - 1, len(usepows), 3))

			imgresult[:, whichcurpow, 0] = np.hstack((kernels_disc[1][0][1:timeaxis_cutoff][::-1], kernels_disc[0][1][:timeaxis_cutoff]))
			imgresult[:, whichcurpow, 1] = np.hstack((kernels_disc[2][1][1:timeaxis_cutoff][::-1], kernels_disc[1][2][:timeaxis_cutoff]))
			imgresult[:, whichcurpow, 2] = np.hstack((kernels_disc[2][0][1:timeaxis_cutoff][::-1], kernels_disc[0][2][:timeaxis_cutoff]))

			results[curpow] = {
				'peakposs': peakposs,
				'peakvals': peakvals,
				'mednposs': mednposs,
			}

	if ploteachstep:
		pdf.close()

	print("imgresult shape: %s" % str(imgresult.shape))

	# the imgresult now needs normalising and inverting
	imgresult_normsrc = imgresult[:, -4:, :] # we ignore the earliest entries because they're typically a bit bonkers
	if np.min(imgresult_normsrc) != np.max(imgresult_normsrc):
		imgresult = np.clip((imgresult - np.max(imgresult_normsrc)) / (np.min(imgresult_normsrc) - np.max(imgresult_normsrc)), 0, 1)
		if method[:3]!='glm':
			imgresult = pow(imgresult, 2)
	print("(min, median, max) of normalised imgresult is (%g, %g, %g)" % (np.min(imgresult), np.median(imgresult), np.max(imgresult)))

	# plot the convergence of the stats versus curpow
	pdf = PdfPages('%s_%s.pdf' % (pdfoutstem, method))

	xticks_pos = range(0, len(usepows)+1, 2)
	xticks_val = [2**int(val) for val in usepows[::2]]

	fig = plt.figure(frameon=False)
	# plot peakposs
	plt.subplot(4,1,1)
	for who1 in range(k):
		for who2 in range(k):
			if who1==who2: continue  # 'peak' usually meaningless for self-self since inhibition is the dominant effect
			plotcol = [None, (0.6)*3][who1==who2]
			plt.plot([results[curpow]['peakposs'][who1, who2] for curpow in usepows], hold=True)#, color=plotcol)
	plt.ylabel("peak positions")
	plt.xticks(xticks_pos, xticks_val)
	plt.xlim(-0.5, len(usepows)-0.5)
	plt.ylim(ymin=0)
	plt.title(method)
	# plot mednposs
	plt.subplot(4,1,2)
	for who1 in range(k):
		for who2 in range(k):
			if who1==who2: continue
			plotcol = [None, (0.6)*3][who1==who2]
			plt.plot([results[curpow]['mednposs'][who1, who2] for curpow in usepows], hold=True)#, color=plotcol)
		
	plt.ylabel("median positions")
	plt.xticks(xticks_pos, xticks_val)
	plt.xlim(-0.5, len(usepows)-0.5)
	plt.ylim(ymin=0)
	# plot peakvals
	plt.subplot(4,1,3)
	for who1 in range(k):
		for who2 in range(k):
			if who1==who2: continue  # 'peak' usually meaningless for self-self since inhibition is the dominant effect
			plotcol = [None, (0.6)*3][who1==who2]
			plt.plot([results[curpow]['peakvals'][who1, who2] for curpow in usepows], hold=True)#, color=plotcol)
	plt.ylabel("peak magnitudes")
	plt.xlabel("log2(num data points)")
	plt.xticks(xticks_pos, xticks_val)
	plt.xlim(-0.5, len(usepows)-0.5)
	plt.subplot(4,1,4)
	plt.imshow(imgresult, interpolation='nearest', origin='lower', aspect='auto')
	plt.xticks(xticks_pos, xticks_val)
	plt.ylabel("lag")
	plt.yticks(map(int, (imgresult.shape[0]-1) * np.array([0, 0.5, 1])), map(str, [maxlag_default, '0', maxlag_default]))
	plt.xlabel("Num calls analysed")
	plt.xlim(-0.5, len(usepows)-0.5)
	pdf.savefig()
	plt.close()
	pdf.close()


def simple_timeline_plot(tstamps, k, outpath, plttitle, indlabels):
	"simple plot of timeline of events"
	nsplit = 6
	maxplotdur = 300

	pdf = PdfPages(outpath)
	fig = plt.figure(frameon=False)

	lastpos  = tstamps[-1][0]
	# for plotting we want to limit to a max of maybe 300 sec? so...
	tstamps_toplot = [item for item in tstamps if item[0]>=lastpos-maxplotdur]
	firstpos = tstamps_toplot[0][0]
	chunksize = (lastpos - firstpos)/float(nsplit)
	for whichsub in range(nsplit):
		plt.subplot(nsplit, 1, whichsub+1)
		chunkstart = firstpos +  whichsub    * chunksize
		chunkend   = firstpos + (whichsub+1) * chunksize
		tstamps_toplot_subsec = [item for item in tstamps_toplot if item[0]>=chunkstart and item[0]<chunkend]
		plt.plot(np.array(map(itemgetter(0), tstamps_toplot_subsec)) - firstpos, map(itemgetter(2), tstamps_toplot_subsec), 'b.')
		plt.ylim(-1, k)
		plt.yticks(range(k), indlabels)
		plt.xlim(chunkstart-firstpos, chunkend-firstpos)
		if whichsub==0:
			plt.title(plttitle)
	pdf.savefig()
	plt.close()
	pdf.close()

################################################################################
# lovely aggregate kernel plots

def plot_aggregate_kernels(k, srcname, anlname, pairingslist, subsetlbl, subsetids, mapbacktoindivid, plotvariant=None, femalenesses=None):
	"""
	Makes a lovely aggregate kernel plot onto a single set of axes.

	k - num "channels" (usually individuals) to expect in the data
	srcname - a component of the input filename, also used in output filename
	anlname - a component of the input filename
	pairingslist - optional. a list of pairs, indicating which channels should be considered to be "partners" of each other. e.g. [(0,1), (2,3)] for two pairs.
	subsetlbl - optional. just goes on the end of the plot's title
	subsetids - optional. list of IDs to be included in the plot. else it includes all IDs.
	mapbacktoindivid - optional. if the data has been expanded so an individ occupies multiple channels, supply here a function that maps the IDs back to true individual IDs - this is so we can determine which things are to be labelled as self-self influences.
	"""

	pairings = np.zeros((k,k), dtype=int)
	if pairingslist != None:
		for inda, indb in pairingslist:
			pairings[inda, indb] = 1
			pairings[indb, inda] = 1

	if np.shape(srcname)==():  # single string, not list
		srcnames = [srcname]
	else:
		srcnames = srcname

	if plotvariant=='sexwise':
		kerneltypes = ['mf', 'ff', 'fm', 'mm']
	else:
		kerneltypes = ['ss', 'sp', 'so']

	kerneltype_cmap = {
		'ss': np.array([1.0, 0.3, 0.3]),
		'sp': np.array([0.6, 0.4, 0.0]),
		'so': np.array([0, 0, 1.0]),
		#
		'ff': np.array([0.5, 0.1, 0.1]),
		'mm': np.array([0.1, 0.1, 0.5]),
		'fm': np.array([0.9, 0.4, 0.2]),
		'mf': np.array([0.2, 0.6, 0.9]),
	}
	kerneltype_tmap = {
		'ss': 'Self-self' ,
		'sp': 'Self-partner',
		'so': 'Self-other',
		'ff': 'F->F',
		'mm': 'M->M',
		'fm': 'F->M',
		'mf': 'M->F',
	}

	kernels = {kerneltype:[] for kerneltype in kerneltypes}

	for whichsrcname, onesrcname in enumerate(srcnames):
		with open('%s/zf4f_glm_stats_%s%s_kernels_timeaxis.csv' % (csvstatssourcedir, onesrcname, anlname), 'rb') as infp:
			rdr = csv.reader(infp)
			for row in rdr:
				kernels_timeaxis = np.array(map(float, row))
		with open('%s/zf4f_glm_stats_%s%s_kernels_discret.csv' % (csvstatssourcedir, onesrcname, anlname), 'rb') as infp:
			rdr = csv.reader(infp)
			for whichrow, row in enumerate(rdr):
				frm = int(row[1]) - 1
				too = int(row[2]) - 1
				kern = map(float, row[3:])

				# here we check if our "subset" of individuals to plot is skipping this one. if not we add it
				if subsetids==None or ((frm, too) in subsetids):
					if mapbacktoindivid == None:
						isself = frm == too
					else:
						isself = mapbacktoindivid(frm) == mapbacktoindivid(too)
					ispartner = pairings[frm, too]

					if plotvariant=='sexwise':
						if ispartner or isself:  # we simply ignore other ones in this variant, for plot clarity
							frmf = femalenesses[frm]
							toof = femalenesses[too]
							[[kernels['mm'], kernels['mf']][toof],
							 [kernels['fm'], kernels['ff']][toof]][frmf].append(kern)
					else:
						[[kernels['so'], kernels['sp']][ispartner], kernels['ss']][isself].append(kern)

	plt.gca().set_frame_on(False)

	for layer in range(2):
		for kerneltype in kerneltypes:
			kernellist = np.array(kernels[kerneltype])
			titlelbl   = kerneltype_tmap[kerneltype]
			colourbase = kerneltype_cmap[kerneltype]

			if len(kernellist)==0:
				continue

			if plotvariant=='allcurves':
				if layer==0:
					plt.axhline(0, linestyle=':', color=[0.1]*3)
					for whichcurve, akernel in enumerate(kernellist):
						if whichcurve==0:
							label = "%s (N=%i)" % (titlelbl, len(kernellist))
						else:
							label = None
						plt.plot(kernels_timeaxis, akernel, linestyle='-', color=colourbase,
							hold=True, label=label, lw=0.05)
			else:

				k_median = np.median(kernellist, axis=0)
				k_upper  = np.percentile(kernellist, 97.5, axis=0)
				k_lower  = np.percentile(kernellist,  2.5, axis=0)

				if layer==0:
					plt.fill_between(kernels_timeaxis, k_lower, k_upper, hold=True, color=0.5+colourbase*0.5, facecolor=0.5+colourbase*0.5, alpha=0.4)
					plt.axhline(0, linestyle=':', color=[0.1]*3)
				elif layer==1:
					plt.plot(kernels_timeaxis, k_median, linestyle='-', color=colourbase, hold=True,
							label="%s (N=%i)" % (titlelbl, len(kernellist)))

	plt.xlabel("Time (s)")
	if anlname=='sof':
		sayanlname=''
	else:
		sayanlname=', %s' % anlname
	if np.shape(srcname)==():  # single string, not list
		saysrcname = srcname
	else:
		saysrcname = "%i datasets" % len(srcname)
	plt.title("Aggregated kernels (%s) %s" % (filename_to_titlename("%s%s" % (saysrcname, sayanlname)), [subsetlbl, ''][subsetlbl==None]))
	plt.legend(frameon=False)

	plt.xlim(0, 4.5)
	plt.ylim(-2, 2)
	plt.axhline(plt.ylim()[0], color='k') # just to add a visual line without having a full frame or grid


def plot_aggregate_kernels_multipdf(outlbl, k, subplotscheme, plotvariant=None, femalenesses=None, subplotlayout=None):
	"""
	outlbl - an identifier appended to the filename
	k - num individual channels in the data
	subplotscheme - a list containing, for each subplot to be plotted, the scheme of arguments to be used to construct it.
	"""
	if subplotlayout==None:
		subplotlayout = (2, 2)
	subplotsperpage = subplotlayout[0] * subplotlayout[1]

	params = {
	   'axes.labelsize': 10,
	   'text.fontsize': 10,
	   'legend.fontsize': 10,
	   'xtick.labelsize': 8,
	   'ytick.labelsize': 8,
	   'text.usetex': False,
	   'figure.figsize': [5.5 * subplotlayout[1], 4 * subplotlayout[0]]
	}
	plt.rcParams.update(params)

	plotvariantlbl = ['_%s' % plotvariant, ''][plotvariant==None]
	pdf = PdfPages('pdf/plot_kernels_ssso%s%s.pdf' % (outlbl, plotvariantlbl))
	fig = plt.figure(frameon=False)
	for whichsrc, (srcname, anlname, pairingslist, subsetlbl, subsetids, mapbacktoindivid) in enumerate(subplotscheme):
		subplotposition = (whichsrc % subplotsperpage) + 1
		if whichsrc!=0 and subplotposition==1:
			pdf.savefig(bbox_inches='tight')
			plt.close()
			fig = plt.figure(frameon=False)

		plt.subplot(subplotlayout[0], subplotlayout[1], subplotposition)
		plot_aggregate_kernels(k, srcname, anlname, pairingslist, subsetlbl, subsetids, mapbacktoindivid, plotvariant=plotvariant, femalenesses=femalenesses)

	pdf.savefig(bbox_inches='tight')
	plt.close()
	pdf.close()

def filename_to_titlename(astr):
	"Instead of plotting a raw filename in a page title, this transforms it to a prettier label which should match up with the paper text"
	# gill/gill_rawfile_4typeall -> Trial II Day 7
	# session2full -> Day 2
	for pattern, retval in [
		('session2full', 'Day 2'),
		('session3full', 'Day 3'),

		('gill/gill_rawfile_1(type|ind|clump)all', 'Trial II Day 1'),
		('gill/gill_rawfile_2(type|ind|clump)all', 'Trial II Day 2'),
		('gill/gill_rawfile_3(type|ind|clump)all', 'Trial II Day 3'),
		('gill/gill_rawfile_4(type|ind|clump)all', 'Trial II Day 7'),
		('gill/gill_rawfile_5(type|ind|clump)all', 'Trial II Day 11'),
		('gill/gill_rawfile_6(type|ind|clump)all', 'Trial II Day 18'),
		('gill/gill_rawfile_7(type|ind|clump)all', 'Trial II Day 20'),

		]:
		if re.match(pattern, astr):
			return retval

	return astr


################################################################################
if __name__ == '__main__':

	netfrags = [
#		compose_tex_network_graph("Example social network plot K=4", [[4,3,2,1], [1,1,1,1], [3,2,4,1], [3,2,4,1]], K=4),
#		compose_tex_network_graph("Example social network plot K=3", [[4,3,1], [1,1,1], [3,2,4]], K=3),
#		compose_tex_network_graph("Example social network plot K=3", [[4,3,1], [1,1,1], [-3,-2,-4]], K=3),
#		None,  # None is newline
	]
	netfrags_abc = []
	netfrags_zf4f = []

	if False:
		print("Running: simple check of what we get from xcorr calc, from limited finite data")
		k = 4
		data1 = [(0.0, 0.1, 0), (1.0, 0.1, 0), (4.0, 0.1, 0)]
		data2 = [(0.1, 0.1, 1), (1.1, 0.1, 1), (4.1, 0.1, 1)]
		data3 = [(0.1, 0.1, 2), (0.9, 0.1, 2), (4.0, 0.1, 2)]
		tstamps = data1 + data2 + data3
		tstamps.sort() # ensure time ordered
		xcorrs = xcorr_pairwise(tstamps, k=k, maxlag=2, smoothdur=0.01)
		# plot
		pdf = PdfPages('pdf/plot_callnets_simplecheck.pdf')
		fig = plt.figure(frameon=False)
		plot_xcorrs(xcorrs, k, fig)
		pdf.savefig()
		plt.close()
		pdf.close()

	if False:
		print("Running: simple check of what we get from xcorr calc, from purely independent random data")
		k = 5
		tstamps = generate_purerandom(k=k, dursecs=6000)
		xcorrs = xcorr_pairwise(tstamps, k=k, smoothdur=0)
		# plot
		pdf = PdfPages('pdf/plot_callnets_indepcheck.pdf')
		fig = plt.figure(frameon=False)
		plot_xcorrs(xcorrs, k, fig)
		pdf.savefig()
		plt.close()
		pdf.close()

	if True:
		print("Running: loading basis functions and plotting them")

		with open("%s/outcsv/zf4f_glm_stats_session2fullsof_kernels_timeaxis.csv" % glmcodedir, 'rb') as infp:
			rdr = csv.reader(infp)
			basisfuncs_t = np.array([map(float, row) for row in rdr]).flatten()
		with open("%s/outcsv/zf4f_glm_stats_session2fullsof_basis.csv" % glmcodedir, 'rb') as infp:
			rdr = csv.reader(infp)
			basisfuncs_r = np.array([map(float, row) for row in rdr]).T
		with open("%s/outcsv/zf4f_glm_stats_session2fullsof_basis_orth.csv" % glmcodedir, 'rb') as infp:
			rdr = csv.reader(infp)
			basisfuncs_o = np.array([map(float, row) for row in rdr])

		pdf = PdfPages('pdf/plot_basisfuncs_r.pdf')
		fig = plt.figure(frameon=False, figsize=[5.5, 3.5])
		plt.axes(frameon=0)
		plt.axhline(0, linestyle=':', color=[0.1]*3)
		plt.plot(basisfuncs_t, basisfuncs_r, '-')
		plt.xlim(0, 4.5)
		plt.xlabel("Lag (s)")
		plt.title("Basis functions used to compose kernels")
		pdf.savefig()
		plt.close()
		pdf.close()

		# here we create a random sample "s":
		basisfuncs_s = np.dot(np.random.normal(size=(10, 16)), basisfuncs_o).T

		pdf = PdfPages('pdf/plot_basisfuncs_s.pdf')
		fig = plt.figure(frameon=False, figsize=[5.5, 3.5])
		plt.axes(frameon=0)
		plt.grid(axis='y')
		plt.axhline(0, linestyle=':', color=[0.1]*3)
		plt.plot(basisfuncs_t, basisfuncs_s, '-')
		plt.ylim(-1, 1)
		plt.xlim(0, 4.5)
		plt.xlabel("Lag (s)")
		plt.title("Random kernel examples, drawn from symmetric Gaussian prior")
		pdf.savefig()
		plt.close()
		pdf.close()

	if True:
		print("Running: timeline plots for empirical data")
		for srclbl, k in [
			('session2a',       4),
			('session2b',       4),
			]:
			with open('%s/zf4f/zcompiled_%s.csv' % (rawdatasourcedir, srclbl), 'rb') as infp:
				rdr = csv.reader(infp)
				tstamps_resim = [(float(row[0]), float(row[1]), int(row[2])) for row in rdr]

			simple_timeline_plot(tstamps_resim, k=k, outpath='pdf/plot_timeline_%s.pdf' % srclbl, plttitle="Calling timeline (excerpt)", indlabels=map(str, range(1,1+k)))

	if True:
		print("Running: timeline plots for resimulated data")
		for srclbl, k in [
			('session2fullsof_resimulated',       4),
			('session2fullsof_resimulated_asolo', 1),
			]:
			with open('%s/zf4f_glm_stats_%s.csv' % (csvstatssourcedir, srclbl), 'rb') as infp:
				rdr = csv.reader(infp)
				tstamps_resim = [(float(row[0]), float(row[1]), int(row[2])) for row in rdr]

			simple_timeline_plot(tstamps_resim, k=k, outpath='pdf/plot_timeline_%s.pdf' % srclbl, plttitle="Resimulated timeline (excerpt)", indlabels=map(str, range(1,1+k)))

	if True:
		print("Running: loading zf4f results and plotting nets / aggregate kernel-plots")
		k = 4
		for srcname in [
			'session2full', 'session3full',
				]:
			if srcname==None:
				netfrags.append(None) # nl
				continue
			with open('%s/zf4f_glm_stats_%s%s_2d.csv' % (csvstatssourcedir, srcname, nlname), 'rb') as infp:
				rdr = csv.DictReader(infp)
				zf4f_glm_peakvals = np.zeros((k,k))
				zf4f_glm_peakposs = np.zeros((k,k))
				for row in rdr:
					frm = int(row['frm'])-1
					too = int(row['too'])-1
					zf4f_glm_peakvals[frm, too] = float(row['peakval'])
					zf4f_glm_peakposs[frm, too] = float(row['peakpos'])

				anetfrag = compose_tex_network_graph("GLM fitted network, %s" % filename_to_titlename(srcname), zf4f_glm_peakvals, K=k)
				netfrags.append(anetfrag)
				if srcname in ['session2full', 'session3full']:
					netfrags_zf4f.append(anetfrag)


		print("Running: loading zf4f results and plotting self-self and self-other plots")

		gill_femalenesses = {k:([1]*k+[0]*k)*4 for k in [1,3,5]}
		gilltype_labels = {
			5: ['Distance', 'Tets', 'Stacks', 'Kackles', 'Whines'],
			3: ['Distance', 'TetStack', 'KackleWhine'],
		}
		gillpairinglist = [(0,1), (2,3), (4,5), (6,7)]
		# here we construct the complicated list of pairings after gill's data is expanded ("inflated") to treat call-type separately
		gillinflations = [3, 5]
		gilltypepairings = {infl:[] for infl in gillinflations}
		gilltype_subsets = {infl:[] for infl in gillinflations}
		def gillmapindtobigid(infl, ind, calltype):
			return ind * infl + calltype
		# NOTE we can't create this from iteration because the variable-scope messes us up!
		gillmapbacktoind  = {
			3: lambda ind: (ind-(ind % 3))/3,
			5: lambda ind: (ind-(ind % 5))/5,
		}
		for infl in gillinflations:
			for calltype in range(infl):
				for calltype_too in range(infl):
					# we create a "subset" list that selects certain channels to be combined onto a single subplot
					gilltype_subsets[infl].append(("%s->%s" % (gilltype_labels[infl][calltype], gilltype_labels[infl][calltype_too]), [
							(frmid, tooid) for frmid in range(calltype, 8 * infl, infl) for tooid in range(calltype_too, 8 * infl, infl)]))

					# we create expanded "pairings" lists that tell us which channels are to be considered "within-pair" even though we've expanded each bird out to be multiple channels
					for frm, too in gillpairinglist:
						gilltypepairings[infl].append((gillmapindtobigid(infl, frm, calltype), gillmapindtobigid(infl, too, calltype_too)))

		for outlbl, k, subplotscheme, femalenesses, subplotlayout in [('', 4, [
						('session2full', 'sof', None, None, None, None),
						('session3full', 'sof', None, None, None, None),
					], None, (1,2)),
					('_resim', 4, [
						('resim/zf4f_glm_stats_session2fullsof_resim', 'sof', None, None, None, None),
						('resim/zf4f_glm_stats_session3fullsof_resim', 'sof', None, None, None, None),
						('resim/zf4f_glm_stats_session2fullsof_resimasolo', 'sof', None, None, None, None),
						('resim/zf4f_glm_stats_session3fullsof_resimasolo', 'sof', None, None, None, None),
					], None, (1,2)),
					] + [
					('_gill_ind', 8, [
						('gill/gill_rawfile_%iindall' % whichgillday, 'sof', gillpairinglist, None, None, None)
					for whichgillday in [1, 4, 5, 6]] # + [
						# this one should summarise them all onto one plot
						#(['gill/gill_rawfile_%iindall' % whichgillday for whichgillday in [1, 2, 3, 4, 5, 6, 7]], 'sof', gillpairinglist, None, None, None)]
					, gill_femalenesses[1], (4,1))
					] + [
					('_gill_%itype' % whichgillday, 40, [
						('gill/gill_rawfile_%itypeall' % whichgillday, 'sof', gilltypepairings[5], asubsetlbl, subsetids, gillmapbacktoind[5])
						for asubsetlbl, subsetids in gilltype_subsets[5]
					]
					, gill_femalenesses[5], (5,5)) for whichgillday in [1, 2, 3, 4, 5, 6, 7]
					] + [
					('_gill_aggtype', 40, [
						(['gill/gill_rawfile_%itypeall' % whichgillday for whichgillday in [1, 2, 3, 4, 5, 6, 7]], 'sof', gilltypepairings[5], asubsetlbl, subsetids, gillmapbacktoind[5])
						for asubsetlbl, subsetids in gilltype_subsets[5]
					]
					, gill_femalenesses[5], (5,5))
					]:

			print("plots underway for '%s'" % outlbl)
			plot_aggregate_kernels_multipdf(outlbl, k, subplotscheme, subplotlayout=subplotlayout)
			plot_aggregate_kernels_multipdf(outlbl, k, subplotscheme, plotvariant='allcurves', subplotlayout=subplotlayout)
			if femalenesses != None:
				plot_aggregate_kernels_multipdf(outlbl, k, subplotscheme, plotvariant='sexwise', femalenesses=femalenesses, subplotlayout=subplotlayout)

	if True:
		print("Creating ABC simulation")
		abcdata = synthesise_abc()
		print("   %i events" % len(abcdata))
		with open(os.path.expanduser('output/data_abc.csv'), 'wb') as csvfp:
			csvfp.write("time,dursecs,individ\n")
			for row in abcdata:
				csvfp.write("%f,%g,%i\n" % row)

	if True:
		print("Loading data for ABC simulation")
		with open(os.path.expanduser('output/data_abc.csv'), 'rb') as infp:
			rdr = csv.DictReader(infp)
			k = 3
			tstamps_abc = [(float(row['time']), float(row['dursecs']), int(row['individ'])) for row in rdr]

		simple_timeline_plot(tstamps_abc, k=3, outpath='pdf/plot_timeline_abc.pdf', plttitle="ABC timeline (excerpt)", indlabels=['A', 'B', 'C'])

		if False:
			# now we load the groundtruth parameters
			with open(os.path.expanduser('output/gtparams_mrp_abc_2d.csv'), 'rb') as csvfp:
				rdr = csv.DictReader(csvfp)
				k = 3
				getters = ['peakval', 'peakpos', 'mednpos']
				gtdata_abc = {getter:np.zeros((k, k)) for getter in getters}
				for row in rdr:
					frm = int(row['frm'])
					too = int(row['too'])
					for getter in getters:
						gtdata_abc[getter][frm, too] = float(row[getter])

				anetfrag = compose_tex_network_graph("Ground-truth influence strengths, ABC model", gtdata_abc['peakval'], K=k, nodelbls=['A', 'B', 'C'])
				netfrags.append(anetfrag)
				netfrags_abc.append(anetfrag)
		else:
			anetfrag = compose_tex_network_graph("Ground-truth influence strengths, ABC model", [[0, 1, 0], [0, 0, 1], [0, 0, 0]], K=k, nodelbls=['A', 'B', 'C'])
			netfrags.append(anetfrag)
			netfrags_abc.append(anetfrag)

	if True:
		print("Running: xcorr calc from ABC simulation")
		xcorrs = xcorr_pairwise(tstamps_abc, k=k, smoothdur=0.05, maxlag=2.5)
		# plot
		pdf = PdfPages('pdf/plot_callnets_abcsim.pdf')
		fig = plt.figure(frameon=False)
		plot_xcorrs(xcorrs, k, fig, labels=['A', 'B', 'C'])
		pdf.savefig()
		plt.close()
		pdf.close()

		thepeaks = find_xcorr_peaks(xcorrs, k)
		anetfrag = compose_tex_network_graph("xcorr peak strengths, ABC data", thepeaks['peakvals'], K=k, nodelbls=['A', 'B', 'C'])
		netfrags.append(anetfrag)
		netfrags_abc.append(anetfrag)

		tt_abc = find_markov_tt(tstamps_abc, k)
		anetfrag = compose_tex_network_graph("empirical Markov model, ABC data", tt_abc, K=k, nodelbls=['A', 'B', 'C'])
		netfrags.append(anetfrag)
		netfrags_abc.append(anetfrag)

	if True:
		print("Running: glm calc from ABC simulation")
		(stats0d, stats1d, stats2d, kernels_disc, timeaxis) = apply_glm_analysis(tstamps_abc, k, workingfolder='working', resimuldur=600)
		# plot in the same way as the xcorr stuff above

		peakvals = [[item['peakval'] for item in row] for row in stats2d]
		print("peakvals:")
		for row in peakvals:
			print(row)
		anetfrag = compose_tex_network_graph("GLM peak strengths, ABC data", peakvals, K=k, nodelbls=['A', 'B', 'C'])
		netfrags.append(anetfrag)
		netfrags_abc.append(anetfrag)


	if True:
		if False:
			print("Running: purerandom through the data-size-convergence analysis")
			tstamps = generate_purerandom(k=k, dursecs=60000)
			datasetlbl = 'indep'
			#gtdata = None
		else:
			print("Running: ABC through the data-size-convergence analysis")
			tstamps = tstamps_abc
			datasetlbl = 'abc'

		datasize_convergence_analysis(tstamps=tstamps, k=3, pdfoutstem='pdf/plot_callnets_sizetest_%s' % datasetlbl, ploteachstep=False, method='xcorr')
		datasize_convergence_analysis(tstamps=tstamps, k=3, pdfoutstem='pdf/plot_callnets_sizetest_%s' % datasetlbl, ploteachstep=False, method='glm_reg')
		datasize_convergence_analysis(tstamps=tstamps, k=3, pdfoutstem='pdf/plot_callnets_sizetest_%s' % datasetlbl, ploteachstep=False, method='glm_unreg')

	plot_network_graphs('pdf/netplot_callnets', netfrags, numcols=4)
	plot_network_graphs('pdf/netplot_callnets_zf4f', netfrags_zf4f, numcols=4)
	plot_network_graphs('pdf/netplot_callnets_abc', netfrags_abc, numcols=4)

