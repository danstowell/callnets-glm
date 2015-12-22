#!/usr/bin/env python

# script to find the log-odds-ratios between model-fits that use exponential and that use softplus

import numpy as np
import csv, os

##################################
# config

csvindir = os.path.expanduser('../code_GLM/outcsv')

sessions = [
	'session2full', 'session3full',
]

nltypes = {
	'exp':'exp',
	'sof':'sof',
}


##################################

print("Log odds ratios (positive values mean evidence in favour of the exp model; negative values favour the softplus model)")

lls = {}
lors = {}

for sessname in sessions:

	lls[sessname] = {}

	for nltype, nlextension in nltypes.items():
		with open('%s/zf4f_glm_stats_%s%s_0d.csv' % (csvindir, sessname, nlextension), 'rb') as infp:
			rdr = csv.DictReader(infp)
			for row in rdr:
				lls[sessname][nltype] = 0 - float(row['neglogli'])

	lors[sessname] = lls[sessname]['exp'] - lls[sessname]['sof']
	print("Log odds ratio for %s: %g" % (sessname, lors[sessname]))

lorsdata = lors.values()
print("(Min, Median, Max) of log odds ratios: (%g, %g, %g)" % (np.min(lorsdata), np.median(lorsdata), np.max(lorsdata)))

print("de-logged, median odds ratio is %g" % np.exp(np.median(lorsdata)))

