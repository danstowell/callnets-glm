# coding=UTF8
# just to plot the nonlinearity functions

import numpy as np
import matplotlib
#matplotlib.use('PDF') # http://www.astrobetter.com/plotting-to-a-file-in-python/
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams.update({'font.size': 6})




nlfuns = {
	'Rectifier':     ('r', '', lambda x: np.maximum(0, x)),
	'Softplus c=1':  ('c', '+', lambda x: np.log(1 + np.exp( 1 * x))/ 1),
	'Softplus c=10': ('b', 'o', lambda x: np.log(1 + np.exp(10 * x))/10),
	'Exponential':   ('g', 'x', lambda x: np.exp(x)),
}

evalpoints = np.linspace(-2, 2, 51)

params = {
   'axes.labelsize': 10,
   'text.fontsize': 10,
   'legend.fontsize': 10,
   'xtick.labelsize': 8,
   'ytick.labelsize': 8,
   'text.usetex': False,
   'figure.figsize': [5.5, 4],
   'lines.markersize': 4,
}
plt.rcParams.update(params)


pdf = PdfPages('pdf/plot_nonlins.pdf')
plt.figure(frameon=False)
plt.axes(frameon=0)
plt.axvline(0, color=[0.6]*3)
plt.axhline(0, color=[0.6]*3)
for nlname, (color, marker, nlfun) in nlfuns.items():
	plt.plot(evalpoints, map(nlfun, evalpoints), hold=True, label=nlname, color=color, marker=marker)
plt.title('Nonlinearities')
plt.xlim(-2, 2)
plt.ylim(-0.1, 2)
plt.legend(loc='upper left', frameon=False)
plt.xlabel('x')
plt.ylabel(u'σ(x)')
pdf.savefig(bbox_inches='tight')
plt.close()
pdf.close()

