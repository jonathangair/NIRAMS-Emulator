import sys
sys.path.insert(1,'/Library/Python/2.7/site-packages')
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import george
import scipy.optimize as op
import emcee
import corner as triangle
import h5py
import math
import time
from progress.bar import Bar

from george.kernels import ExpSquaredKernel

if ((len(sys.argv) != 3)):
    print("Arguments: 1. site ID")
    print("           2. number of posterior samples to use (max 180000)")
    exit()
site=int(sys.argv[1])
numpost=int(sys.argv[2])

#datafile=r'NIRAMS_Emulator_Debug_GWBGroup%i.h5' % (site)
#datafile=r'OldRes/NIRAMS_Emulator_GWBGroup%i.h5' % (site)
datafile=r'NIRAMS_Emulator_GWBGroup%i.h5' % (site)
hf=h5py.File(datafile,'r')
e = 'calibration' in hf
if e is False:
	print("No calibration data found, exiting!")
	hf.close()
	exit()
g2=hf.get('cell list')
celllist=np.array(g2.get('cells in group'))
#for i in range(celllist.shape[0]):
#	print i," ",celllist[i][0]," ",celllist[i][1]
goodcols=np.array(g2.get('cells with good NIRAMS data'))
g3=hf.get('emulator')
kertype=np.array(g3.get('Kernel type'))
if (kertype != 'ExpSquared'):
	print("Unrecognized kernel type, exiting!")
	exit()
kerpar=np.array(g3.get('Kernel parameters'))
groupav=np.array(g3.get('group average'))
x=np.array(g3.get('x'))
y=np.array(g3.get('y'))
yerr=np.array(g3.get('yerr'))
kernel=100.0*ExpSquaredKernel([1.08,0.0015,0.00005,0.15],ndim=4)
gp=george.GP(kernel)
gp.kernel.set_parameter_vector(kerpar)
#print x
#print gp.kernel.get_parameter_vector()
#exit()
ct=0
while (~(goodcols[ct])):
	ct+=1
gp.compute(x, yerr[:,ct])
mu, cov = gp.predict(groupav,x)
#print np.linalg.inv(cov)
g4=hf.get('calibration')
obserror=np.array(g4.get('observation error'))
postsamps=np.array(g4.get('posterior samples'))
postsamps=postsamps[72000:,:]
numfullpost=postsamps.shape[0]
#numpost=10000
idx=np.random.randint(numfullpost,size=numpost)
postsamps=postsamps[idx,:]
npost=np.array(g4.get('production samples'))
nwalk=np.array(g4.get('nwalkers'))
groupavprob13=0.0
groupavprob37=0.0
groupavprob50=0.0
grouppostav=0.0
prob13=np.zeros((y.T).shape[0])
prob37=np.zeros((y.T).shape[0])
prob50=np.zeros((y.T).shape[0])
cellpostavs=np.zeros((y.T).shape[0])

bar = Bar('Processing', max=numpost, suffix='%(index)d/%(max)d - %(percent).1f%% - %(eta)ds')
for p in postsamps:
	#mu, cov = gp.predict(groupav,[np.array(p)])
	muarr, covarr = gp.predict(groupav,[np.array(p),[0.225,0.2,0.02,1.1]])
	mu=muarr[0]
	cov=np.diag(covarr)[0]
	modcov=cov+obserror*np.identity(1)
	#print mu," ",cov," ",modcov
	cov=modcov
	if (cov > 0.0):
		std=np.sqrt(cov)
	else:
		print("Error, cov = ",cov," p = ",p)
	std=1.0
	grouppostav+=mu
	thisprob=0.5*(1.0-math.erf((13.-mu)/std))
	groupavprob13+=thisprob
	thisprob=0.5*(1.0-math.erf((37.-mu)/std))
	groupavprob37+=thisprob
	thisprob=0.5*(1.0-math.erf((50.-mu)/std))
	groupavprob50+=thisprob
	for i in range((y.T).shape[0]):
		if (goodcols[i]):
			col=(y.T)[i,:]
			muarr, covarr = gp.predict(col,[np.array(p),[0.225,0.2,0.02,1.1]])
			mu=muarr[0]
			cov=np.diag(covarr)[0]
			modcov=cov+obserror*np.identity(1)
			#print mu," ",cov," ",modcov
			cov=modcov
			if (cov > 0.0):
				std=np.sqrt(cov)
			else:
				print("Error, cov = ",cov," p = ",p)
			std=1.0
			cellpostavs[i]+=mu
			thisprob=0.5*(1.0-math.erf((13.-mu)/std))
			prob13[i]+=thisprob
			thisprob=0.5*(1.0-math.erf((37.-mu)/std))
			prob37[i]+=thisprob
			thisprob=0.5*(1.0-math.erf((50.-mu)/std))
			prob50[i]+=thisprob
	bar.next()
bar.finish()
print("For site ",site,": Group average concentration ",grouppostav/numpost," Group average probability of exceedances: ",groupavprob13/numpost," ",groupavprob37/numpost," ",groupavprob50/numpost)
problow=1.0
probmed=1.0
probhigh=1.0
for i in range((y.T).shape[0]):
	problow*=(1.0-prob13[i]/numpost)
	probmed*=(1.0-prob37[i]/numpost)
	probhigh*=(1.0-prob50[i]/numpost)

print("One or more cell probability of exceedances: ",1.0-problow," ",1.0-probmed," ",1.0-probhigh) 
print("Individual cell averages and probabilities of exceedance:")
print("CellID xcoord ycoord cellav ExProb13 ExProb37 ExProb50")
for i in range((y.T).shape[0]):
	print(i," ",celllist[i][0]," ",celllist[i][1]," ",cellpostavs[i]/numpost," ",prob13[i]/numpost," ",prob37[i]/numpost," ",prob50[i]/numpost)
hf.close()
exit()

h5out = h5py.File(outname,'w')
g1=h5out.create_group('parameters')
g2=h5out.create_group('cell list')
g2.create_dataset('cells in group',data=outcelllist)
g2.create_dataset('cells with good NIRAMS data',data=goodcols)

g3=h5out.create_group('emulator')
g1.create_dataset('GWB Group',data=sitename)
g1.create_dataset('GWB Group ID Number',data=site)
g1.create_dataset('Number of cells',data=ncells)

fitOptGP('NIRAMS_GWBGroup%i_YearAv_%i_OptHyp'%(site,year),showplots,1)
g3.create_dataset('Kernel type',data='ExpSquared')
g3.create_dataset('Kernel parameters',data=gp.kernel.vector)
g3.create_dataset('emulation type',data='six year average, centred at %i'%(year))
g3.create_dataset('x',data=x)
g3.create_dataset('y',data=y)
g3.create_dataset('yerr',data=yerr)
g3.create_dataset('group average',data=groupav)
mu, cov = gp.predict(groupav,x)
g3.create_dataset('Inverse covariance for training data',data=np.linalg.inv(cov))

obserror=0.05
nwalkers, ndim = 36, 4
nburnin=2000
nproduction=5000
ranges=np.array([[x1min,x1max],[x2min,x2max],[x3min,x3max],[x4min,x4max]])
if (calibrate):
	g4=h5out.create_group('calibration')
	g4.create_dataset('calibration data',data=calibval)
	g4.create_dataset('observation error',data=obserror)
	g4.create_dataset('sampler',data='emcee')
	g4.create_dataset('nwalkers',data=nwalkers)
	g4.create_dataset('burn in samples',data=nburnin)
	g4.create_dataset('production samples',data=nproduction)
	g4.create_dataset('prior ranges',data=ranges)
	SamplePosterior('NIRAMS_GWBGroup%i_YearAv_%i_OptHyp_Posterior'%(site,year),showplots)

h5out.close()

