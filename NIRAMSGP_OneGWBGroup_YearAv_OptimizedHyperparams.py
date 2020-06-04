import sys
sys.path.insert(1,'/Library/Python/2.7/site-packages')
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import george
import scipy.optimize as op
import emcee
import corner as triangle
import pygtc
import h5py

from george.kernels import ExpSquaredKernel

def nll(p):
    gp.kernel[:]=p
    ll=0.0
    ct=0
    for col in y.T:
      if (goodcols[ct]):
      	gp.compute(x, yerr[:,ct])
      	ll+=gp.lnlikelihood(col,quiet=True)
      	ct+=1
    #col=(y.T)[38,:]
    #ll+=gp.lnlikelihood(col,quiet=True)
    return -ll if np.isfinite(ll) else 1e25

def grad_nll(p):
    gp.kernel[:]=p
    grad_ll=0.0
    ct=0
    for col in y.T:
      if (goodcols[ct]):
      	gp.compute(x, yerr[:,ct])
      	grad_ll-=gp.grad_lnlikelihood(col, quiet=True)
      	ct+=1
    return grad_ll

def lnprior(p):
    #print "p = ",p
    #if (p.shape[0]==1):
    #	q=p[0]
    #else:
    #	q=p
    #print "q = ",q
    #print q.shape
    #if np.any((p[0] < 0.15)+(p[0] > 0.35)):
    q=p
    if ((q[0] < x1min) or (q[0] > x1max)):
        #print "Prob 0"
        return -np.inf
    #if np.any((p[1] < 0.1)+(p[1] > 0.3)):
    if ((q[1] < x2min) or (q[1] > x2max)):
        #print "Prob 1"
        return -np.inf
    #if np.any((p[2] > 0.04)+(p [2] < 0)):
    if ((q[2] > x3max) or (q [2] < x3min)):
        #print "Prob 2"
        return -np.inf
    #if np.any((p[3] < 0.9)+(p[3] > 1.3)):
    if ((q[3] < x4min) or (q[3] > x4max)):
        #print "Prob 3"
        return -np.inf
    lnprior=0.
    return lnprior

def lnlike(p):
    prival=lnprior(p)
    if (prival == -np.inf):
        return -np.inf
    points=[]
    points.append(p)
    points.append([0.225,0.2,0.02,1.1])
    muarr, covarr = gp.predict(groupav,np.array(points))
    mu=muarr[0]
    cov=np.diag(covarr)[0]
    #std = np.sqrt(np.diag(cov))
    modcov=cov+obserror*np.identity(1)
    if (modcov < 0.0):
	return -np.inf
    loglike=-0.5*((calibval-mu)*(calibval-mu)/modcov+np.log(modcov))
    #loglike=-0.5*(np.dot(obsdata[:,1]-mu,np.dot(np.linalg.inv(modcov),obsdata[:,1]-mu))+np.log(np.linalg.det(modcov)))
    #loglike=-0.5*(np.linalg.norm(obsdata[:,1]-mu)*np.linalg.norm(obsdata[:,1]-mu)/(std*std+obserror*obserror)+np.log(np.linalg.det(modcov)))
    #print p,mu,cov,loglike
    return loglike

def obslnprob(p):
    prival=lnprior(p)
    if (prival == -np.inf):
        return -np.inf
    val=prival+lnlike(p)
    return val

def makefig(x1v,mu,std,thisdatax,thisdatay,thisdatayerr,title,filename,showplots):
	plt.plot(x1v,mu,'k-',linewidth=2)
	plt.plot(x1v,mu+std,'k-',linewidth=1)
	plt.plot(x1v,mu-std,'k-',linewidth=1)
	plt.fill_between(x1v,mu-std,mu+std,facecolor='grey',alpha=0.5)
	plt.errorbar(thisdatax,thisdatay,thisdatayerr,fmt='o')
	plt.xlabel('x')
	plt.ylabel('y')
	plt.title(title)
	plt.savefig(filename)
	if (showplots):
		plt.show()
	plt.close()

def makeplots(col,whicherrcol,basename,showplots):
	ngrid=50
	x1v=np.linspace(0.15, 0.35, ngrid)
	x2v=fixx2v*np.ones_like(x1v)
	x3v=fixx3v*np.ones_like(x1v)
	x4v=fixx4v*np.ones_like(x1v)
	xpred = np.array(zip(x1v,x2v,x3v,x4v))
	mu, cov = gp.predict(col,xpred)
	std = np.sqrt(np.diag(cov))
	thisdatax=x[(x[:,1]==fixx2v) & (x[:,2]==fixx3v) & (x[:,3]==fixx4v)][:,0]
	thisdatay=col[(x[:,1]==fixx2v) & (x[:,2]==fixx3v) & (x[:,3]==fixx4v)]
	thisdatayerr=yerr[(x[:,1]==fixx2v) & (x[:,2]==fixx3v) & (x[:,3]==fixx4v),whicherrcol]
	filename='{base}_x1.png'.format(base=basename)
	makefig(x1v,mu,std,thisdatax,thisdatay,thisdatayerr,'Organic_N',filename,showplots)

	ngrid=50
	x2v=np.linspace(0.1, 0.3, ngrid)
	x1v=fixx1v*np.ones_like(x2v)
	x3v=fixx3v*np.ones_like(x2v)
	x4v=fixx4v*np.ones_like(x2v)
	xpred = np.array(zip(x1v,x2v,x3v,x4v))
	mu, cov = gp.predict(col,xpred)
	std = np.sqrt(np.diag(cov))

	thisdatax=x[(x[:,0]==fixx1v) & (x[:,2]==fixx3v) & (x[:,3]==fixx4v)][:,1]
	thisdatay=col[(x[:,0]==fixx1v) & (x[:,2]==fixx3v) & (x[:,3]==fixx4v)]
	thisdatayerr=yerr[(x[:,0]==fixx1v) & (x[:,2]==fixx3v) & (x[:,3]==fixx4v),whicherrcol]
	
	filename='{base}_x2.png'.format(base=basename)
	makefig(x2v,mu,std,thisdatax,thisdatay,thisdatayerr,'Mineralisation',filename,showplots)

	ngrid=50
	x3v=np.linspace(0.005, 0.035, ngrid)
	x2v=fixx2v*np.ones_like(x3v)
	x1v=fixx1v*np.ones_like(x3v)
	x4v=fixx4v*np.ones_like(x3v)
	xpred = np.array(zip(x1v,x2v,x3v,x4v))
	mu, cov = gp.predict(col,xpred)
	std = np.sqrt(np.diag(cov))

	thisdatax=x[(x[:,1]==fixx2v) & (x[:,0]==fixx1v) & (x[:,3]==fixx4v)][:,2]
	thisdatay=col[(x[:,1]==fixx2v) & (x[:,0]==fixx1v) & (x[:,3]==fixx4v)]
	thisdatayerr=yerr[(x[:,1]==fixx2v) & (x[:,0]==fixx1v) & (x[:,3]==fixx4v),whicherrcol]

	filename='{base}_x3.png'.format(base=basename)
	makefig(x3v,mu,std,thisdatax,thisdatay,thisdatayerr,'Denitrification',filename,showplots)

	ngrid=50
	x4v=np.linspace(0.95, 1.25, ngrid)
	x2v=fixx2v*np.ones_like(x4v)
	x3v=fixx3v*np.ones_like(x4v)
	x1v=fixx1v*np.ones_like(x4v)
	xpred = np.array(zip(x1v,x2v,x3v,x4v))
	mu, cov = gp.predict(col,xpred)
	std = np.sqrt(np.diag(cov))

	thisdatax=x[(x[:,1]==fixx2v) & (x[:,2]==fixx3v) & (x[:,0]==fixx1v)][:,3]
	thisdatay=col[(x[:,1]==fixx2v) & (x[:,2]==fixx3v) & (x[:,0]==fixx1v)]
	thisdatayerr=yerr[(x[:,1]==fixx2v) & (x[:,2]==fixx3v) & (x[:,0]==fixx1v),whicherrcol]

	filename='{base}_x4.png'.format(base=basename)
	makefig(x4v,mu,std,thisdatax,thisdatay,thisdatayerr,'Leaching',filename,showplots)

def fitOptGP(filebasename,showplots,nplots):
	ct=0
	while (~(goodcols[ct])):
		ct+=1
	gp.compute(x, yerr[:,ct])
	#print(gp.kernel.vector)
	#print(gp.lnlikelihood(y,quiet=True))
	#print(gp.grad_lnlikelihood(y,quiet=True))
	#gp.kernel[:]=np.log([66.0,0.059,0.0014,0.000063,0.35])
	p0=gp.kernel.vector
	#print(p0)
	#print(np.exp(gp.kernel.vector))

	#results=op.minimize(nll,p0,method='Nelder-Mead')
	results=op.minimize(nll,p0,method='Powell')
	#results=op.minimize(nll,p0,jac=grad_nll)
	gp.kernel[:]=results.x
	print(np.exp(gp.kernel.vector))

	for i in range(nplots):
		idx=np.random.randint((y.T).shape[0],size=1)
		while (~(goodcols[idx[0]])):
			idx=np.random.randint((y.T).shape[0],size=1)
		selcols=(y.T)[idx,:]
		filename='{base}_cell{cell}'.format(base=filebasename,cell=idx[0])
		for col in selcols:
			makeplots(col,idx,filename,showplots)
	filename='{base}_GWBGroupAverage'.format(base=filebasename)
	#print "Doing group"
	makeplots(groupav,0,filename,showplots)

def SamplePosterior(filebasename,showplots):
    #print(np.exp(gp.kernel.vector))
    #obsll= lambda *args: -lnlike(*args)
    #results=op.minimize(obsll,[fixx1v,fixx2v,fixx3v,fixx4v])
    #print(results.x)
    #exit()
    sampler = emcee.EnsembleSampler(nwalkers, ndim, obslnprob)
    #p0 = [results.x + 1e-4 * np.random.randn(ndim)
    #      for i in range(nwalkers)]
    p0 = [[fixx1v,fixx2v,fixx3v,fixx4v] + 1e-4 * np.random.randn(ndim)
          for i in range(nwalkers)]
    #print p0,obslnprob(p0[0],y)
    print("Running burn-in")
    p0, _, _ = sampler.run_mcmc(p0, nburnin)
    print("Running production chain")
    sampler.run_mcmc(p0, nproduction)
    samples=sampler.flatchain
    g4.create_dataset('posterior samples',data=samples)
    labels=[r"$x_1$",r"$x_2$",r"$x_3$",r"$x_4$"]
    #fig = triangle.corner(samples[:, 2:], labels=labels)
    #print np.array(samples)
    #fig = triangle.corner(np.array(samples), bins=100, labels=labels)
    nburn=nburnin*nwalkers
    fig=pygtc.plotGTC(chains=[samples[nburn:,:]])
    fig.savefig("%s_Triangle.png"%(filebasename),dpi=150)
    if showplots:
        fig.show()

if ((len(sys.argv) != 3)):
    print "Arguments: 1. site ID"
    print "           2. show plots on screen (1) or not (0)?"
    exit()
site=int(sys.argv[1])
showplots=int(sys.argv[2])

datafile=r'../../data/IndividualGWBGroupData_2017/GWBGroupData_YearAvsProc_%i' % (site)
calibfile=r'../data/NewData2017/GWBody6YrMeasuredAv.txt'
outname=r'NIRAMS_Emulator_GWBGroup%i.h5' % (site)
f=open(datafile,"r")
rawdata=np.loadtxt(f)
f.close()

f=open(calibfile,"r")
#calibdata=np.loadtxt(f,dtype={'names': ['site label', 'site name', 'Average concentration'], 'formats': ['i2', 'S', 'f6']})
calibdata=[line.rstrip('\n') for line in f]
thiscalibdata=((calibdata[site-1]).rstrip('\r')).split(' ')
sitename=thiscalibdata[1]
calibrate=True
if (thiscalibdata[2]=='NA'):
	print "Calibration data not available, only emulating"
	calibrate=False
if (calibrate):
	calibval=np.float(thiscalibdata[2])
f.close()
#print calibval

cellfile=r'../../data/IndividualGWBGroupData_2017/GWBGroupCellList_%i' % (site)
f=open(cellfile,"r")
celllist=[line.rstrip('\n') for line in f]
celllist.pop(0)
thiscelllist=[]
for x in celllist:
	thiscelllist.append(x.split(' '))
outcelllist=(np.array(thiscelllist)).astype(int)

print 'Processing site number %i, name %s' % (site,sitename)

#transdata=np.transpose(np.array(rawdata))
year=2012

data=rawdata[rawdata[:,4]==year]
x=data[:,0:4]
y=data[:,5:]
yerr=0.001 * np.ones_like(y)

goodcols=~np.all(np.isnan(y),axis=0)

yerr[np.isnan(y)]=np.inf
y[np.isnan(y)]=0.0
#print(x)
#print(y)
groupav=np.array(np.sum(y.T,axis=0)/((y.T).shape[0]))
groupav.reshape((1,(y.T).shape[1]))
#print groupav

ncells=y.shape[1]

kernel=100.0*ExpSquaredKernel([1.08,0.0015,0.00005,0.15],4)
#kernel=ExpSquaredKernel([0.03,0.01,0.0001,0.1],4)
#kernel=ExpSquaredKernel([10,10,10,10],4)
#kernel=ExpSquaredKernel(0.0001,4)
gp=george.GP(kernel)

fixx1v=0.25
fixx2v=0.2
fixx3v=0.02
fixx4v=1.1

#fixx1v=0.225
#fixx2v=0.22
#fixx3v=0.017
#fixx4v=1.06

x1min=0.1
x2min=0.1
x3min=0.005
x4min=0.95
#x1min=0.
#x2min=0.
#x3min=0.
#x4min=0.

x1max=0.35
x2max=0.3
x3max=0.035
x4max=1.25
#x1max=1.
#x2max=1.
#x3max=0.15
#x4max=10.

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
#print mu
#print np.sqrt(np.diag(cov))
#xpred=[[0.225,0.2,0.02,1.1],[0.275,0.2,0.02,1.1],[0.29,0.2,0.02,1.1]]
#mu, cov = gp.predict(groupav,xpred)
#print mu
#print np.sqrt(np.diag(cov))
#exit()

g3.create_dataset('Inverse covariance for training data',data=np.linalg.inv(cov))

obserror=0.01
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

