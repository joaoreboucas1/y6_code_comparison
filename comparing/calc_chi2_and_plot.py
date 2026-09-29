#from astropy.io import fits
import os
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA

#directory for plots and logs
plotdir = "plots"
os.makedirs(plotdir, exist_ok=True)

# New files:
datavfiles = {
	"cosmolike": "../cosmolike_lighthouse/LCDM_lighthouse.modelvector",
	"cocoa": "../cocoa/LCDM.modelvector",
	"cosmosis": "COSMOSIS.modelvector",
}
dvs = {name: np.genfromtxt(f)[:,1] for name, f in datavfiles.items()}

# pairs to compare: (d1, d2)
pairs = [("cosmolike", "cocoa"), ("cosmolike", "cosmosis"), ("cocoa", "cosmosis")]

#use this covariance for redmagic lens sample
covfile = "DESY6.cov"

ndata = dvs["cosmolike"].shape[0]
m = np.genfromtxt("DESY6.mask")[:,1]

ind1 = np.where(m)
ind0 = np.where(m-1.0)
data = np.genfromtxt(covfile)
cov = np.zeros((ndata,ndata))
for i in range(0,data.shape[0]):
	if int(data[i,1]) >= ndata or int(data[i,0]) >= ndata: continue
	cov[int(data[i,0]),int(data[i,1])] = data[i,2]
	cov[int(data[i,1]),int(data[i,0])] = data[i,2]
	if (int(data[i,0])-int(data[i,1])):
		cov[int(data[i,0]),int(data[i,1])]*= m[int(data[i,0])]*m[int(data[i,1])]
		cov[int(data[i,1]),int(data[i,0])]*= m[int(data[i,0])]*m[int(data[i,1])]
#s now contains sqrt(cov[i,i])
s = np.sqrt(np.diag(cov))
ind = np.arange(0,ndata)

NUM_SOURCE_BINS = 4
NUM_LENS_BINS = 6
NUM_SHEAR_BLOCKS = NUM_SOURCE_BINS * (NUM_SOURCE_BINS + 1) // 2
NUM_ANG_BINS = 26
BLOCK_SIZE = NUM_ANG_BINS
GGL_START = 2*NUM_SHEAR_BLOCKS*BLOCK_SIZE
GC_START = 2*NUM_SHEAR_BLOCKS*BLOCK_SIZE + NUM_SOURCE_BINS*NUM_LENS_BINS*BLOCK_SIZE

nxip = 260
nxim = 260
nw = 156
nggl = ndata - nw - nxip -nxim

def compare(name1, d1, name2, d2):
	logfile = open(f"{plotdir}/{name1}_vs_{name2}.log", "w")
	def log(msg):
		print(msg)
		print(msg, file=logfile)

	log(f"===== {name1} vs {name2} =====")
	log("chi2 calculated using Y6 extension scale cuts")
	inv = LA.inv(cov)
	chi = ((d1 - d2)*m).T@inv@((d1-d2)*m)
	log("3x2pt: Delta chi2 = %f" %(chi))

	shear1 = d1[:2*NUM_SHEAR_BLOCKS*BLOCK_SIZE]
	ggl1 = d1[2*NUM_SHEAR_BLOCKS*BLOCK_SIZE:2*NUM_SHEAR_BLOCKS*BLOCK_SIZE + NUM_SOURCE_BINS*NUM_LENS_BINS*BLOCK_SIZE]
	gc1 = d1[2*NUM_SHEAR_BLOCKS*BLOCK_SIZE + NUM_SOURCE_BINS*NUM_LENS_BINS*BLOCK_SIZE:]
	twoxtwo_1 = np.concatenate((ggl1, gc1))
	# 26 bins between 2.5 and ~1000, cut > 250

	inv = LA.inv(cov[GGL_START:, GGL_START:])
	chi2 = 0
	for i in range(0, len(twoxtwo_1)):
		for j in range(0, len(twoxtwo_1)):
			chi2 +=(d1[GGL_START+i]-d2[GGL_START+i])*inv[i,j]*(d1[GGL_START+j]-d2[GGL_START+j])*m[GGL_START+i]*m[GGL_START+j]
	log(f"2x2pt: Delta chi2 = {chi2}")

	chi = 0.0
	inv = LA.inv(cov[GGL_START:GC_START,GGL_START:GC_START])
	for i in range(0, len(ggl1)):
		for j in range(0, len(ggl1)):
			chi +=(d1[GGL_START+i]-d2[GGL_START+i])*inv[i,j]*(d1[GGL_START+j]-d2[GGL_START+j])*m[GGL_START+i]*m[GGL_START+j]
	log("ggl: Delta chi2 = %f" %(chi))

	chi = 0.0
	inv = LA.inv(cov[GC_START:,GC_START:])
	for i in range(0, len(gc1)):
		for j in range(0, len(gc1)):
			chi +=(d1[GC_START+i]-d2[GC_START+i])*inv[i,j]*(d1[GC_START+j]-d2[GC_START+j])*m[GC_START+i]*m[GC_START+j]
	log("wtheta: Delta chi2 = %f" %(chi))

	chi = 0.0
	inv = LA.inv(cov[:GGL_START,:GGL_START])
	for i in range(0, len(shear1)):
		for j in range(0, len(shear1)):
			chi +=(d1[i]-d2[i])*inv[i,j]*(d1[j]-d2[j])*m[i]*m[j]
	log("shear: Delta chi2 = %f" %(chi))
	logfile.close()

	plt.figure(figsize=(8,8), dpi=400)
	fs = 18
	plt.subplot(4,2,1)
	plt.yscale('log')
	plt.ylim(2.e-7,1.2e-4)
	plt.xlim(0,nxip-1)
	#plt.title(r'$\xi_+$')
	plt.ylabel(r'$\xi_+$', fontsize = fs)
	plt.errorbar(ind,d1,s,marker='o', color='gray',linestyle = '',markersize = 0.5,alpha = 0.3)
	plt.plot(ind,d1,marker='o', color='r',linestyle = '',markersize = 1.5, label='left: data vector 1')


	plt.subplot(4,2,2)
	plt.ylim(-0.1,0.1)
	plt.plot([0,1000],[0,0],linestyle ='--',color='k')
	plt.xlim(0,nxip-1)
	plt.ylabel(r'$\frac{\xi_{+,2}-\xi_{+,1}}{\xi_{+,2}}$', fontsize = fs)
	plt.fill_between(ind,-s/np.abs(d1),s/np.abs(d1),step='mid', color='gray',alpha = 0.3,linewidth = 0, label=r'$1\sigma$ error (Y6 covariance)')
	plt.plot(ind[ind0],(d2[ind0]-d1[ind0])/d2[ind0],marker='x', color='k',linestyle = '',markersize = 1.0, label='right: removed by scale cuts')
	plt.plot(ind[ind1],(d2[ind1]-d1[ind1])/d2[ind1],marker='o', color='r',linestyle = '',markersize = 1.0, label='right: kept by scale cuts')

	plt.subplot(4,2,3)
	plt.yscale('log')
	plt.ylim(2.e-7,6.e-5)
	plt.xlim(nxip,nxip+nxim-1)
	#plt.title(r'$\xi_-$')
	plt.ylabel(r'$\xi_-$', fontsize = fs)
	plt.errorbar(ind,d1,s,marker='o', color='gray',linestyle = '',markersize = 0.5,alpha = 0.3)
	plt.plot(ind,d1,marker='o', color='r',linestyle = '',markersize = 1.5)

	plt.subplot(4,2,4)
	plt.ylim(-0.1,0.1)
	plt.plot([0,1000],[0,0],linestyle ='--',color='k')
	plt.xlim(nxip,nxip+nxim-1)
	plt.ylabel(r'$\frac{\xi_{-,2}-\xi_{-,1}}{\xi_{-,2}}$', fontsize = fs)
	plt.fill_between(ind,-s/np.abs(d1),s/np.abs(d1),step='mid', color='gray',alpha = 0.3,linewidth = 0)
	plt.plot(ind[ind0],(d2[ind0]-d1[ind0])/d2[ind0],marker='x', color='k',linestyle = '',markersize = 1.0)
	plt.plot(ind[ind1],(d2[ind1]-d1[ind1])/d2[ind1],marker='o', color='r',linestyle = '',markersize = 1.0)

	plt.subplot(4,2,5)

	plt.yscale('log')
	plt.ylim(2.e-6,2.5e-3)
	plt.xlim(nxip+nxim,nxip+nxim+nggl-1)
	#plt.title(r'$\gamma_t$')
	plt.ylabel(r'$\gamma_t$', fontsize = fs)
	plt.errorbar(ind,d1,s,marker='o', color='gray',linestyle = '',markersize = 0.5,alpha = 0.3)
	plt.plot(ind,d1,marker='o', color='r',linestyle = '',markersize = 1.5)
	#plt.plot(ind,d3,linestyle = '-')

	plt.subplot(4,2,6)
	plt.ylim(-0.2,0.2)
	plt.plot([0,1000],[0,0],linestyle ='--',color='k')
	plt.xlim(nxip+nxim,nxip+nxim+nggl-1)
	#plt.title(r'$\gamma_t$')
	plt.ylabel(r'$\frac{\gamma_{t,2}-\gamma_{t,1}}{\gamma_{t,2}}$', fontsize = fs)
	plt.fill_between(ind,-s/np.abs(d1),s/np.abs(d1),step='mid', color='gray',alpha = 0.3,linewidth = 0)
	plt.plot(ind[ind0],(d2[ind0]-d1[ind0])/d2[ind0],marker='x', color='k',linestyle = '',markersize = 1.0)
	plt.plot(ind[ind1],(d2[ind1]-d1[ind1])/d2[ind1],marker='o', color='r',linestyle = '',markersize = 1.0)


	plt.subplot(4,2,7)
	plt.yscale('log')
	plt.ylim(1.e-4,0.6)
	plt.xlim(nxip+nxim+nggl,ndata)
	#plt.title(r'$w$')
	plt.ylabel(r'$w$', fontsize = fs)
	plt.xlabel(r'bin number', fontsize = fs)
	plt.errorbar(ind,d1,s,marker='o', color='gray',linestyle = '',markersize = 0.5,alpha = 0.3)
	plt.plot(ind,d1,marker='o', color='r',linestyle = '',markersize = 1.5)
	#plt.plot(ind,d3,linestyle = '-')


	plt.subplot(4,2,8)
	plt.ylim(-0.1,0.1)
	plt.plot([0,1000],[0,0],linestyle ='--',color='k')
	plt.xlim(nxip+nxim+nggl,ndata)
	#plt.title(r'$w$')
	plt.xlabel(r'bin number', fontsize = 18)
	plt.ylabel(r'$\frac{w_2-w_1}{w_2}$', fontsize = fs)
	plt.fill_between(ind,-s/np.abs(d1),s/np.abs(d1),step='mid', color='gray',alpha = 0.3,linewidth = 0)
	plt.plot(ind[ind0],(d2[ind0]-d1[ind0])/d2[ind0],marker='x', color='k',linestyle = '',markersize = 1.0)
	plt.plot(ind[ind1],(d2[ind1]-d1[ind1])/d2[ind1],marker='o', color='r',linestyle = '',markersize = 1.0)

	plt.suptitle(f"1: {name1.capitalize()}, 2: {name2.capitalize()}", fontsize = fs)
	plt.figlegend(loc='upper center', bbox_to_anchor=(0.5, 0.95), ncol=2, fontsize=10, markerscale=4)
	plt.tight_layout(rect=(0, 0, 1, 0.9))
	plt.savefig(f"{plotdir}/{name1}_vs_{name2}.pdf",dpi=400)
	plt.close()

for name1, name2 in pairs:
	compare(name1, dvs[name1], name2, dvs[name2])
