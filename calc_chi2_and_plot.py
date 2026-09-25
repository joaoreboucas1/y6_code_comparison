#from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA

#add location of data vector file for plotting
plotfile = "plots/test.png"

# New files:
datavfile2 = "cocoa/LCDM_no_prefactor_no_BMAG.modelvector"
datavfile1 = "cosmolike_lighthouse/LCDM_lighthouse_no_BMAG.modelvector"

d1 = np.genfromtxt(datavfile1)[:,1]
d2 = np.genfromtxt(datavfile2)[:,1]

#use this covariance for redmagic lens sample
covfile = "./cosmosis/DESY6.cov"

ndata = d1.shape[0]
m = np.genfromtxt("./DESY6.mask")[:,1]

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
print(d2/d1*m)
print("chi2 calculated using Y6 extension scale cuts")
inv = LA.inv(cov)
chi = ((d1 - d2)*m).T@inv@((d1-d2)*m)
print("3x2pt: Delta chi2 = %f" %(chi))

NUM_SOURCE_BINS = 4
NUM_LENS_BINS = 6
NUM_SHEAR_BLOCKS = NUM_SOURCE_BINS * (NUM_SOURCE_BINS + 1) // 2
NUM_ANG_BINS = 26
BLOCK_SIZE = NUM_ANG_BINS
SHEAR_SLICE = slice(0, 2*NUM_SHEAR_BLOCKS*NUM_ANG_BINS)
TWOXTWOPT_SLICE = slice(2*NUM_SHEAR_BLOCKS*NUM_ANG_BINS)
GGL_START = 2*NUM_SHEAR_BLOCKS*BLOCK_SIZE
GC_START = 2*NUM_SHEAR_BLOCKS*BLOCK_SIZE + NUM_SOURCE_BINS*NUM_LENS_BINS*BLOCK_SIZE
GGL_SLICE = slice(2*NUM_SHEAR_BLOCKS*NUM_ANG_BINS, 2*NUM_SHEAR_BLOCKS*BLOCK_SIZE + NUM_SOURCE_BINS*NUM_LENS_BINS*BLOCK_SIZE)
GC_SLICE = slice(2*NUM_SHEAR_BLOCKS*BLOCK_SIZE + NUM_SOURCE_BINS*NUM_LENS_BINS*BLOCK_SIZE)
shear1 = d1[:2*NUM_SHEAR_BLOCKS*BLOCK_SIZE]
ggl1 = d1[2*NUM_SHEAR_BLOCKS*BLOCK_SIZE:2*NUM_SHEAR_BLOCKS*BLOCK_SIZE + NUM_SOURCE_BINS*NUM_LENS_BINS*BLOCK_SIZE]
gc1 = d1[2*NUM_SHEAR_BLOCKS*BLOCK_SIZE + NUM_SOURCE_BINS*NUM_LENS_BINS*BLOCK_SIZE:]
twoxtwo_1 = np.concatenate((ggl1, gc1))
shear2 = d2[:2*NUM_SHEAR_BLOCKS*BLOCK_SIZE]
ggl2 = d2[2*NUM_SHEAR_BLOCKS*BLOCK_SIZE:2*NUM_SHEAR_BLOCKS*BLOCK_SIZE + NUM_SOURCE_BINS*NUM_LENS_BINS*BLOCK_SIZE]
gc2 = d2[2*NUM_SHEAR_BLOCKS*BLOCK_SIZE + NUM_SOURCE_BINS*NUM_LENS_BINS*BLOCK_SIZE:]
twoxtwo_2 = np.concatenate((ggl2, gc2))
# 26 bins between 2.5 and ~1000, cut > 250

inv = LA.inv(cov[GGL_START:, GGL_START:])
chi2 = 0
for i in range(0, len(twoxtwo_1)):
	for j in range(0, len(twoxtwo_1)):
		chi2 +=(d1[GGL_START+i]-d2[GGL_START+i])*inv[i,j]*(d1[GGL_START+j]-d2[GGL_START+j])*m[GGL_START+i]*m[GGL_START+j]
print(f"2x2pt: Delta chi2 = {chi2}")

chi = 0.0
inv = LA.inv(cov[GC_START:,GC_START:])
for i in range(0, len(gc1)):
	for j in range(0, len(gc1)):
		chi +=(d1[GC_START+i]-d2[GC_START+i])*inv[i,j]*(d1[GC_START+j]-d2[GC_START+j])*m[GC_START+i]*m[GC_START+j]
print("wtheta: Delta chi2 = %f" %(chi))

chi = 0.0
inv = LA.inv(cov[:GC_START,:GC_START])
for i in range(0, len(shear1)):
	for j in range(0, len(shear1)):
		chi +=(d1[i]-d2[i])*inv[i,j]*(d1[j]-d2[j])*m[i]*m[j]
print("shear: Delta chi2 = %f" %(chi))

nxip = 260
nxim = 260
nw = 156
nggl = ndata - nw - nxip -nxim
s = np.sqrt(np.diag(cov))

plt.figure(figsize=(8,8), dpi=400)
fs = 18
plt.subplot(4,2,1)
plt.yscale('log')
plt.ylim(2.e-7,1.2e-4)
plt.xlim(0,nxip-1)
#plt.title(r'$\xi_+$')
plt.ylabel(r'$\xi_+$', fontsize = fs)
plt.errorbar(ind,d1,s,marker='o', color='k',linestyle = '',markersize = 0.5,alpha = 0.25)
plt.plot(ind,d1,marker='o', color='r',linestyle = '',markersize = 1.5)


plt.subplot(4,2,2)
plt.ylim(-0.1,0.1)
plt.plot([0,1000],[0,0],linestyle ='--',color='k')
plt.xlim(0,nxip-1)
plt.ylabel(r'(d2-d1)/d2', fontsize = fs)
plt.errorbar(ind,d1*0,s/np.abs(d1),marker='o', color='k',linestyle = '',markersize = 0.0,alpha = 0.1)
plt.plot(ind[ind0],(d2[ind0]-d1[ind0])/d2[ind0],marker='x', color='k',linestyle = '',markersize = 1.0)
plt.plot(ind[ind1],(d2[ind1]-d1[ind1])/d2[ind1],marker='o', color='r',linestyle = '',markersize = 1.0)

plt.subplot(4,2,3)
plt.yscale('log')
plt.ylim(2.e-7,6.e-5)
plt.xlim(nxip,nxip+nxim-1)
#plt.title(r'$\xi_-$')
plt.ylabel(r'$\xi_-$', fontsize = fs)
plt.errorbar(ind,d1,s,marker='o', color='k',linestyle = '',markersize = 0.5,alpha = 0.25)
plt.plot(ind,d1,marker='o', color='r',linestyle = '',markersize = 1.5)

plt.subplot(4,2,4)
plt.ylim(-0.1,0.1)
plt.plot([0,1000],[0,0],linestyle ='--',color='k')
plt.xlim(nxip,nxip+nxim-1)
plt.ylabel(r'(d2-d1)/d2', fontsize = fs)
plt.errorbar(ind,d1*0,s/np.abs(d1),marker='o', color='k',linestyle = '',markersize = 0.0,alpha = 0.1)
plt.plot(ind[ind0],(d2[ind0]-d1[ind0])/d2[ind0],marker='x', color='k',linestyle = '',markersize = 1.0)
plt.plot(ind[ind1],(d2[ind1]-d1[ind1])/d2[ind1],marker='o', color='r',linestyle = '',markersize = 1.0)

plt.subplot(4,2,5)

plt.yscale('log')
plt.ylim(2.e-6,2.5e-3)
plt.xlim(nxip+nxim,nxip+nxim+nggl-1)
#plt.title(r'$\gamma_t$')
plt.ylabel(r'$\gamma_t$', fontsize = fs)
plt.errorbar(ind,d1,s,marker='o', color='k',linestyle = '',markersize = 0.5,alpha = 0.2)
plt.plot(ind,d1,marker='o', color='r',linestyle = '',markersize = 1.5)
#plt.plot(ind,d3,linestyle = '-')

plt.subplot(4,2,6)
plt.ylim(-0.2,0.2)
plt.plot([0,1000],[0,0],linestyle ='--',color='k')
plt.xlim(nxip+nxim,nxip+nxim+nggl-1)
#plt.title(r'$\gamma_t$')
plt.ylabel(r'(d2-d1)/d2', fontsize = fs)
plt.errorbar(ind,d1*0,s/np.abs(d1),marker='o', color='k',linestyle = '',markersize = 0.0,alpha = 0.1)
plt.plot(ind[ind0],(d2[ind0]-d1[ind0])/d2[ind0],marker='x', color='k',linestyle = '',markersize = 1.0)
plt.plot(ind[ind1],(d2[ind1]-d1[ind1])/d2[ind1],marker='o', color='r',linestyle = '',markersize = 1.0)


plt.subplot(4,2,7)
plt.yscale('log')
plt.ylim(1.e-4,0.6)
plt.xlim(nxip+nxim+nggl,ndata)
#plt.title(r'$w$')
plt.ylabel(r'$w$', fontsize = fs)
plt.xlabel(r'bin number', fontsize = fs)
plt.errorbar(ind,d1,s,marker='o', color='k',linestyle = '',markersize = 0.5,alpha = 0.4)
plt.plot(ind,d1,marker='o', color='r',linestyle = '',markersize = 1.5)
#plt.plot(ind,d3,linestyle = '-')


plt.subplot(4,2,8)
plt.ylim(-0.1,0.1)
plt.plot([0,1000],[0,0],linestyle ='--',color='k')
plt.xlim(nxip+nxim+nggl,ndata)
#plt.title(r'$w$')
plt.xlabel(r'bin number', fontsize = 18)
plt.ylabel(r'(d2-d1)/d2', fontsize = fs)
plt.errorbar(ind,d1*0,s/np.abs(d1),marker='o', color='k',linestyle = '',markersize = 0.0,alpha = 0.1)
plt.plot(ind[ind0],(d2[ind0]-d1[ind0])/d2[ind0],marker='x', color='k',linestyle = '',markersize = 1.0)
plt.plot(ind[ind1],(d2[ind1]-d1[ind1])/d2[ind1],marker='o', color='r',linestyle = '',markersize = 1.0)

plt.tight_layout()
plt.savefig(plotfile,dpi=400)

