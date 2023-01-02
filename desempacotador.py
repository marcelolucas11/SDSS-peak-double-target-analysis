import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

filelist =  np.loadtxt('lsst_thaisa_supersample.csv', usecols=6, dtype='str', unpack=True, delimiter='; ')
Z = np.loadtxt('lsst_thaisa_supersample.csv', usecols=7, dtype='float', unpack=True, delimiter='; ')
plate, mjd, fiberID = np.loadtxt('lsst_thaisa_supersample.csv', usecols=(3,4,5), dtype='int', unpack=True, delimiter='; ')

for i in range(216,Z.shape[0]):
	
	if Z[i] < 0.4:
		hdul = fits.open('..DR16Q/other_spectra/' + filelist[i])
				
	else:
		hdul = fits.open('..DR16Q/z_subsample/' + filelist[i])
		
	flux = hdul[1].data['FLUX'] * 1e-17
	wl = 10**hdul[1].data['LOGLAM']
	
	rest_wl = 10**hdul[1].data['LOGLAM'] / (1. + Z[i])
	
	fig, ax = plt.subplots(1,1, figsize=(12,4))

	fig.suptitle(filelist[i])
	ax.plot(rest_wl, flux, color='firebrick')
	ax.axvline(x=5007., linestyle='--', color='black', label='[O III]')
	ax.axvline(x=4959., linestyle='--', color='grey', label='[O III]')
	ax.axvline(x=4861., linestyle='--', color='royalblue', label='H$\\beta$')
	ax.axvline(x=4340., linestyle='--', color='purple', label='H$\gamma$')
	ax.axvline(x=6563., linestyle='--', color='red', label='H$\\alpha$')
	ax.axvline(x=6553., linestyle='--', color='green',label'N II'=)
	ax.axvline(x=6607., linestyle='--', color='green',label'N II'=)
	ax.set_xlabel('Rest Wavelength [$\AA$]')
	ax.set_ylabel('Flux Density [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]')

	plt.legend()
	plt.tight_layout()
	plt.show(block=True)
	print('File Number: ', i)
	

