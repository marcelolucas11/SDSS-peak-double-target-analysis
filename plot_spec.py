## code for plot spec of AGN 
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import lineid_plot

plt.close('all')

hdul = fits.open('spec-1624-53386-0032.fits')
flux = hdul[1].data['FLUX'] * 1e-17
wl = 10**hdul[1].data['LOGLAM']

z = (1. + 0.007548481)
 
rest_wl = 10**hdul[1].data['LOGLAM'] 

fig, ax = plt.subplots(1,1, figsize=(12,4))

fig.suptitle('spec-1624-53386-0032.fits')
ax.plot(rest_wl, flux, color='royalblue')

#linhas
line_wave = np.array([4340,4861,4959,5007,6553,6563,6607])*z
#line_wave = lines
line_label1 = ['H$\gamma$','H$\\beta$','[O$\,$III]','[O$\,$III]','[N$\,$II]','H$\\alpha$','[N$\,$II]']

#lineid_plot.plot_line_ids(rest_wl, flux, line_wave, line_label1,ax=ax)
ax.set_xlabel('Rest Wavelength [$\AA$]')
ax.set_ylabel('Flux Density [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]')

# Lines
fluxmax = np.nanmax(flux)
for i in range(len(line_wave)):
    wli = line_wave[i]
    fluxmaxi = flux[np.argmin(abs(rest_wl - wli))]

    ax.vlines(wli, ymin=fluxmaxi, ymax=fluxmax, 
              color='black', ls=':', lw=1, alpha=.9)

# Text
#line_wave2 = np.array([4340,4861,4959,5007,6563])*z
#line_label2 = ['H$\gamma$','H$\\beta$','[O$\,$III]$\lambda$4959','[O$\,$III]$\lambda$5007','[N$\,$II]+H$\\alpha$']
line_wave2 = np.array([4340,4861,4983,6563])*z
line_label2 = ['H$\gamma$','H$\\beta$','[O$\,$III]','[N$\,$II]+H$\\alpha$']

for i in range(len(line_wave2)):
    wli2 = line_wave2[i]
    fluxmaxi2 = fluxmax * 1.05

    labeli2 = line_label2[i]
    #ax.text(wli2, fluxmaxi2,labeli2, 
    #        color='black', fontsize='x-small', 
    #        horizontalalignment='center', verticalalignment='bottom',
    #        rotation='vertical', rotation_mode='default')
    ax.text(wli2, fluxmaxi2,labeli2, 
            color='black', fontsize='x-small', 
            horizontalalignment='center', verticalalignment='bottom',
            rotation='horizontal')



#plt.legend()
plt.tight_layout()
plt.show()
# Save figure
#plt.savefig('teste.pdf', bbox_inches='tight')

print(line_wave)


