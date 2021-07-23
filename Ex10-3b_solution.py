import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import cm, rcParams,rc
# number 3b)
# My local settings to make fonts pretty, optional
LocalMac = True
if LocalMac:
    import os
    os.environ["PATH"] += os.pathsep + '/usr/local/texlive/2019/bin/x86_64-darwin'

fontsize=7
font = {'family' : 'serif',
        'serif': ['Computer Modern'],
        'size'   : fontsize}
rcParams['mathtext.fontset']= 'cm'
rc('text', usetex=True)
rc('font', **font)


def b_n(n):
	return 2*n - 1./3. + 4./(405.*n) + 46./(25515.*n**2)

def sersic_profile(r, n, mu_e, R_e):
	I_e = 10**(-mu_e/2.5)
	return I_e * np.exp( -b_n(n)*((r/R_e)**(1./n) - 1.))

# Load data from txt-file
# radius[arcsec], counts, mag[mag/"^2]
rad,counts,mag = np.loadtxt('galaxy.txt', delimiter='    ', usecols=(1, 2,3 ), unpack=True)
# Fix parameters

mu_e = 19.9012 # mag/arcsec^2
R_e  = 30.1617   # arcsec

# Define arrays
R_array = np.logspace(-2.,3.,500)
n_array = np.array([1.,2.,3.,5.,6.])

# Calculate sersic profiles for different n
sersic_profiles_n = np.array([sersic_profile(rad, n, mu_e, R_e) for n in n_array])
sersic_profile_4  = sersic_profile(rad, 4, mu_e, R_e)
mag_sersic4 = -2.5*np.log10(sersic_profile_4)
residual_n4 = mag -mag_sersic4
std_n4 = np.sqrt(np.mean(abs(residual_n4 - residual_n4.mean())**2))

mag_sersic_n = -2.5*np.log10( sersic_profiles_n )

residuals = mag[np.newaxis,:] - mag_sersic_n
std = np.sqrt(np.mean(abs(residuals - residuals.mean(axis=1)[:,np.newaxis])**2, axis=1))
# Make the plot
figsizex = 5
f = plt.figure( figsize=(figsizex,figsizex*0.6), dpi=200 )

# define colors for lines
jet    = cm.get_cmap('jet')
colors = jet(np.linspace(0.7, 1, len(n_array)))

plt.axes( [0.1,0.15,0.85,0.75] ) # [x,y,size_x, size_y] 
for i in range(len(sersic_profiles_n)):
	plt.plot(rad, -2.5*np.log10(sersic_profiles_n[i]),
			ls = '-',lw=1., color=colors[i], label=r'$n=%.1f, \sigma=%.2f$' %(n_array[i],std[i]))

plt.plot(rad, -2.5*np.log10(sersic_profile_4),
			ls = '-',lw=1., color='b', label=r'$n=4.0, \sigma=%.2f$' %(std_n4) )

plt.plot(rad,mag,'ko', ms=1.5)
plt.xscale('log')
#plt.xlim((3e-1,200))
plt.xlim((2.,130))

plt.gca().invert_yaxis()
plt.ylim((28,14))
plt.xlabel(r'$R\,[\mathrm{arcsec}]$', size=fontsize)
plt.ylabel(r'$\mathrm{Surface\ brightness\ [mag/arcsec^2]}$', size=fontsize)
plt.legend(loc='upper right', ncol=3, prop = {'size':fontsize})
plt.text(70,17,r'$\mu_e=%.1f$' %(mu_e))

plt.savefig('deVaucouleurs_profiles_n4_zoom_mu_e_%d_std.pdf' %(int(10*mu_e)))
plt.show()
