import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import cm, rcParams,rc
# number 3a)

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

mu_e = 20.6  # mag/arcsec^2
#R_e  = 30.   # arcsec
n = 4

# Define arrays
R_array = np.logspace(-2.,3.,500)
#n_array = np.array([0.2,0.6,1.,2.,3.,4.,6.,8.,10.])
R_e_array = np.array([25.,30.,35.])
# Calculate sersic profiles for different Re
sersic_profiles_Re = np.array([sersic_profile(R_array, n, mu_e, R_e) for R_e in R_e_array])

# Make the plot
figsizex = 5
f = plt.figure( figsize=(figsizex,figsizex*0.6), dpi=200 )

# define colors for lines
jet    = cm.get_cmap('jet')
colors = jet(np.linspace(0.7, 1, len(R_e_array)))

plt.axes( [0.1,0.15,0.85,0.75] ) # [x,y,size_x, size_y] 
for i in range(len(sersic_profiles_Re)):
	plt.plot(R_array, -2.5*np.log10(sersic_profiles_Re[i]),
			ls = '-',lw=1., color=colors[i], label=r'$R_e=%.1f$' %(R_e_array[i]))

plt.plot(rad,mag,'ko', ms=1.5)
plt.xscale('log')
plt.xlim((3e-1,200))

plt.gca().invert_yaxis()
plt.ylim((28,14))
plt.xlabel(r'$R\,[\mathrm{arcsec}]$', size=fontsize)
plt.ylabel(r'$\mathrm{Surface\ brightness\ [mag/arcsec^2]}$', size=fontsize)
plt.legend(loc='upper right', ncol=3, prop = {'size':fontsize})
plt.text(70,17,r'$\mu_e=%.1f$' %(mu_e))
plt.savefig('deVaucouleurs_profiles_mu_e_%d.pdf' %(int(10*mu_e)))
plt.show()
