#! /usr/bin/env python
# Survey a synthetic population - SC 20190131

"""
  dsa2kpop_sigtoa.asc:  Population synthesis by Tyler Cohen + Paul Demorest 20190612

  Pop.shape is (37233 rows, 16 cols)
  Cols are {P Pdot DM t_scatt W_ms GL GB S1400 L1400 SPIDX Dist X Y Z sigtoa3 sigtoa10}
  Col IDs  {0 1    2  3       4    5  6  7     8     9     10 11 12 13 14     15}

  # P in ms
  # X**2 + (Y-8.5)**2 + Z**2 = Dist**2
  # S1400 = L1400/Dist**2 ; in mJy
  # GL, GB are in degrees. GL=(-180,180) GB=(-90,90)

  test_pop.asc:    1% random sample of das2kpop for tests
  NG12p5_pop.asc:  NANOGrav 12.5yr pulsars in same format, Michael Lam 20190615
"""

import sys
from numpy import *
from matplotlib.pylab import *
from astropy import units as u
from astropy.coordinates import SkyCoord
from matplotlib import gridspec

matplotlib.pylab.rcParams['text.usetex'] = True
matplotlib.pylab.rcParams['font.family'] = 'serif'

# Change 10^0 and 10^1 axis labels, via MTL 20190709
from matplotlib.ticker import FuncFormatter

def nolog(x,pos):
    y = np.log10(x)
    if y == 1:
        return "$\hfill 10$"
    elif y == 0:
        return "$\hfill 1$"
    elif y == -1:
        return "$\hfill 0.1$"
    elif y == -2:
        return "$\hfill 0.01$"
    return "$\hfill 10^{%i}$" % np.log10(x)

myformatter = FuncFormatter(nolog)


# column IDs
kGL = 5
kGB = 6
kS1400 = 7
kX = 11
kY = 12

# filenames
fsimpop = 'dsa2kpop_sigtoa.asc'
ftestpop = 'test_pop.asc'
fnanopop = 'NG12p5_pop.asc'

def defsurvey(tel):
	# Survey parameters for notional GB, AO, DSA surveys
	aopars = {}
	aopars['Name'] = 'Arecibo'
	aopars['Dec1'] = 0.0
	aopars['Dec2'] = 35.0
	aopars['Smin'] = 0.02   # PSRCAT empirical, 20 uJy

	gbpars = {}
	gbpars['Name'] = 'GBT'
	gbpars['Dec1'] = -30.0
	gbpars['Dec2'] = 90.0
	gbpars['Smin'] = 0.15   # PSRCAT empirical, 150 uJy

	dsapars = {}
	dsapars['Name'] = 'DSA 2000'
	dsapars['Dec1'] = -30.0
	dsapars['Dec2'] = 85.0
	dsapars['Smin'] = 0.02

	if (tel=='AO'):
		return(aopars)
	elif (tel=='GB'):
		return(gbpars)
	elif (tel=='DSA'):
		return(dsapars)
	else:
		raise Exception('Unknown survey')


def dosurvey(pop, surveypars):
	# Use surveypars to select from synthetic population
	print("Selecting for survey %s" % (surveypars['Name']))

	popgl = pop[:,kGL]
	popgb = pop[:,kGB]
	popcoord = SkyCoord(popgl, popgb, frame='galactic', unit='deg')

	declim1 = surveypars['Dec1']*u.degree
	declim2 = surveypars['Dec2']*u.degree
	is_in_sky = where(logical_and(
		(popcoord.icrs.dec>declim1),
		(popcoord.icrs.dec<declim2)
		))
	in_sky = pop[is_in_sky]
	print("%d pulsars in Dec range" % len(in_sky))

	sky_flux = in_sky[:,kS1400]
	bright = where(sky_flux > surveypars['Smin'])
	detected = in_sky[bright]
	print("%d pulsars detected" % len(detected))

	return(detected)


def calcseps_crude(psrpop):
	# What are the sampled angular separations?
	popgl = psrpop[:,kGL]
	popgb = psrpop[:,kGB]
	pcoord = SkyCoord(popgl, popgb, frame='galactic', unit='deg')

	# astropy.Separation for vectors is calculated on a per-element basis.
	# Just brute-force iterate?
	angseps = []
	for i in range(pcoord.size):
		print("Separations to item %4d of %4d" % (i+1, pcoord.size), end='\r')
		for j in range(i):
			angseps.append(float(pcoord[i].separation(pcoord[j])/u.degree))

	print("Proper count; use weight=1.0 (non-default) in plotangsep")
	return(array(angseps))


def calcseps(psrpop):
	# What are the sampled angular separations?
	popgl = psrpop[:,kGL]
	popgb = psrpop[:,kGB]
	pcoord = SkyCoord(popgl, popgb, frame='galactic', unit='deg')

	# astropy.Separation for vectors is calculated on a per-element basis.
	angseps = []
	for i in range(1,pcoord.size):
		print("Separations to item %4d of %4d" % (i+1, pcoord.size), end='\r')
		# Create a skycoord list that is rotated step by step
		rotated = SkyCoord(roll(popgl,i), roll(popgb,i), frame='galactic', unit='deg')
		# Append list of separations as astropy Quantities
		angseps.append(ndarray.tolist(pcoord.separation(rotated)/u.degree))
		# Concatenate into one long list before returning
	print("Double counted! Use weight = 0.5 (default) in plotangsep")
	return(concatenate(angseps))


def Hammer(GLdeg,GBdeg):
	# Hammer projection
	lrad = GLdeg * pi/180.0
	brad = GBdeg * pi/180.0
	dterm = sqrt((1+ cos(brad)*cos(lrad/2))/2.0)
	xHammer = 2.0*cos(brad)*sin(lrad/2)/dterm * 180.0/pi
	yHammer = sin(brad)/dterm * 180.0/pi
	return(xHammer,yHammer)


def plotxy(pop_detected, pop_all='', name='Survey', color='orange'):
	# Assumes figure is set up and axes are available
	if(pop_all!=''):
		plot(pop_all[:,kX], pop_all[:,kY], 'b,')
	plot(pop_detected[:,kX], pop_detected[:,kY],
		color=color, marker='*', linestyle='', label=name)
	xlabel('X (kpc)')
	ylabel('Y (kpc)')
	legend()
	return(0)

def plotlb(pop_detected, pop_all='', name='Survey', color='orange'):
	# Assumes figure is set up and axes are available
	if(pop_all!=''):
		px,py = Hammer(pop_all[:,kGL], pop_all[:,kGB])
		plot(px, py, 'b,')
	pdx,pdy = Hammer(pop_detected[:,kGL], pop_detected[:,kGB])
	plot(pdx, pdy, color=color, marker='*', linestyle='', label=name)
	xlabel('Gal_l (deg)')
	ylabel('Gal_b (deg)')
	legend()
	return(0)

def plotangsep(angseps, name='Survey separations', color='#b5323a', weight=0.5):
	weights = ones(angseps.size)*weight
	hist(angseps, bins=45, label=name, color=color, rwidth=0.9, weights=weights, log=True)
	legend()
	#xlabel("Angle (degrees)")
	#ylabel("Counts")
	return(0)

def genHDcurve(thetadeg):
	# Return an array with HD curve for 0...180 degrees
	# x = (1-cos theta)/2
	# HD = 3/2 x ln x - x/4 + 1/2 + 1/2 delta(x)
	theta = thetadeg*math.pi/180.0
	x = (1 - cos(theta))/2.0
	HD = 1.5*x*log(x) - 0.25*x + 0.5
	return(HD)


def makeA2020plot():
	theta = linspace(0.0001, 180.000, 900)
	hd = genHDcurve(theta)

	pngao = loadtxt('NGAO.asc')
	pnggb = loadtxt('NGGB.asc')
	pnano = loadtxt(fnanopop)
	psynth = loadtxt('dummy200.txt')

	s4 = calcseps(psynth)
	s3 = calcseps(pnano)
	s2 = calcseps(pngao)
	s1 = calcseps(pnggb)

	gs = gridspec.GridSpec(2,1, height_ratios=[1,3])

	ax2 = subplot(gs[1])
	plotangsep(s4, name='200 pulsars (synthetic)', color='#0B486B')
	plotangsep(s3, name='NANOGrav 12.5 yr', color='#3B8686')
	plotangsep(s2, name='NANOGrav Arecibo-only', color='#79BD9A')
	plotangsep(s1, name='NANOGrav GBT-only', color='#CFF09E')

	#ax2.set_xticklabels(arange(15,190,15))
	ax2.yaxis.set_major_formatter(myformatter)
	ax2.set_ylabel('Counts', fontsize=16)
	ax2.set_xlabel("Angle (degrees)", fontsize=16)

	ax1 = subplot(gs[0], sharex=ax2)
	ax1.plot(theta, hd, ls='-', color='#0B486B', label='Hellings-Downs correlation')
	ax1.legend(fontsize=14)
	ax1.plot(theta, zeros(theta.size), ls='--', color='#b5323a')
	ax1.set_ylabel('Correlation', fontsize=16)
	# ax1.set_xticklabels('')
	ax1.tick_params(axis='x', labelbottom=False)

	subplots_adjust(hspace=.0)
	subplots_adjust(top=0.95, bottom=0.11, left=0.11, right=0.95)


def main():

	print('''
	# Read in the population [ fsimpop | ftestpop | fnanopop ]
	pop = loadtxt(fnanopop)

	# Define the survey [ AO | GB | DSA ]
	surveypars = defsurvey('AO')

	# Filter by survey parameters
	foundpop = dosurvey(pop, surveypars)	# 433 detected

	# Calculate and plot separations
	seps = calcseps(foundpop)
	plotangsep(seps)
	''')

	return(0)

if __name__ == '__main__':
    main()

