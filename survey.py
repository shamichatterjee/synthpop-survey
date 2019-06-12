#! /usr/bin/env python
# Survey a synthetic population - SC 20190131

"""
  Population synthesis by Tyler Cohen + Paul Demorest 20190612
  Pop.shape is (37233 rows, 16 cols)
  Cols are {P Pdot DM t_scatt W_ms GL GB S1400 L1400 SPIDX Dist X Y Z sigtoa3 sigtoa10}
  Col IDs  {0 1    2  3       4    5  6  7     8     9     10 11 12 13 14     15}
  Note: X**2 + (Y-8.5)**2 + Z**2 = Dist**2
  Note: S1400 = L1400/Dist**2 ; in mJy
  Note: GL, GB are in degrees. GL=(-180,180) GB=(-90,90)
"""

import sys
from numpy import *
from matplotlib.pylab import *
from astropy import units as u
from astropy.coordinates import SkyCoord

kGL = 5
kGB = 6
kS1400 = 7
kX = 11
kY = 12


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
		print("Separations to item %4d of %4d" % (i, pcoord.size), end='\r')
		for j in range(i):
			angseps.append(float(pcoord[i].separation(pcoord[j])/u.degree))

	return(angseps)


def calcseps(psrpop):
	# What are the sampled angular separations?
	popgl = psrpop[:,kGL]
	popgb = psrpop[:,kGB]
	pcoord = SkyCoord(popgl, popgb, frame='galactic', unit='deg')

	# astropy.Separation for vectors is calculated on a per-element basis.
	angseps = []
	for i in range(1,pcoord.size):
		print("Separations to item %4d of %4d" % (i, pcoord.size), end='\r')
		# Create a skycoord list that is rotated step by step
		rotated = SkyCoord(roll(popgl,i), roll(popgb,i), frame='galactic', unit='deg')
		# Append list of separations as astropy Quantities
		angseps.append(ndarray.tolist(pcoord.separation(rotated)/u.degree))
		# Concatenate into one long list before returning
	print("")
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

def plotangsep(angseps, name='Survey separations', color='#b5323a'):
	hist(angseps, bins=45, label=name, color=color, rwidth=0.9)
	legend()
	xlabel("Angle (degrees)")
	ylabel("Counts")
	return(0)


def main():
	# Read in the entire population
	pop = loadtxt('dsa2kpop_sigtoa.asc')
	# pop = loadtxt('test_pop.asc')

	# Define the survey
	surveypars = defsurvey('AO')
	# Filter by survey parameters
	aopop = dosurvey(pop, surveypars)	# 433 detected
	aoseps = calcseps(aopop)

	#surveypars = defsurvey('GB')
	#gbpop = dosurvey(pop, surveypars)	# 315 detected
	#gbseps = calcseps(gbpop)

	#surveypars = defsurvey('DSA')
	#dsapop = dosurvey(pop, surveypars)  # 1904 detected
	#dsaseps = calcseps(dsapop)

	return(0)


if __name__ == '__main__':
    main()

