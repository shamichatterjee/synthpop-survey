from numpy import *
import matplotlib.pylab as plt
from numpy.random import random_sample

plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'


col1 = '#3b4aab'
col2 = '#3a69a6'
col3 = '#3a86a1'
col4 = '#3a9c99'

col1t = '#3b4aab55'
col2t = '#3a69a640'
col3t = '#3a86a125'
col4t = '#3a9c9910'

psrdates = loadtxt('starts.txt', usecols=(2,))
lastdate = int(sorted(psrdates)[-1])+2

known = 25 # 25 known pulsars, not being timed
knowndates = concatenate((psrdates, 2020.0*ones(known)))

numpsr = len(psrdates) + known
addnum = 200 - numpsr  # Want to get to 200 pulsars by 2030, have numpsr already
newdates = (2026.0 - lastdate)*random_sample(addnum) + lastdate

alldates = concatenate((knowndates, newdates))

bins = arange(2004,2035,1)

def datarelease(ax, yr, num, str):
	boxht = 180
#	ax.arrow(yr, num+40, 0, -40, head_width=0.5, head_length=2.5, length_includes_head=True)
	ax.arrow(yr-0.2, boxht, 0, num-boxht, head_width=0.5, head_length=2.5, length_includes_head=True)
#	ax.text(yr, num+43, str, rotation='vertical', va='bottom', ha='center', bbox=dict(facecolor='white', alpha=1.0))
	ax.text(yr-0.2, boxht-50, str, rotation='vertical', va='bottom', ha='center', bbox=dict(facecolor='white', alpha=1.0))
	return(0)


fig = plt.figure()
ax = fig.add_subplot(1,1,1)

ax.hist(alldates, bins=bins, cumulative=True, histtype='stepfilled', color=col4t)
ax.hist(knowndates, bins=bins, cumulative=True, histtype='stepfilled', color=col2t)
ax.hist(psrdates, bins=bins, cumulative=True, histtype='stepfilled', color=col1t)

ax.hist(alldates, bins=bins, cumulative=True, histtype='step', color=col4, label='Future goal')
ax.hist(knowndates, bins=bins, cumulative=True, histtype='step', color=col2, label='Already known')
ax.hist(psrdates, bins=bins, cumulative=True, histtype='step', color=col1, label='Currently timed')

ax.set_xlabel('Year')
ax.set_ylabel('Number of MSPs timed + known + future')
# ax.legend(loc=2)

datarelease(ax, 2010, 17, r'5 yr')
datarelease(ax, 2013.5, 37, r'9 yr')
datarelease(ax, 2016, 45, r'11 yr')
datarelease(ax, 2017.5, 48, r'12.5 yr')
datarelease(ax, 2019, 67, r'14 yr')

ax.text(2019, 180, 'NANOGrav Data Releases', ha='right', va='bottom', bbox=dict(facecolor='white', alpha=1.0))

ax.text(2033, 25, r'{\Large Timed by NANOGrav}', ha='right', bbox=dict(facecolor='white', alpha=0.7, ec='None'))
ax.text(2033, len(psrdates)+known/2.0, r'{\large Already known}', ha='right', va='center', bbox=dict(facecolor='white', alpha=0.7, ec='None'))
ax.text(2033, numpsr+25, r'{\large Future surveys}', ha='right', bbox=dict(facecolor='white', alpha=0.7, ec='None'))


# savefig('obsbydate.pdf')
