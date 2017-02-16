'''
## timematch.py
##
## matches times of one sparse dataset data2 to one full regular dataset data1 in a xx:xx:xx format
## assumes data2 has same length as times2
##
## df 2017-01-24
## last modified df 2017-01-25
'''

import numpy as np

def convertdecimal(timelis):
	''' convert time to decimal'''
	newlis = []
	for item in timelis:
		time  = item.split(':')
		timed = float(time[0]) + float(time[1])/60. + float(time[2])/(60.*60)
		newlis.append(timed)
	return newlis

def timematch(ti1, ti2, da1, da2):
	''' matching times and data '''
	ti1     = np.array(ti1)
	datapad = np.zeros_like(ti1)
	timepad = np.zeros_like(ti1)
	proffwhm= []
	for i, item in enumerate(ti2):
		## overwrites if there are two data2 points closest to the same data1 point 
		arg = np.argmin((ti1-item)**2)
		datapad[arg]   = da2[i]
		timepad[arg]   = item
		proffwhm.append(da1[i])
	return datapad, timepad, proffwhm

if __name__ == "__main__":
	''' reads data set '''

	## profiler times and data here
	times1 = np.loadtxt("times1", dtype='str')		
	data1  = np.loadtxt("data_mass") #data_mass

	## image times and data here
	times2 = np.loadtxt("times2", dtype='str') #times3
	data2  = np.loadtxt("80EE_closed_12") #data12open #80EE

	dtimes1 = convertdecimal(times1)
	dtimes2 = convertdecimal(times2)
	newdata, newtime, profdata = timematch(dtimes1, dtimes2, data1, data2)

	## outptuts array that is zeros except where there is data in data2 to compare with data1.
	## ouptputs decimal time for sanity check
	np.savetxt('output.txt', np.array([newdata, newtime]).T, newline='\n')

	## outputs values of data1 at times closest to data2
	np.savetxt('short.txt', np.array(profdata).T , newline='\n')

	
	
	
