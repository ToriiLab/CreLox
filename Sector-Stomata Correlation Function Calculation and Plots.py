#Aggregates the stomata and random point distance histogram data, calculates correlation function, and plots.

import numpy as np;
import matplotlib.pyplot as pyplot;

#Load your sector and random point distance histograms. Array structure should be:

#Sector_Data is a 3xnxN array

#The first dimension is stomata distance histograms, random point histograms, and normalized ratio data.

#The second dimension is the length of the number of distance bins n

#The third dimension is the total number of distance histograms N

EPFL9_Sector_Data = np.load('EPFL9_sector_data.npy', encoding = 'latin1');
EPF1_Sector_Data = np.load('EPF1_sector_data.npy', encoding = 'latin1');
Control_Sector_Data = np.load('Control_sector_data.npy', encoding = 'latin1');

#Make sure to use the same bin intervals as in your histogramming

distance_fun = np.logspace(np.log10(15),np.log10(2500),15);
distance_bins = np.insert(distance_fun,0,0);

#The Sector-Stomata Correlation Function is calculated using two numbers:

#The numerator, which is the number of stomatal distances within a bin
#The denominator, which is the expected value of normalized random points within a bin

#The correlation function is then the (numerator/denominator) - 1

#Calculate the numerator and denominator, normalized for the number of points,
#and then calculate the correlation function


EPF1_num = (np.sum(EPF1_Sector_Data[0],axis=0)/np.sum(np.sum(EPF1_Sector_Data[0],axis=0))); #denominators in this line and below line are normalization factors
EPF1_den = (np.sum(EPF1_Sector_Data[1],axis=0)/np.sum(np.sum(EPF1_Sector_Data[1],axis=0)));
EPF1_y = EPF1_num/EPF1_den - 1;

EPFL9_num = (np.sum(EPFL9_Sector_Data[0],axis=0)/np.sum(np.sum(EPFL9_Sector_Data[0],axis=0)));
EPFL9_den = (np.sum(EPFL9_Sector_Data[1],axis=0)/np.sum(np.sum(EPFL9_Sector_Data[1],axis=0)));
EPFL9_y = EPFL9_num/EPFL9_den - 1;

Control_num = (np.sum(Control_Sector_Data[0],axis=0)/np.sum(np.sum(Control_Sector_Data[0],axis=0)));
Control_den = (np.sum(Control_Sector_Data[1],axis=0)/np.sum(np.sum(Control_Sector_Data[1],axis=0)));
Control_y = Control_num/Control_den - 1;

#Calculate confidence intervals either through standard error propagation
#or load bootstrap resampling data

total_epf1_err = np.std(np.sum(EPF1_Sector_Data[0],axis=1));
total_epfl9_err = np.std(np.sum(EPFL9_Sector_Data[0],axis=1));
total_Control_err = np.std(np.sum(Control_Sector_Data[0],axis=1));


individual_epf1_err = np.std(EPF1_Sector_Data[0],axis =0);
individual_epfl9_err = np.std(EPFL9_Sector_Data[0],axis =0);
individual_Control_err = np.std(Control_Sector_Data[0],axis =0);

epf1_std = 1.96*np.absolute(EPF1_y ) * np.sqrt( (individual_epf1_err/(np.mean(EPF1_Sector_Data[0],axis=0)+1))**2 + (total_epf1_err / np.mean(np.sum(EPF1_Sector_Data[0],axis=1)))**2 );
epfl9_std = 1.96*np.absolute(EPFL9_y ) * np.sqrt( (individual_epfl9_err/np.mean(EPFL9_Sector_Data[0],axis=0))**2 + (total_epfl9_err / np.mean(np.sum(EPFL9_Sector_Data[0],axis=1)))**2 );
Control_std = 1.96*np.absolute(Control_y ) * np.sqrt( (individual_Control_err/np.mean(Control_Sector_Data[0],axis=0))**2 + (total_Control_err / np.mean(np.sum(Control_Sector_Data[1],axis=0)))**2 );

epf1_err = epf1_std / np.sqrt(len(EPF1_Sector_Data[0]))
epfl9_err = epfl9_std / np.sqrt(len(EPFL9_Sector_Data[0]))
Control_err = Control_std / np.sqrt(len(Control_Sector_Data[0]))

#Upper and lower bounds of confidence interval

epf1_high = EPF1_y + epf1_err;
epf1_low = EPF1_y - epf1_err;

epfl9_high = EPFL9_y + epfl9_err;
epfl9_low = EPFL9_y - epfl9_err;

Control_high = Control_y + Control_err;
Control_low = Control_y - Control_err;

#Plot correlation functions with confidence intervals


fig, ax = pyplot.subplots();

ax.fill_between(distance_fun, epf1_high,epf1_low, color='turquoise', alpha = 0.75);
ax.plot(distance_fun, EPF1_y, color='darkblue',linestyle='--', marker='.', label='EPF1');
ax.fill_between(distance_fun, epfl9_high,epfl9_low, color='crimson', alpha = 0.6);
ax.plot(distance_fun, EPFL9_y, color='crimson',linestyle='--', marker='.', label='STOMAGEN');
ax.fill_between(distance_fun, Control_high,Control_low, color='darkslategrey', alpha = 0.5);
ax.plot(distance_fun, Control_y, color='k',linestyle='--',marker='.', label='Empty Vector');
ax.plot(np.arange(0,301), np.zeros(np.shape(np.arange(0,301))), 'k-',alpha=0.75);
ax.legend(loc='lower right')
pyplot.title('Stomata-Sector Correlation with 95% Confidence Interval');
pyplot.xlabel('Distance (microns)');
pyplot.ylabel('Correlation');
pyplot.xlim(0,300);
pyplot.ylim(-1.05,1);
pyplot.savefig('Sector Correlation Plot.eps');
pyplot.savefig('Sector Correlation Plot.pdf');

pyplot.show();
