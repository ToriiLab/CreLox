#Generates sample point probablistic distributions relative to a sector
# and calculates the sector-point correlation function relative to a
#uniformly random point distribution.


import numpy as np;
import matplotlib.pyplot as pyplot;
import random;
import math;
import shapely.geometry as shg;

'''
Begin section for defining necessary functions.
'''

#All probability distributions and sector boundaries are radial distributions, but
#it is necessary to export them in x-y coordinates.
    
#Creates a Point Distrubtion along the circumference of a circle
#to serve as a Sector, where r is the radius and n is the number of points to create.
#Returns the radial coordinates in x-y form in a nx2 array.

def PointsInCircum(r,n):
    return np.asarray([(math.cos(2*math.pi/n*x)*r,math.sin(2*math.pi/n*x)*r) for x in range(0,n)])

#Creates a radial point distrubtion that is positively correlated at closer ranges,
#and linearly decreases in a triangular probability distribution.
#Points are randomly generated from this probability distribution, rather than uniformly along the range.
#R1 is the radius at which the probability distribution begins, R2 is the length of the distribution,
#and n is the number of points to generate from this probability distribution.
#The probability distribution is radial, but x and y coordinates are returned.
#Returns the radial coordinates in x-y form in a nx2 array.
    
def PositivePoints(radius1,radius2,n):
    sample_list = [];
    for i in range(0,n):
        theta = random.uniform(0,n);
        radius = radius1 + np.random.triangular(0,0,radius2);
        sample_list.append([math.cos(2*math.pi/n*theta)*radius,math.sin(2*math.pi/n*theta)*radius]);
    return np.asarray(sample_list);
    
#Create a uniform random distribution that should have zero correlation.
#R1 is the lower endpoint and R2 is the length of the distribution.
#N is the number of random points to generate based off this probability distribution.
#Returns the radial coordinates in x-y form in a nx2 array.
    
def ZeroCorrPoints(radius1,radius2,n):
    sample_list = [];
    for i in range(0,n):
        theta = random.uniform(0,n);
        radius = radius1 + np.random.uniform(0,radius2);
        sample_list.append([math.cos(2*math.pi/n*theta)*radius,math.sin(2*math.pi/n*theta)*radius]);
    return np.asarray(sample_list);    

#Creates the polygon defined by the sector outline and computes the nearest distance
#between each individual point of your sapmle distribution and the sector.
#points_list is the point distribution, and outline is the sector. 
#Returns a nx2 array where n is the length of points_list, which must also be an nx2 array.
    
def computeSectorDistance(points_list,outline):    
    point_count = len(points_list);
    poly = shg.Polygon(outline);
    distance_list = np.zeros(point_count)
    for i in range(0,point_count):
        point = shg.Point(points_list[i]);
        distance_list[i] = poly.distance(point);
    distance_list = np.asarray(distance_list);
    distances = distance_list[distance_list!=0];
    return distances;


'''
End of function definition section.
'''    
    

'''
Begin of procedural section.
'''


#Create Sample Distributions and Sectors

#Create your sector outline

sector = PointsInCircum(150,100);
sector = sector + 1250; #Shift the center of the sector

#Create a uniformly spaced set of points extending outward from the sector circle
#to demonstrate regions of negative correlation

negative_sample = [];

for i in np.arange(0,5):
    ring = PointsInCircum(250+i*200,100);
    ring = ring+1250;
    negative_sample.append(ring);
    
sector = np.asarray(sector);

#Create sample distribution with close ranged clusters to demonstrate positive correlation


positive_sample = PositivePoints(150,1050,500);
positive_sample = positive_sample+1250


#Create sample distribution with uniform probability distribution

zero_sample = ZeroCorrPoints(150,1050,500);
zero_sample = zero_sample+1250;

#Plot the distributions

#First distribution

for i in np.arange(0,5):
    pyplot.scatter(negative_sample[i][:,0],negative_sample[i][:,1],marker='.',color='k');

pyplot.scatter(sector[:,0],sector[:,1],marker='.',color='g');
pyplot.axis('equal');
pyplot.xlim(0,2500);
pyplot.ylim(0,2500);
pyplot.savefig('Inhibition Plot.eps');
pyplot.show(); 

#Second distribution

pyplot.plot(positive_sample[:,0], positive_sample[:,1], 'k.');
pyplot.scatter(sector[:,0],sector[:,1],marker='.',color='g');
pyplot.axis('equal');
pyplot.xlim(0,2500);
pyplot.ylim(0,2500);
pyplot.savefig('Clustered Plot.eps');
pyplot.show(); 

#Third distribution

pyplot.plot(zero_sample[:,0], zero_sample[:,1], 'k.');
pyplot.scatter(sector[:,0],sector[:,1],marker='.',color='g');
pyplot.axis('equal');
pyplot.xlim(0,2500);
pyplot.ylim(0,2500);
pyplot.savefig('Zero Plot.eps');
pyplot.show(); 

#Generate a set of random distributions, where Ntrials is the number of distributions
#generated independently

Ntrials = 250;

random_sets = []

for i in range(0,Ntrials):
    trial = ZeroCorrPoints(150,1050,500);
    trial = trial+1250;
    random_sets.append(trial);

random_sets = np.asarray(random_sets);
    
#Calculate the distances between individual points in each distribution and the sector.
#Then calculate the correlation function estimator and plot the result.

random_distances = [];

distance_bins = np.arange(-49,1052,100); #Set bins and function for histogramming
distance_fun = np.delete(distance_bins,0);

#For Correlation of First Distribution

negative_sample = np.reshape(negative_sample, (500,2))

sector_distances = computeSectorDistance(negative_sample,sector);

random_distances = [];

for i in range(0,Ntrials):
    rand_dist = computeSectorDistance(random_sets[i],sector);
    random_distances.append(rand_dist);
    
random_distances = np.asarray(random_distances);
random_distances = np.reshape(random_distances, 500*Ntrials);

sector_counts,bin_data,etc = pyplot.hist(sector_distances, bins = distance_bins, alpha = 0.5, color = 'b',histtype='step');    
random_counts,bin_data,etc = pyplot.hist(random_distances, bins = distance_bins, alpha = 0.5, color = 'b',histtype='step');  
  
correlation = sector_counts/(random_counts/Ntrials) - 1; #dividing by Ntrials is for normalization

fig = pyplot.figure();

pyplot.plot(distance_fun,correlation,'k.-');
pyplot.xlabel('Distance (microns)');
pyplot.ylabel('Correlation');
pyplot.title('Inhibitory Correlation');
pyplot.xlim(0,1200);
pyplot.ylim(-1.5,1.5);
pyplot.axhline(0)
pyplot.savefig('Inhibition Correlation.eps');
pyplot.show();

#For Correlation of Second Distribution

sector_distances = computeSectorDistance(positive_sample,sector);

random_distances = [];

sector_counts,bin_data,etc = pyplot.hist(sector_distances, bins = distance_bins, alpha = 0.5, color = 'b',histtype='step');    
  
correlation = sector_counts/(random_counts/Ntrials) - 1; #dividing by Ntrials is for normalization

fig = pyplot.figure();

pyplot.plot(distance_fun,correlation,'k.-');
pyplot.xlabel('Distance (microns)');
pyplot.ylabel('Correlation');
pyplot.title('Clustered Correlation');
pyplot.xlim(0,1200);
pyplot.ylim(-1.5,1.5);
pyplot.axhline(0)
pyplot.savefig('Clustered Correlation.eps');
pyplot.show();

#For Correlation of Second Distribution

sector_distances = computeSectorDistance(zero_sample,sector);

random_distances = [];

sector_counts,bin_data,etc = pyplot.hist(sector_distances, bins = distance_bins, alpha = 0.5, color = 'b',histtype='step');    
  
correlation = sector_counts/(random_counts/Ntrials) - 1; #dividing by Ntrials is for normalization

fig = pyplot.figure();

pyplot.plot(distance_fun,correlation,'k.-');
pyplot.xlabel('Distance (microns)');
pyplot.ylabel('Correlation');
pyplot.title('Zero Correlation');
pyplot.xlim(0,1200);
pyplot.ylim(-1.5,1.5);
pyplot.axhline(0)
pyplot.savefig('Zero Correlation.eps');
pyplot.show();