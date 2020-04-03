# SPACE Pipeline

How to use the scripts provided to analyze stomatal density and stomatal correlation relative to a region of interest (in this case, a GFP-expressing mosaic sector).

## File Preqrequisites:

The data you analyze must be in .xlsx files. Each .xlsx file must contain the complete spatial data for one cotyledon. For completeness, the following requirements must be met for an .xlsx file:

1. A worksheet that has the XY spatial coordinates of a cotyledon outline titled Cotyledon Outline.
2. A worksheet that has the XY spatial coordinates of the stomatal distribution titled Stomatal Positions.
3. A worksheet for each region of interest (in this case, a mosaic sector) titled Sector 1, Sector 2, etc.

Consult the sample worksheet for reference of format.

## Multiple File Analysis Prerequisites:

To save time on analyzing each dataset individually, the scripts are written so that you can load and analyze multiple files via a for loop. To do this, your data set file directories are stored in a filelist that is read by the scripts.

Consult the sample filelist for reference of format.

## Stomatal Density Analysis

To analyze stomatal density, use the script Stomatal Density Calculations.py. This script is capable of analyzing three distinct regions:

1. The stomatal density of your sectors within the cotyledon.
2. The stomatal density inside an range extended beyond the sectors, excluding the area of the original sectors themselves.
3. The remainder of the cotyledon.

You can set the range you want to calculate. The output is three npy files that store your stomatal densities in the corresponding regions.

## Spatial Correlation Analysis

First use the script Stomata and Random Point Histograms.py to calculate and histogram the distances between stomata and sectors, and the distances between randomly generated point distrubitions and sectors. You can set a filter for upper and lower bound of sector areas to analyze if desired. The output is an array that stores your histograms.

Then load these histograms with Sector-Stomata Correlation Function Calculation and Plots.py. Correlation function is calculated and plotted.

## Sample Distribution Correlation

This file is to illustrate the concept of spatial correlation between a region of interest and a given point distribution. No files or prior input is needed to run, but it can be modified if you want to see what different point distributions look like. In the case of modification, though, you must be careful that your distance bins are appropriate for capturing the shape of your distribution, or your results will be unrepresentative.

## Systems

Scripts were written and run on Spyder.

## Acknowledgments

Majority of code was written by Scott Zeng, with initial foundation for data extraction and correlation analysis written by Emily Lo, in the Torii Lab. Primary contact on the corresponding paper is Keiko Torii.

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

