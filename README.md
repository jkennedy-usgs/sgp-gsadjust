[![Build Status](https://travis-ci.org/jkennedy-usgs/sgp-gsadjust.svg?branch=master)](https://travis-ci.org/jkennedy-usgs/sgp-gsadjust)
[![Build status](https://ci.appveyor.com/api/projects/status/mb5yw6ksgwmfu18d?svg=true)](https://ci.appveyor.com/project/jkennedy-usgs/sgp-gsadjust)

**This is the development branch of GSadjust. Code in this repository has not necessarily been reviewed by USGS.**

**The approved preliminary software release is at** [https://github.com/usgs/sgp-gsadjust](https://github.com/usgs/sgp-gsadjust "USGS GSadjust repo")

# GSadjust
GSadjust is a cross-platform graphical interface for processing relative gravity surveys. It provides an interface for data selection, drift evaluation and correction, network adjustment, and integrating data from modern relative (Scintrex, ZLS) and absolute (Micro-g Lacoste) gravity meters. This software is a product of the USGS Southwest Gravity Program (http://go.usa.gov/xqBnQ).

GSadjust is public-domain software written for Python 3.5 and the PyQt5 framework. Although basic functionality has been tested on Mac and Linux machines, there may be unknown issues on those operating systems. 

Bug reports by email (jkennedy@usgs.gov) or by submitting an issue through GitHub are welcome. 

**Online help is available at [https://jkennedy-usgs.github.io/sgp-gsadjust/](https://jkennedy-usgs.github.io/sgp-gsadjust/)**

**Installing GSadjust**

**NEW! (2019-03-25):** Windows and Mac installers are [now available](https://github.com/jkennedy-usgs/sgp-gsadjust/releases/tag/v0.9-alpha). This installs a standalone Python environment with everything needed for GSadjust. On Windows, GSadjust is added to the Start menu and optionally to the Desktop. On Mac, after mounting the .dmg, GSadjust can be copied to the Applications folder or elsewhere. The Python installation should not interfere with other previously installed Python environements. Admin rights are NOT required on Windows if installed to the default user directory (i.e., anywhere other than the C:\Program Files directory).

On Linux, the most reliable way to install GSadjust is to use Conda to create a custom Python environment (these instructions can also be used on Windows or Mac):

1. Download and install the appropriate Miniconda installer for your platform (Windows, Mac, Linux, 32 or (probably) 64-bit):
[https://conda.io/miniconda.html](https://conda.io/miniconda.html "Miniconda download"). Any Python version (e.g., 2.7 or 3.6) is okay. When installing on Windows, select the "Add Anaconda to my PATH environment variable" option, or, in step 3, launch the command window by choosing "Anaconda (64-bit)" from the Windows Start Menu. Alternatively, the Miniconda executables (i.e., conda.exe) can be manually added to the Windows path.

2. Download the GSadjust repository as a zip file using the green "Clone or download" button near the top of this page. Unzip the directory in a location of your choice.

3. Open a terminal in the new directory (sgp-gsadjust or sgp-gsadjust-master) and create a new conda environment:

    ```
    conda env create --file environment.yml
    ```

    This creates a new conda environment named py35 with the appropriate python modules. 
    
    Alternatively, install Python 3.5+ directly and/or in a different virtual environment. The following packages are required:
    
    ```
    matplotlib numpy pyqt=5 scipy networkx

    ```
    
4a. If Miniconda (conda.exe) is installed on your path, double-click gsadjust.bat to launch GSadjust.

4b. If Miniconda is _not_ on your path, open the Anaconda command window (or re-use the command window from step 3), navigate to the sgp-gsadjust-master directory, and type:

    ```
    activate py35
    python gsadjust\gsadjust.py
    ```

## Data selection interface
* Select/deselect samples and stations
* Plot gravity, earth tide correction, tilt, etc.
* Imports both Scintrex and Burris data formats

![data tab](https://cloud.githubusercontent.com/assets/9754449/22989620/d5923194-f373-11e6-80a6-1d106a1623f7.png)

## Drift correction
* Allows for Roman (1946) method, continuous drift correction (Kennedy and Ferre, 2016, doi: 10.1093/gji/ggv493), or determined from least-squares solution
* Graphical interface shows drift rate vs. time

![drift_tab](https://user-images.githubusercontent.com/9754449/38323886-6d0c9a68-37f3-11e8-96ab-17f7cf4eff05.PNG)

## Network adjustment
* Individual delta-g's can be enabled/diabled
* Uses numpy or Gravnet (Hwang and others, 2002) for network adjustment
* Import absolute-gravity data directly from .project.txt files
* Adjustment diagnostics: residual histogram, measured vs. adjusted, etc.

![adjust_tab](https://cloud.githubusercontent.com/assets/9754449/22990492/9d34be7c-f376-11e6-92e9-f3b3e85a808c.png)

## Gravity change
* Calculates gravity change between series
* Shows time-series plot of gravity change at all stations
* Shows map view with user-selectable base map

![gravity_change](https://user-images.githubusercontent.com/9754449/56771106-9b3a6300-676a-11e9-95a6-12fd0da14810.png)
---
Originally based on 

PyGrav
Version 1.0  
21 Jan. 2016  
[pyGrav website](http://basileh.github.io/pyGrav/)

Hector, B. (1,2) and Hinderer, J.(1): *pyGrav, a Python-based program for handling and 
processing relative gravity data*, Computers & Geosciences, [doi:10.1016/j.cageo.2016.03.010](http://dx.doi.org/10.1016/j.cageo.2016.03.010), 2016.

(1): Institut de Physique du Globe de Strasbourg  UMR 7516 CNRS/Université de Strasbourg, 5 rue Descartes 67084 Strasbourg, France

(2): CNES/IRD/UJF-Grenoble 1/CNRS/G-INP, LTHE, UMR 5564, 38041 Grenoble, France

---
For convenience GSadjust is distributed with a Windows executable for Gravnet software (Hwang and others, 2002; not to be confused with Gravnet software used by the National Geospatial-Intelligence Agency). USGS assumes no responsibility for the accuracy of the network-adjustment results generated by Gravnet. 

Hwang, C., Wang, C., & Lee, L. (2002). Adjustment of relative gravity measurements using weighted and datum-free constraints. Computers & Geosciences, 28, 1005�1015.

---
GSadjust is public domain software. No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material. Please refer to LICENSE.md and DISCLAIMER.md for complete information.

