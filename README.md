# Automated Fracture Detector (Modified Version). 2025
# Authors: Losarcos, J. M., Bernardi, M. I.,and Hryb D. E.
# Corresponding author: Juan Manuel Losarcos, University of Illinois Chicago (jlosa@uic.edu)  

This repository contains a single MATLAB script (AutomatedFractureDetector_v1_1.m) that implements an automated method for extracting digitized fractures from images of fractured rocks using the complex shearlet transform. The method is adapted from the manuscript titled *"An automated fracture trace detection technique using the complex shearlet transform", Prabhakaran et al. (2019).
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16696176.svg)](https://doi.org/10.5281/zenodo.16696176)


---

## Contents

- **AutomatedFractureDetector_v1_1.m**  
  Main script that runs the full fracture extraction workflow:  
  → Generates shearlet-based ridge ensemble from input rock images and parameters selection 
  → Applies sigmoid thresholding to filter probable fractures  
  → Segments and skeletonizes ridges  
  → Fits polylines to fracture traces  
  → Computes orientation and length statistics  
  → Exports results to Excel and generates rose diagrams and histograms  

  This version includes additional interactive tools for scale calibration and north orientation correction.

---

## Folder Structure

- **StudyCases/**  
  Contains the three original images used in the study at different scales: local scale, outcrop scale, and regional scale.

- **Dependencies/**  
  Contains all required MATLAB dependencies and precompiled MEX files, including:

  - CoSHREM Toolbox: Please refer to : "Reisenhofer, R., Kiefer, J., and King, E. J.: Shearlet-based         detection of flame fronts, Experiments in Fluids, 57, 41, https://doi.org/10.1007/s00348-016-2128-6,     2016" for a detailed description of the CoSHREM toolbox.  
    (http://www.math.uni-bremen.de/cda/software.html)  
  - Geom2D Toolbox  
    (https://www.mathworks.com/matlabcentral/fileexchange/7844-geom2d)  
  - Douglas-Peucker algorithm  
    (https://www.mathworks.com/matlabcentral/fileexchange/61046-douglas-peucker-algorithm)  
  - ShearLab 3D functions  
    (http://www.shearlab.org/)  
  - yapuls.m from YAWTb  
    (http://sites.uclouvain.be/ispgroup/yawtb/)

---

## Notes on Modifications

This script consolidates all steps into a single executable file, with simplified outputs and no dependency on shapefile generation. Additional interactive modules allow users to calibrate real-world scale and magnetic north for each image.

---

## License

SPDX-License-Identifier: GPL-3.0-only
Copyright (c) 2025 Losarcos, J. M.; Bernardi, M. I.; and Hryb, D. E.
Portions Copyright (c) 2019 Rahul Prabhakaran, TU Delft
This program is free software: you can redistribute it and/or modify it
under the terms of the GNU GPL v3.0 as published by the Free Software Foundation.
See the LICENSE file for full text. Distributed WITHOUT ANY WARRANTY.

---

## For more information about the full functionality of the original multi-script workflow, please refer to the original manuscript:

## Prabhakaran, R., Bruna, P.-O., Bertotti, G., & Smeulders, D. (2019). An automated fracture trace detection technique using the complex shearlet transform. Solid Earth, 10(6), 2137–2166. https://doi.org/10.5194/se-10-2137-2019
