# lightIonSummerSchool_forStudents_2025

## Introduction

This course is an interactive session to learn about how to perform a simple EIC data analysis. The complicated part that associated with the ePIC experiment has been taken care of, where the pseudo data file has been processed into a light weight flat root tree. 

The pseudo data is based on Sartre coherent VM simulation and BeAGLE incoherent VM simulation. The beam momentum is 18 x 137.5 GeV with electron scattering off an ion. 

Due to the Interaction Region design at ePIC (IP6 at the EIC), there is a 25 mrad crossing angle of the hadron beam. In addition, realistic beam effects, e.g., angular divergence and momentum spread, have been implemented to both beams. 

At the end of exercise, one should be able to perform a basic imaging measurement at the EIC using ePIC detector.

Pre-requisite:
* basic ROOT macros

### Get the repo:

```git clone https://github.com/KongTu/lightIonSummerSchool_forStudents_2025.git```

```cd lightIonSummerSchool_forStudents_2025```

### Download the pseudo data from google drive:
https://drive.google.com/file/d/1ofgtVJOQZuD-fsEqRfIrLZhoGxKIj9dp/view?usp=sharing

### Running the macros (in order):

``` root -l analyze.C+```

``` root -l getImage.C```

### Understanding the tree

#### mcp 
* MC particles that have the information of px, py, pz, mass, and status. Stable particles in the final-state have a status = 1.

#### particles
* Reconstructed charged particles in main detector and B0.

#### clusters_eemc
* Reconstructed cluster in backward EM calorimeter

#### clusters_zdc
* Reconstructed cluster in zdc

#### hit_rp
* reconstructed hits at RP but no tracking has been done.

#### hit_omd
* reconstructed hits at omd but no tracking has been done.

### Understanding the example macros (with basic structures)

#### analyze.C 
* Use this to process all data (loop over events, particles, etc.) and save histograms (results)

#### getImage.C
* use this to perform Fourier Transformation from "t" to "b"


