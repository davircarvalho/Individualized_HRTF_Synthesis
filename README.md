# EAC - Individualized HRTF Synthesis

<p align="left">
  <a href="https://github.com/davircarvalho/Individualized_HRTF_Synthesis/releases/" target="_blank">
    <img alt="GitHub release" src="https://img.shields.io/github/v/release/davircarvalho/Individualized_HRTF_Synthesis?include_prereleases&style=flat-square">
  </a>

  <a href="https://github.com/davircarvalho/Individualized_HRTF_Synthesis/commits/master" target="_blank">
    <img src="https://img.shields.io/github/last-commit/davircarvalho/Individualized_HRTF_Synthesis?style=flat-square" alt="GitHub last commit">
  </a>

  <a href="https://github.com/davircarvalho/Individualized_HRTF_Synthesis/issues" target="_blank">
    <img src="https://img.shields.io/github/issues/davircarvalho/Individualized_HRTF_Synthesis?style=flat-square&color=red" alt="GitHub issues">
  </a>

  <a href="https://github.com/davircarvalho/Individualized_HRTF_Synthesis/blob/master/LICENSE" target="_blank">
    <img alt="LICENSE" src="https://img.shields.io/github/license/davircarvalho/Individualized_HRTF_Synthesis?style=flat-square&color=yellow">
  <a/>

</p>
<hr>



**Synthesis of individualized HRTFs based on Neural Networks, Principal Component Analysis and anthropometry.**


*This repository was developed as part of a bachelor thesis project for Acoustic Engineering (EAC) @ Federal University of Santa Maria - Brazil.*

The thesis can be found [here](https://drive.google.com/file/d/1JVDNxQreYzg7jfauMwFnK3Sg21aB5ri5/view?usp=sharing)

## Individualized HRTF app 

MATLAB app that can be used to generate individualized HRTFs with SOFA or HeSuVi extensions, using the proposed model.


## Auraliza app

MATLAB app to create real-time virtual auditory scenes using SOFA SimpleFreeField HRTFs and n-channel audio inputs.

Latest updates for this app can be found [here](https://github.com/davircarvalho/Auralization_Engine)


## Main - (PCA based)

Contains the pre-processing and post-processing, regression model and HRTF reconstruction routines.


> *It may be necessary to ajust paths for your local directories.*


## Functions 

Toolbox to work with SOFA HRTFs and general functions used by the main routines.
