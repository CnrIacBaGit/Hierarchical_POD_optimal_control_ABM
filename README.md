# Hierarchical_POD_optimal_control_ABM

Hierarchical clustering and dimensional reduction for optimal control of large-scale agent-based models

Angela Monti, Istituto per le Applicazioni del Calcolo "M. Picone" (IAC), CNR, Bari, Italy. Mail: angela.monti@cnr.it  
Fasma Diele, Istituto per le Applicazioni del Calcolo "M. Picone" (IAC), CNR, Bari, Italy  
Dante Kalise, Department of Mathematics, Imperial College London, UK  

This repository contains MATLAB routines (tested with version R2024b) for the simulation and optimal control of large-scale agent-based models using hierarchical clustering and model reduction techniques.

The implemented framework combines:
- hierarchical clustering strategies for agent aggregation,
- dimensionality reduction techniques,
- optimal control algorithms for large-scale dynamical systems.

The repository includes:

- script_full_uncontrolled.m : simulation of the full agent-based system without control
- script_full_controlled.m : optimal control applied to the full agent-based system  
- script_cluster.m : optimal control applied to the reduced dynamics by clustering
- script_pod.m : optimal control applied to the reduced model obtained by POD 
- script_complete.m : optimal control applied to the two-level reduced model (clustering + POD) 

Data files (.mat): datasets for reproducibility of numerical experiments.

The routines have been implemented and developed by Angela Monti and Fasma Diele. They can be used under the conditions of CC-BY-NC 2.0. When utilizing this codebase, please cite the following publication:

A. Monti, F. Diele, D. Kalise, Hierarchical clustering and dimensional reduction for optimal control of large-scale agent-based models, arXiv preprint, 2025, https://doi.org/10.48550/arXiv.2507.19644.

The complete description of the model, control framework, and numerical methods is available in the cited manuscript.

The development and implementation of the model and routines have been made possible thanks to the National Recovery and Resilience Plan (NRRP), Mission 4 Component 2 Investment 1.4—Call for tender No. 3138 of 16 December 2021, rectified by Decree No. 3175 of 18 December 2021 of the Italian Ministry of University and Research, funded by the European Union—NextGenerationEU; Award Number: Project code CN 00000033, Concession Decree No. 1034 of 17 June 2022 adopted by the Italian Ministry of University and Research, CUP B83C22002930006, Project title “National Biodiversity Future Center” (NBFC).

A. Monti acknowledges support from a scholarship (borsa per l’estero) granted by the Istituto Nazionale di Alta Matematica (INdAM) to carry out a research stay at Imperial College London.
