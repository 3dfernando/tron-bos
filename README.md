# TRON-BOS
TRON-BOS (Telecentric Rotating Nozzle Background Oriented Schlieren) is a technique to obtain volumetric density fields through parallel-ray density tomography. 

# Pre-Requisites
This code currently requires:
* Matlab R2022a or later (likely works in earlier versions)
* The OSMODI solver code [(See here our open-source project)](https://github.com/3dfernando/pressure-osmosis)
  * This solver is GPU-native (Requires CUDA GPU Computing Toolkit and Matlab maxcuda capability if a version for your OS is not available)
  * You can compile the CPU-only solver OSMODI_CPU if you have trouble with CUDA (it may take you a few hours to figure it out!)
* A code to track the displacement of your speckle pattern (a PIV code/program can also be used) 

# Publications
Please refer to the following publications to understand the experimental setup and details:
* Zigunov (2025) ["Time-averaged density tomography of non-axisymmetric jets with a rotating nozzle and telecentric, single-camera BOS"](https://link.springer.com/article/10.1007/s00348-025-04052-7)
* Zigunov, Wagner and Namatsaliuk (2026) ["High-resolution, quantitative density tomography using a rotating nozzle"](https://drive.google.com/file/d/1jYes9pHc0Qot5JNKzwqnKZXtaBotmhhR/view?usp=sharing)

Considering citing the above publications if you use this code or setup ideas in your work.

# Workflow
The figure below summarizes the workflow within this technique:
![alt text](https://github.com/3dfernando/tron-bos/blob/main/img/Workflow.PNG "Workflow")


# Sample Files
Please check our Zenodo Repository for a set of sample raw data used with this code:
[10.5281/zenodo.18118110](https://doi.org/10.5281/zenodo.18118110)


Please contact [fzigunov@syr.edu](fzigunov@syr.edu) if you have any questions or trouble using this work!

