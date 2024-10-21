This repository contains four scripts which use the ProDy library to analyze molecular dynamics of 7FAY, the SARS-CoV2 main protease.

gnm-fluctuations.py uses ProDy to calculate the squared C-alpha fluctuations using a Gaussian Network Model (GNM). Requires system input of protein for parsePDB and prints the maximum squared fluctuation calculated amongst protein residues.

anm-fluctuations.py functuions very similarly to gnm-fluctuations.py, but utilizes the Anisotropic Network Model (ANM).

fluctuation-diff.py calculates the absolute difference in squared fluctuations for each C-alpha between GNM and ANM analyses and prints the max difference. Then, it normalizes squared fluctuations such that the maximum squared flux is set to 1.0, and prints the max absolute difference with the normalized values. This script also generates two plots which compare the GNM and ANM analyses using squared fluctuations and normalized squared fluctuations:

![sq_fluct](https://github.com/user-attachments/assets/d2d24267-e053-4864-b7f2-7cb0d7cfd977)
![norm_sq_fluct](https://github.com/user-attachments/assets/f9994e75-329e-4085-8b91-b5d64e91198c)
