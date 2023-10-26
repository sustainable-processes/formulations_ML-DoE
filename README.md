# ML-guided MaxProQQ Designs for Liquid Formulations

We developed a (weighted-) space filling design for liquid formulations. These codes and files support our manuscript with this title. 

The ML-guided DoE is in two main codes. Firstly, the `Formulations_DoE.R` script is used to generate a MaxProQQ design and read our experimental data. Then, we use the `PhaseStability_ML.ipynb` Python notebook to train the best performing phase-stability classifier over this experimental data and use the classifier to predict feasible regions of the design space. A file outputted from this notebook (`RestrictedDesign.csv`) needs to be read into the R script which uses the `MaxProAugment()` function to then suggest the next batch of experiments to perform. 

The CSV files on this repository are required dependencies to run parts of the code. And finally, one can load an experimental set and import a previously tuned model in the phase stability Jupyer notebook to investigate the results further.

