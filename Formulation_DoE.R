### Design of Experiments for Liquid Formulations ### 

## This code presents a model-weighted space filling DoE methodology to prepare shampoo formulations.  ## 
## The model-component is related to a phase stability classifier to drive the formulations towards regions of stability. ## 
## The machine learning component is done in Python in a separate jupyter notebook "PhaseStability_ML.ipynb"
## The space-filling component utilises the MaxProQQ package written by Joseph, Gul & Ba, 2020. ## 

# Clearing workspace and setting working directory ---- 

rm(list=ls()) # This clears the workspace, ensuring that any variables saved in the environment do not effect our current code. 
setwd("/Users/ac2349/R/DoE") # Setting working directory, change as necessary.

# Installing/ loading libraries required to run this code ---- 

#install.packages('MaxPro') # one-time installation 
#install.packages('caret')
#install.packages('purrr')

library(caret)
library(MaxPro)
library(purrr)

# Function definitions for programme ---- 

# The MaxPro DoE package generates variables on the interval [0,1]
# Need to be able to inter convert between these variables and concentrations

# convert design variables [0,1] to concentrations 
conc_converter <- function(x,LB,UB) {
    conc <- x*(UB - LB) + LB
    return(conc)
}

# convert concentrations to [0,1]
invert_conc <- function(conc,LB,UB) {
    x <- (conc - LB)/(UB - LB)
    return(x)
}

# Pre-defined concentration bounds for formulations ingredients 
S_LB <- 8.0
S_UB <- 13.0

P_LB <- 1.0
P_UB <- 3.0

T_LB <- 1.0
T_UB <- 5.0

# Sample density = 1/sum(weight_frac_ingredient/density_ingredient) 
# Below function takes a dataframe with the ingredient mass % in a particular sample and returns its density

calc_density <- function(df,species_dict=ingredients_info) {
    active_ingredients <- names(df)[df != 0]
    active_conc <- df[, active_ingredients]
    
    vol_frac <- 0.0 # initalise variable 
    
    for (ingredient in colnames(active_conc)) { # iterate over each active ingredient 
        ingredient_density <- as.numeric(species_dict[species_dict['name'] == ingredient][2])
        vol_frac <-  vol_frac + (active_conc[,ingredient]/100)/ingredient_density 
    } # equivalent to summing over all ingredients 
    return(1.0/vol_frac) # returns the sample density 
}

ingredients_info <- read.csv('SpeciesDictionary.csv') # csv file contains densities for ingredients.

# Takes a candidate set of experiments generated via the MaxPro package and converts them to a format consistent with the csv file sent to the OT-2 robot.
Cand_to_Dataset <- function(Cand_matrix, surf_pairs, OT_example_csv, n_samples, polymer, thickener, calc_density=FALSE) {
    
    # First input is a matrix with 5 variables - the candidate experiments. 
    # First four variables correspond to concentrations: C_S1, C_S2, C_P, C_T
    # Fifth variable is a discrete factor representing the choice of surfactant pair. 
    
    # Create a data frame placeholder to populate with the candidate experiments 
    DoE_cols <- colnames(OT_example_csv) # the column names should match the csv file sent to the Opentrons robot.
    Cand_DoE <- data.frame(matrix(nrow=n_samples, ncol=length(DoE_cols))) # create an empty dataframe
    colnames(Cand_DoE) <- DoE_cols # rename the data frame columns to match with the target csv 
    
    # Populate the ID and Sample rows - currently the same, but created 2 separate rows in case there are repeat experiments performed.
    Cand_DoE[,"ID"] <- 1:n_samples
    Cand_DoE[,"Sample"] <- 1:n_samples
    
    Cand_df <- data.frame(Cand_matrix) # convert candidate matrix to data frame 
    
    # Populate all remaining NA values as 0. Non-active ingredient concentrations will = 0.
    Cand_DoE[is.na(Cand_DoE)] = 0 
    
    # The surfactant pairs have been encoded into 66 nominal factors - need to convert between these and the selected surfactants.
    S_factor = Cand_matrix[,5]  # read the factor for each sample in the candidate set - column 5, the last column corresponds this categorical factor. Columns 1-4 are concentrations.
    
    # read off S1, S2
    S1 = surf_pairs[S_factor,2] 
    S2 = surf_pairs[S_factor,3]
    
    surf_names <- colnames(OT_example_csv[,3:14])
    
    for (s in surf_names) {
        S1_idx <- which(S1 %in% s)
        S2_idx <- which(S2 %in% s)
        Cand_DoE[S1_idx,3:14][s] <- conc_converter(Cand_df[S1_idx,1], S_LB, S_UB)
        Cand_DoE[S2_idx,3:14][s] <- conc_converter(Cand_df[S2_idx,2], S_LB, S_UB)
    }
    
    Cand_DoE[polymer] <- conc_converter(Cand_df[,3], P_LB, P_UB)
    Cand_DoE[thickener] <- conc_converter(Cand_df[,4], T_LB, T_UB)
    
    Cand_DoE[,"Water"] <- 100 - rowSums(Cand_DoE[,3:20]) # Mass % must sum to 100. The remaining content in the formulations is water.
    
    
    if (calc_density == FALSE) {
        
        Cand_DoE[,"Sample Density"] <- 1.00 # assuming the sample density = 1 (this is a very good approximation & avoids doing the row-wise operation, over possibly 100000s candidate samples)
    
    } else  {
        
        # Extract out the ingredient names from the candidate dataframe.
        ingredient_list <- colnames(Cand_DoE[,3:21]) # indices 3:21 correspond to the 19 ingredients, including water.
        
        sample_density <- c() # Generate an empty vector to store the density values of the sample.
        
        # Compute the sample density using the pre-defined density calculation function above.
        
        for (idx in 1:n_samples) { # iterate through the samples 1-by-1
            x <- Cand_DoE[idx, ingredient_list]
            sample_density <- append(sample_density, calc_density(x)) # append the density to the vector
        }
        
        Cand_DoE[,"Sample Density"] <- sample_density
        
    }

    Cand_DoE[,3:22] <- round(Cand_DoE[,3:22], digits=2) 
    
    return(Cand_DoE)
    
}

# Candidate DoE points ---- 

set.seed(42) # fix the random seed so to ensure reproducibility between code runs.

N <- 360000 # number of candidate points, equivalent to all combinations if the conc. space was discretised into 0.5 w/w% intervals.
C <- CandPoints(N, p_cont=4, l_nom=66) # 4 continuous variables (S1, S2, P, T concentrations) & 66 possible surfactant pairs from 12 surf ingredients.
Cand_df <- data.frame(C) # convert candidate matrix to data frame 

# Write out the candidate design to a .csv file to be read by Pandas - need this to filter into the restricted candidate set based on phase stability predictions 
write.csv(Cand_df, 'CandDesign.csv')  # Note design here refers to the DoE design with variables in the interval [0,1]

OT_example_csv <- read.csv('OT-2_DoE_Example.csv', check.names = FALSE)
surf_pairs <- read.csv('SurfactantPairs.csv') # csv file to read which factor corresponds to which surfactant pair.

# ENSURE CORRECT POLYMER, THICKENER SELECTED!! 

polymer = "Dehyquart CC7 Benz"
thickener = "Arlypon TT"

Cand_DoE <- Cand_to_Dataset(C, surf_pairs, OT_example_csv, N, polymer, thickener)  # Converting the MaxPro candidate points to the format of the formulations dataset
#View(Cand_DoE)

# Read the formulations data and features ---- 
my_data <- read.csv("PhDFormulationsDataset_2023.csv", check.names=FALSE) # read the master dataset containing the samples according to their mass composition & their properties.
#my_data <- rbind(my_data[523:582, ], my_data[619:666, ]) # select only the columns with data for the current (P, T) combination.
my_data <- my_data[751:822, ]

surf_conc_expt <- my_data[2:13] # extract out the data for the 12 surfactants (Texapon SB 3 KC ... Dehyquart A-CA)
surf_conc_expt <- surf_conc_expt[complete.cases(surf_conc_expt), ] # drop any fully empty rows (sometimes erroneously saved into the dataset csv). 

# Extract the surfactant concentrations from the candidate DoE 
surf_conc_cand <- Cand_DoE[,3:14] # indices 3:14 correspond to Texapon SB 3 KC ... Dehyquart A-CA 

# Extract the concentration for the conditioning polymer (P) and thickener (T) - ensure to change ingredient name according to the choice of P,T used in the experiment.
poly_thick_expt <- data.frame(cbind(my_data[polymer], my_data[thickener])) # ARE THE CHOSEN P, T CORRECT?
poly_thick_expt <- poly_thick_expt[complete.cases(poly_thick_expt), ]
colnames(poly_thick_expt) <- c("P", "T") # rename the columns in the P,T dataframe.

# extract out the P, T concentrations for the candidate set 
poly_thick_cand <- data.frame(cbind(Cand_DoE[][polymer], Cand_DoE[][thickener]))
colnames(poly_thick_cand) <- c("P", "T") # rename the columns in the P,T dataframe.

surf_features <- read.csv("Surf_FGs.csv", row.names=1, check.names = FALSE) # read the data file containing the scaled FG count for each surfactant and its average chain length.
# check.names = FALSE so it doesn't change how the ingredient names are written - important for indexing.

# Pre-processing: calculate and scale the sample features ---- 

# Calculate the sample features via matrix multiplication of the surfactant conc. and features + append the P, T concentrations

expt_sample_features <- as.matrix(surf_conc_expt) %*% as.matrix(surf_features) # matrix dot product between the surfactant concentrations and their features results in a (#n_samples, #n_features) matrix.
expt_sample_features <- cbind(expt_sample_features, poly_thick_expt) # append the P, T concentrations onto the features

cand_sample_features <- as.matrix(surf_conc_cand) %*% as.matrix(surf_features) # matrix dot product between the surfactant concentrations and their features results in a (#n_samples, #n_features) matrix.
cand_sample_features <- cbind(cand_sample_features, poly_thick_cand) # append the P, T concentrations onto the features

# Feature scaling 

preProc <- preProcess(expt_sample_features, method = c("range")) # use min-max scaling in the caret package to pre-process the data.
X_expt <- predict(preProc, expt_sample_features) 

X_cand <- predict(preProc, cand_sample_features) # ensured the candidate set has been pre-processed in the same way. 

# Extract the response variable (phase stability) from the experimental dataset

y <- as.factor(my_data$Stability_Test[complete.cases(my_data$Stability_Test)]) # read the output (stability) from the formulations dataset and it must be coded as a factor in R.

X_expt <- round(X_expt, digits=2) 

featurePlot(X_expt, y, plot="strip")  # Plot the features grouped by stability response for feature exploration/visualisation

df <- cbind(X_expt, y) # add the output to the input features, X, to generate the overall dataset, df.
df <- as.data.frame(df) # this will be used for machine learning.

# Machine Learning part in Python, write out csv files and read back in csv files ---- 

write.csv(df, "ExperimentalSet.csv") # experimental data with the surfactant features, polymer/thickener concentrations & phase stability responses 

write.csv(X_cand, "CandidateSet.csv") # Note set refers to the candidate set of sample experiments in the concentrations space. 

# A variety of machine learning models are trained and tuned in Python on the stability classification problem 
# The most promising model at any given iteration (of the formulations workflow, i.e. week-on-week) is used to predict the phase stability of the candidate points
# Depending on the predicted likelihood of preparing stable formulations, the candidate design is filtered into a restricted set of possible experiments. This is what MaxProAugment chooses from. 

Restricted_Cand_Design <- read.csv('RestrictedCandDesign.csv')
Restricted_Cand_Design <- Restricted_Cand_Design[,2:6]

# Convert the experimental data into the candidate set format ---- 

Sample_ID <- c(1:nrow(surf_conc_expt))
df_surf_expt <- cbind(Sample_ID, surf_conc_expt)
expt_surf_pairs <- df_surf_expt %>% split(.$Sample_ID) %>% map(~names(.x)[!!.x][-1])

# Create an empty dataframe with these column names.
Exist_Design <- data.frame(matrix(nrow=nrow(df_surf_expt), ncol=length(Cand_df)))
colnames(Exist_Design) = colnames(Cand_df)

for (idx in 1:nrow(surf_conc_expt)) {
    
    S1 = unlist(expt_surf_pairs[idx])[1]
    S2 = unlist(expt_surf_pairs[idx])[2]
    
    Exist_Design[idx,1] <- invert_conc(surf_conc_expt[idx,][S1], S_LB, S_UB)
    Exist_Design[idx,2] <- invert_conc(surf_conc_expt[idx,][S2], S_LB, S_UB)
    
    Exist_Design[idx,3] <- invert_conc(poly_thick_expt[idx, ]["P"], P_LB, P_UB)
    Exist_Design[idx,4] <- invert_conc(poly_thick_expt[idx, ]["T"], T_LB, T_UB)
    
    Exist_Design[idx,5] <- which(surf_pairs[,2] == S1 & surf_pairs[,3] == S2)
}

# Write out Existing Design 

write.csv(Exist_Design, "P3_T1_Design.csv")

# Generate new sample points ---- 

nNew <- 36

#New_Design <- MaxProAugment(as.matrix(Exist_Design), Restricted_Cand_Design, nNew, p_nom=1, l_nom=66)$Design # run this line for ML-guided DoE
New_Design <- MaxProAugment(as.matrix(Exist_Design), C, nNew, p_nom=1, l_nom=66)$Design  # run this line for purely space filling DoE 

#new_samples <- New_Design[(nrow(New_Design) - nNew + 1):nrow(New_Design), ]
#new_samples 

Expt_DoE <- Cand_to_Dataset(New_Design, surf_pairs, OT_example_csv, nrow(New_Design), polymer, thickener, calc_density = TRUE)  # Converting the MaxPro candidate points to the format of the formulations dataset
View(Expt_DoE)

path_out = '/Users/ac2349/R/DoE/SuggestedExperiments/' # change this path to the folder where you want to save the suggested experiments csv file.
fileName = paste(path_out, 'P3-T1_OT_DoE_0410-09-2023.csv', sep = '')
write.csv(Expt_DoE, fileName, row.names = FALSE)


# Initialise dataset for new polymer, thickener combination ----

set.seed(42) # fix the random seed so to ensure reproducibility between code runs.

N_P3_T1_init <- 24
P3_T1_init_design <- CandPoints(N_P3_T1_init, p_cont=4, l_nom=66) # 4 continuous variables (S1, S2, P, T concentrations) & 66 possible surfactant pairs from 12 surf ingredients.

# ENSURE CORRECT POLYMER, THICKENER SELECTED!!

polymer = "Dehyquart CC7 Benz"
thickener = "Arlypon TT"

P3_T1_init <- Cand_to_Dataset(P3_T1_init_design, surf_pairs, OT_example_csv, N_P3_T1_init, polymer, thickener, calc_density = TRUE)
View(P3_T1_init)

path_out = '/Users/ac2349/R/DoE/SuggestedExperiments/'
fileName = paste(path_out, 'MasterDataset_OT_DoE_P3_T1_init.csv', sep = '')
write.csv(P3_T1_init, fileName, row.names = FALSE)