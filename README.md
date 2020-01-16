# scpred-workflow
Nextflow pipeline for [scpred-cli](https://github.com/ebi-gene-expression-group/scpred-cli) package. 
Workflow can be run for training/evaluating a model or for predicting labels on new data.   

## Evaluating the Model

The goal of this branch of the workflow is to train the prediction model and evaluate its performance. 
Method has to be set to "evaluation" in the nextflow.config file.

As an input, a training 10X dataset is given. 

The SingleCellExperiment object is then splited into its training and test subsets. Eigenvalue decomposition and feature selection are then performed with the train subset to train the model. Once trained, cell labels are predicted on the test subset. 

As an output, the predicted cell labels, the prediction probabilities, and a confusion matrix comparing predicted and actual annotations. 

## Predicting cell labels 

This branch of the workflow is for actual cell label prediction with query data. 
Method has to be set to "prediction" in the nextflow.config file.

As an input, a query dataset in the 10X format is provided along with the training dataset. 
The model is trained with the entire traing dataset (following the same eigenvalue decomposition and feature selection as in the evaluation mode) and applied to the query dataset. 
As an output, a table with the predicted cell labels of the query dataset.  

