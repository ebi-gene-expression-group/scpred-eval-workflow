#!/usr/bin/env nextflow 

///////// training and testing model /////////

// extract training and prediction datasets

TRAINING_DATA = Channel.fromPath(params.training_10x_dir)
TRAINING_METADATA = Channel.fromPath(params.metadata_file)
process read_training_data{
  conda "${baseDir}/envs/dropletutils.yaml"

  errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
  maxRetries 10
  memory { 16.GB * task.attempt }

  input: 
    file(training_10X_data) from TRAINING_DATA
    file(training_metadata) from  TRAINING_METADATA

  output:
    file("training_sce.rds") into TRAINING_SCE

  """
  dropletutils-read-10x-counts.R\
            --samples ${training_10X_data}\
            --col-names ${params.col_names}\
            --metadata-files ${training_metadata}\
            --cell-id-column ${params.cell_id_col_name}\
            --output-object-file training_sce.rds
  """
}

TRAIN_TEST_SPLIT = Channel.create()
PROCESS_TRAIN_SCE = Channel.create()

// map to re-direct input into correct channel 
method = params.method
channels = ["evaluation":0, "prediction":1]
TRAINING_SCE.choice(TRAIN_TEST_SPLIT, PROCESS_TRAIN_SCE){channels[method]}

process eval_train_test_split{
  conda "${baseDir}/envs/scpred.yaml"

  errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
  maxRetries 10
  memory { 16.GB * task.attempt }

  input:
    file(training_sce) from TRAIN_TEST_SPLIT

  output:
    file("training_matrix.rds") into EVAL_TRAINING_MATRIX
    file("test_matrix.rds") into EVAL_TEST_MATRIX
    file("training_labels.txt") into EVAL_TRAINING_LABELS
    file("test_labels.txt") into EVAL_TEST_LABELS

  """
  scpred_train_test_split.R\
            --input-sce-object ${training_sce}\
            --normalised-counts-slot ${params.normalised_counts_slot}\
            --training-matrix training_matrix.rds\
            --test-matrix test_matrix.rds\
            --training-labels training_labels.txt\
            --test-labels test_labels.txt\
            --cell-types-column ${params.cell_types_col_name}\
  """
}

process eval_eigen_decompose{
  conda "${baseDir}/envs/scpred.yaml"

  errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
  maxRetries 10
  memory { 16.GB * task.attempt }

  input:
    file(training_matrix) from EVAL_TRAINING_MATRIX
    file(training_labels) from EVAL_TRAINING_LABELS

  output:
    file("scpred_training_object.rds") into SCPRED_TRAINING_OBJECT 

  """
  scpred_eigen_decomp.R\
            --training-matrix ${training_matrix}\
            --log-transform ${params.log_transform}\
            --training-labels ${training_labels}\
            --output-path scpred_training_object.rds 
  """
}

process eval_get_features{
  conda "${baseDir}/envs/scpred.yaml"

  errorStrategy { task.attempt < 5  ? 'retry' : 'finish' }   
  maxRetries 5
  memory { 16.GB * task.attempt }

  input:
    file(scpred_training_object) from SCPRED_TRAINING_OBJECT

  output:
    file("scpred_training_features.rds") into SCPRED_TRAINING_FEATURES

  """
  scpred_get_feature_space.R\
          --input-object ${scpred_training_object}\
          --prediction-column ${params.cell_types_col_name}\
          --output-path scpred_training_features.rds\
          --eigenvalue-plot-path ${params.eigenvalue_plot_path}
  """
}

process eval_train_model{
  conda "${baseDir}/envs/scpred.yaml"

  errorStrategy { task.attempt < 5  ? 'retry' : 'finish' }   
  maxRetries 5
  memory { 16.GB * task.attempt }

  input:
    file(scpred_training_features) from SCPRED_TRAINING_FEATURES

  output:
    file("eval_trained_model.rds") into EVAL_TRAINED_MODEL

  """
  scpred_train_model.R\
          --input-object ${scpred_training_features}\
          --model ${params.model}\
          --output-path eval_trained_model.rds\
          --train-probs-plot ${params.train_probs_plot_path}
  """
}

process eval_predict_labels{
  publishDir "${baseDir}/data/evaluation_outputs", mode: 'copy'
  conda "${baseDir}/envs/scpred.yaml"

  errorStrategy { task.attempt < 5  ? 'retry' : 'finish' }   
  maxRetries 5
  memory { 16.GB * task.attempt }

  input:
    file(eval_trained_model) from EVAL_TRAINED_MODEL
    file(test_matrix) from EVAL_TEST_MATRIX
    file(test_labels) from EVAL_TEST_LABELS

  """
  scpred_predict.R\
          --input-object ${eval_trained_model}\
          --pred-data ${test_matrix}\
          --test-labels ${test_labels}\
          --cell-types-column ${params.cell_types_col_name}\
          --output-path ${params.model_predictions_path}\
          --plot-path ${params.prediction_probs_path}\
          --confusion-table ${params.confusion_table_path}
  """
}


///////// predictions for new data /////////

// supply value to channel depending on selected mode 
PREDICTION_DATA = Channel.create()
if(method == 'prediction'){
  PREDICTION_DATA = Channel.fromPath(params.prediction_10x_dir)
}

process read_prediction_data{
  conda "${baseDir}/envs/dropletutils.yaml"

  errorStrategy { task.attempt < 5  ? 'retry' : 'finish' }   
  maxRetries 5
  memory { 16.GB * task.attempt }

  input:
    file(prediction_10X_data) from PREDICTION_DATA

  output: 
    file("prediction_sce.rds") into PREDICTION_SCE

    """
    dropletutils-read-10x-counts.R\
            --samples ${prediction_10X_data}\
            --col-names ${params.col_names}\
            --output-object-file prediction_sce.rds
    """
}

process pred_process_training_sce{
  conda "${baseDir}/envs/scpred.yaml"

  errorStrategy { task.attempt < 5  ? 'retry' : 'finish' }   
  maxRetries 5
  memory { 16.GB * task.attempt }

  input:
    file(training_sce) from PROCESS_TRAIN_SCE

  output:
    file("training_matrix.rds") into PRED_TRAINING_MATRIX
    file("training_labels.txt") into PRED_TRAINING_LABELS

  """
  scpred_preprocess_data.R\
            --input-sce-object ${training_sce}\
            --normalised-counts-slot ${params.normalised_counts_slot}\
            --output-matrix-object training_matrix.rds\
            --output-labels training_labels.txt 
  """
}

process pred_process_prediction_sce{
  conda "${baseDir}/envs/scpred.yaml"

  errorStrategy { task.attempt < 5  ? 'retry' : 'finish' }   
  maxRetries 5
  memory { 16.GB * task.attempt }

  input:
    file(prediction_sce) from PREDICTION_SCE

  output:
    file("prediction_matrix.rds") into PREDICTION_MAT

  """
  scpred_preprocess_data.R\
            --input-sce-object ${prediction_sce}\
            --normalised-counts-slot ${params.normalised_counts_slot}\
            --output-matrix-object prediction_matrix.rds\
  """
}

process pred_eigen_decompose{
  conda "${baseDir}/envs/scpred.yaml"

  errorStrategy { task.attempt < 5  ? 'retry' : 'finish' }   
  maxRetries 5
  memory { 16.GB * task.attempt }

  input:
    file(training_matrix) from PRED_TRAINING_MATRIX
    file(training_labels) from PRED_TRAINING_LABELS

  output:
    file("scpred_training_object.rds") into PRED_TRAINING_OBJECT

  """
  scpred_eigen_decomp.R\
            --training-matrix ${training_matrix}\
            --log-transform ${params.log_transform}\
            --training-labels ${training_labels}\
            --output-path scpred_training_object.rds 
  """
}

process pred_get_features{
  conda "${baseDir}/envs/scpred.yaml"

  errorStrategy { task.attempt < 5  ? 'retry' : 'finish' }   
  maxRetries 5
  memory { 16.GB * task.attempt }

  input:
    file(scpred_training_object) from PRED_TRAINING_OBJECT

  output:
    file("scpred_training_features.rds") into PRED_TRAINING_FEATURES

  """
  scpred_get_feature_space.R\
          --input-object ${scpred_training_object}\
          --prediction-column ${params.cell_types_col_name}\
          --output-path scpred_training_features.rds\
          --eigenvalue-plot-path ${params.eigenvalue_plot_path}
  """
}

process pred_train_model{
  publishDir "${baseDir}/data/prediction_outputs", mode: 'copy'
  conda "${baseDir}/envs/scpred.yaml"

  errorStrategy { task.attempt < 5  ? 'retry' : 'finish' }   
  maxRetries 5
  memory { 16.GB * task.attempt }

  input:
    file(scpred_training_features) from PRED_TRAINING_FEATURES

  output:
    file("pred_trained_model.rds") into PRED_TRAINED_MODEL
    file("${params.train_probs_plot_path}") into PRED_TRAIN_PROBS


  """
  scpred_train_model.R\
          --input-object ${scpred_training_features}\
          --model ${params.model}\
          --output-path pred_trained_model.rds\
          --train-probs-plot ${params.train_probs_plot_path}
  """
}

process pred_predict_labels{
  publishDir "${baseDir}/data/prediction_outputs", mode: 'copy'
  conda "${baseDir}/envs/scpred.yaml"

  errorStrategy { task.attempt < 5  ? 'retry' : 'finish' }   
  maxRetries 5
  memory { 16.GB * task.attempt }

  input:
    file(pred_trained_model) from PRED_TRAINED_MODEL
    file(prediction_matrix) from PREDICTION_MAT

  output:
    file("${params.model_predictions_path}") into PRED_MODEL_PREDICTIONS

  """
  scpred_predict.R\
          --input-object ${pred_trained_model}\
          --pred-data ${prediction_matrix}\
          --output-path ${params.model_predictions_path}\
  """
}

process get_pred_output{
  publishDir "${params.results_dir}", mode: 'copy' // send the table into 'outer' workflow
  conda "${baseDir}/envs/scpred.yaml"

  errorStrategy { task.attempt < 5  ? 'retry' : 'finish' }   
  maxRetries 5
  memory { 16.GB * task.attempt }

  input:
    file(model_predicitons) from PRED_MODEL_PREDICTIONS
  output:
    file("scpred_output.txt") into FINAL_TABLE

  """
  scpred_get_std_output.R\
          --predictions-file ${model_predicitons}\
          --output-table scpred_output.txt
  """
}
