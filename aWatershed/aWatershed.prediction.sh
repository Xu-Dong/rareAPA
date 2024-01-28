#!/bin/bash

# Run using Watershed exact inference
model="Watershed_exact"  # Can take on "Watershed_exact", "Watershed_approximate", "RIVER"
number_of_dimensions="3" # Can take on any real number greater than or equal to one
training_input_file="../train04.N2pair.watershed.input.txt" # Input file txt_file_path=./"${folder_name}.txt"
# training_input_file="../pred.testdata_fortraining.txt"
prediction_input_file="../pred.data/pred.allgene.data.txt"
dirich=0.1
p_thresh=0.005
output_prefix="allgene_dirich"$dirich"_thresh"$p_thresh"_"$model
/usr/bin/Rscript predict_watershed.R --training_input $training_input_file --prediction_input $prediction_input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix --model_name $model -p $dirich --binary_pvalue_threshold $p_thresh
