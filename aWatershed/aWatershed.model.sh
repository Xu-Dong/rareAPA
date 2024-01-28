#!/bin/bash

# Run using Watershed exact inference
model="Watershed_exact"  # Can take on "Watershed_exact", "Watershed_approximate", "RIVER"
number_of_dimensions="3" # Can take on any real number greater than or equal to one
input_file="../train04.N2pair.input.txt" # Input file txt_file_path=./"${folder_name}.txt"
dirich=0.1
p_thresh=0.005
p_fraction=0.05
output_prefix="train04_dirich"$dirich"_thresh"$p_thresh"_pfrac"$p_fraction"_"$model
/usr/bin/Rscript evaluate_watershed.R --input $input_file --number_dimensions $number_of_dimensions --output_prefix $output_prefix --model_name $model -p $dirich --pvalue_threshold $p_thresh --pvalue_fraction $p_fraction --l2_prior_parameter NA
