#!/bin/bash

DIH="$1"

# convenience
#options="-num_res 1 -x_label $\chi{_1}$ -y_label $\chi{_2}$ -contour -cmap nipy_spectral_r -levels 250 -font_size 20"
options="-num_res 1 -x_label $\chi{_1}$ -y_label $\chi{_2}$ -font_size 20 -alpha 0.5"


# Run
/dors/meilerlab/data/belle6/miniforge3/envs/brownbp1/bin/python plot_rama.py -input_file $DIH -output_file $DIH.png $options
