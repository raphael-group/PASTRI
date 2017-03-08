set -e 

input_file=$1
clustering_file=$2
out_dir=$3/
prefix=$4
num_iters=$5
out=$out_dir$prefix

num_clusters_sc=`python get_num_clusters.py $clustering_file`
pp_matrix_file=./PPMatrix/GammaMatrix${num_clusters_sc}.txt

### Run Importance Sampling
python importance_sampling_enum.py $input_file $pp_matrix_file $num_iters ${out}.beta_params ${out} ${out}.true

