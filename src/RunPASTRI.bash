set -e 

input_file=$1
clustering_file=$2
out_dir=$3/
prefix=$4
num_iters=$5
out=$out_dir$prefix

num_clusters_sc=`python get_num_clusters.py $clustering_file`
pp_matrix_file=./PPMatrix/Matrix${num_clusters_sc}.txt

### Run Importance Sampling
python -O multi_importance_sampling.py $input_file $pp_matrix_file $num_iters $clustering_file ${out} 

