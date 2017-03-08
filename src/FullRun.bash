set -e 

input_file=$1
out_dir=$2/
prefix=$3
num_iters=$4
sciclone_dir=$5
out=$out_dir$prefix

sciclone_dir=../../ancestree-dirichlet-clustering/clustering/

#### Cluster
pushd $sciclone_dir
python runSciClone_clustering.py --file_name $input_file  --output_dir ${out_dir} --prefix $prefix
popd
##
##echo "Ran Clustering"
##
num_clusters_sc=`python get_num_clusters.py ${out}.beta_params`
##
pp_matrix_file=./PPMatrix/GammaMatrix${num_clusters_sc}.txt
echo "Passed Clustering"
#
### Run Importance Sampling
#{ time python importance_sampling_distance.py ${out}.input $pp_matrix_file $num_iters ${out}.beta_params ${out} ${out}.true 2> /dev/null ; } 2> ${out}.IS.time
python importance_sampling_enum.py $input_file $pp_matrix_file $num_iters ${out}.beta_params ${out} ${out}.true
#
echo "Passed Sampling"

