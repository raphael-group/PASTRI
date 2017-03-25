code_dir='../src'
out_dir=`pwd`
pushd $code_dir
    bash RunPASTRI.bash ${out_dir}/Example.input ${out_dir}/Example.beta_params $out_dir Example 5000
popd

