# this script will run FastSMC on a simulated data as described in the paper (in FILES/FASTSMC_EXAMPLE/)
# parameters can be changed if desired

cd ../ASMC_BUILD_DIR/ || exit

./FastSMC_exe --inFileRoot ../FILES/FASTSMC_EXAMPLE/example \
  --outFileRoot ../c++_example/FastSMC_output_example \
  --decodingQuantFile ../FILES/FASTSMC_EXAMPLE/example.decodingQuantities.gz \
  --mode array \
  --time 50 \
  --min_m 1.5 \
  --segmentLength \
  --GERMLINE \
  --perPairPosteriorMeans \
  --perPairMAP \
  --noConditionalAgeEstimates \
  --bin

# Binary output file can be converted with the following command line
echo 'Showing first lines of the binary output...'
./convertBinary_exe ../c++_example/FastSMC_output_example.1.1.FastSMC.bibd.gz | head
