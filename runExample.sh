# build ASMC
sh build.sh

# prepare it to run on the example array data
sh prepare.sh EXAMPLE/exampleFile.array

# analyze the array data in 10 batches
for n in {1..10}; do
	sh decode.sh EXAMPLE/exampleFile.array array $n 10
done

# merge the results of the 10 batches
sh merge.sh EXAMPLE/exampleFile.array 10

# if you want to clean up:
# sh clean.sh
