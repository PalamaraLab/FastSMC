module load gcc/6.2.0
module load python/3.6.0

python3 get_csfs.py -D $1 -d $2 -n $3 -o $4 > /dev/null

