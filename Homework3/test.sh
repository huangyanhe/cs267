echo "interactive session starts"
echo "=========================="
export GASNET_MAX_SEGSIZE=2G
export UPCXX_SEGMENT_MB=128
salloc -N 1 -A mp309 -t 10:00 -q debug --qos=interactive -C haswell srun -N 1 -n 2 ./kmer_hash /global/project/projectdirs/mp309/cs267-spr2018/hw3-datasets/test.txt verbose
#echo "check the output"
#echo "================"
#cat test*.dat | sort > my_solution.txt 
#diff my_solution.txt /global/project/projectdirs/mp309/cs267-spr2018/hw3-datasets/human-chr14-synthetic_solution.txt
