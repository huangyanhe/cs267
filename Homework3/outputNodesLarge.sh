echo "interactive session starts"
echo "=========================="
export GASNET_MAX_SEGSIZE=4G
export UPCXX_SEGMENT_MB=4000
salloc -N 1 -A mp309 -t 10:00 -q debug --qos=interactive -C haswell srun -N 1 -n 32 ./kmer_hash /global/project/projectdirs/mp309/cs267-spr2018/hw3-datasets/large.txt verbose
export GASNET_MAX_SEGSIZE=3G
export UPCXX_SEGMENT_MB=2048
salloc -N 2 -A mp309 -t 10:00 -q debug --qos=interactive -C haswell srun -N 2 -n 64 ./kmer_hash /global/project/projectdirs/mp309/cs267-spr2018/hw3-datasets/large.txt verbose
export UPCXX_SEGMENT_MB=1024
salloc -N 4 -A mp309 -t 10:00 -q debug --qos=interactive -C haswell srun -N 4 -n 128 ./kmer_hash /global/project/projectdirs/mp309/cs267-spr2018/hw3-datasets/large.txt verbose
export UPCXX_SEGMENT_MB=512
salloc -N 8 -A mp309 -t 10:00 -q debug --qos=interactive -C haswell srun -N 8 -n 256 ./kmer_hash /global/project/projectdirs/mp309/cs267-spr2018/hw3-datasets/large.txt verbose
#echo "================"
#cat test*.dat | sort > my_solution.txt 
#diff my_solution.txt /global/project/projectdirs/mp309/cs267-spr2018/hw3-datasets/test_solution.txt
