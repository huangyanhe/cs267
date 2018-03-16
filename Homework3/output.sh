echo "interactive session starts"
echo "=========================="
export GASNET_MAX_SEGSIZE=4G
export UPCXX_SEGMENT_MB=4000
salloc -N 1 -A mp309 -t 10:00 -q debug --qos=interactive -C haswell srun -N 1 -n 1 ./kmer_hash /global/project/projectdirs/mp309/cs267-spr2018/hw3-datasets/human-chr14-synthetic.txt output single_node.txt
export GASNET_MAX_SEGSIZE=3G
export UPCXX_SEGMENT_MB=2048
salloc -N 1 -A mp309 -t 10:00 -q debug --qos=interactive -C haswell srun -N 1 -n 2 ./kmer_hash /global/project/projectdirs/mp309/cs267-spr2018/hw3-datasets/human-chr14-synthetic.txt output single_node.txt
export UPCXX_SEGMENT_MB=1024
salloc -N 1 -A mp309 -t 10:00 -q debug --qos=interactive -C haswell srun -N 1 -n 4 ./kmer_hash /global/project/projectdirs/mp309/cs267-spr2018/hw3-datasets/human-chr14-synthetic.txt output single_node.txt
export UPCXX_SEGMENT_MB=512
salloc -N 1 -A mp309 -t 10:00 -q debug --qos=interactive -C haswell srun -N 1 -n 8 ./kmer_hash /global/project/projectdirs/mp309/cs267-spr2018/hw3-datasets/human-chr14-synthetic.txt output single_node.txt
export UPCXX_SEGMENT_MB=256
salloc -N 1 -A mp309 -t 10:00 -q debug --qos=interactive -C haswell srun -N 1 -n 16 ./kmer_hash /global/project/projectdirs/mp309/cs267-spr2018/hw3-datasets/human-chr14-synthetic.txt output single_node.txt
export UPCXX_SEGMENT_MB=128
salloc -N 1 -A mp309 -t 10:00 -q debug --qos=interactive -C haswell srun -N 1 -n 32 ./kmer_hash /global/project/projectdirs/mp309/cs267-spr2018/hw3-datasets/human-chr14-synthetic.txt output single_node.txt
#echo "check the output"
#echo "================"
#cat test*.dat | sort > my_solution.txt 
#diff my_solution.txt /global/project/projectdirs/mp309/cs267-spr2018/hw3-datasets/test_solution.txt
