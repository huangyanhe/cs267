echo "interactive session starts"
echo "=========================="
salloc -N 1 -A mp309 -t 10:00 -q debug --qos=interactive -C haswell srun -N 1 -n 1 ./kmer_hash /global/project/projectdirs/mp309/cs267-spr2018/hw3-datasets/test.txt verbose
#echo "check the output"
#echo "================"
#cat test*.dat | sort > my_solution.txt 
#diff my_solution.txt /global/project/projectdirs/mp309/cs267-spr2018/hw3-datasets/test_solution.txt
