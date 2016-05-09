#!/bin/bash

workdir="/home/hosts/disk17/yael/temp-sahar-log-dup-extraction"
src_dir="/home/hosts/disk17/yael/temp-sahar-log"

function print_usage() 
{
	echo "#################################################"
	echo "# main.sh:	generate histogram files from fastq"
	echo "# usage:		main.sh workdir"
	echo "#"
	echo "# workdir - the directory which will the program use for temporary files or outputs"
	echo "#################################################"
}

# Check for parameters
if [ $# -eq 0 ]
	then
		print_usage
		exit 1
fi

# validate that input dir exists
if [ ! -d "$1" ]; then
	echo "no such directory $1"
	print_usage
	exit
fi

workdir=$1

echo "-------------------------------------------------"
echo "           generating histogram files"
echo "-------------------------------------------------"
echo "workdir: $1"
echo "-------------------------------------------------"
echo "Processing Files:"

# Go over each file and generate it's histogram file
for file in `ls -1 *cutadapt* | grep _1.fastq`
do	
	echo "-------------------------------------------------"
	second=`echo $file | sed 's/_1.fast/_2.fast/g'`
	core_name=`echo $file | gawk -F '_cutadapt' '{print $1"_cutadapt"}'`
	short=`echo $file | gawk -F '_cutadapt' '{print $1}'`
	
	echo "params:"
	echo "file: $file" 
	echo "second: $second"
	echo "short: $short"
	echo "core_name: $core_name"

	echo ""	
	echo "running perl script..."

	# cluster path
	#perl /cs/icore/yaelal/progs/exclude-dups-from-fastq-files.pl $file $second $core_name $workdir >$workdir/$short-dup-stat.info

	# local path
	perl /home/users/yael/sahar-clash/progs/exclude-dups-from-fastq-files.pl $file $second $core_name $workdir > $workdir/$short-dup-stat.info
	
	original1=`wc -l $file| gawk '{print $1/4}'`
	new=`paste $file $second | gawk '{if (NR%4 ==2) print $0}' | sort -u |wc -l| gawk '{print $1}'`

	echo ""
	echo "results:"
	echo $original1 $new $file

	echo ""
	echo "generating plot:"

	# local path
	python /home/users/amirbar/lab/duplicates/plot_graph.py -i $workdir/$short-dup-stat.info -e 2 -o $workdir/$short-dup-stat.png
done

echo "-------------------------------------------------"
echo "done!"
