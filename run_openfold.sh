#!/usr/bin/bash

# this script precomputes the protein alignments and then runs OpenFold

function usage {
	echo Usage: $0 [-p prefix] [-m {1-5}] input.fasta
	echo "where the option -m sets the AlphaFold model to be used (default: 1)."
	echo "If the input file is not given, it's read from the stdin. Then it's best to give the prefix."
	echo "The prefix may contain a path to an output directory (e.g.. output/path/prefix)."
	exit 1
}

function exit_abnormal() {
	usage
	exit 1
}

model=1 # use the first model by default
defaultprefix="prefix"
prefix=$defaultprefix

while [ $# -gt 0 ] && [ "$1" != "--" ]; do
	while getopts ":p:m:" option; do
		case "${option}" in
			p)
				prefix=${OPTARG}
				;;
			m)
				model=${OPTARG}
				;;
			:)
				echo "Error: -${OPTARG} requires an argument."
				exit_abnormal
				;;
			*)
				exit_abnormal
				;;
		esac
	done
	shift $((OPTIND-1))
	while [ $# -gt 0 ] && ! [[ "$1" =~ ^- ]]; do
		fasta=$1
		shift
	done
done

if [ "$1" == "--" ]; then
  shift
  fasta=$1
fi

scriptdir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

aligndir=$PWD/alignment_dir
if [ ! -d ${aligndir} ] && [ ! -h ${aligndir} ]
then
	mkdir -p ${aligndir}
fi

# read the sequence from the given file or stdin
# also read the prefix from the file if not given on the command line

tmpfasta=${aligndir}/tmp.fasta
rm ${tmpfasta} 2> /dev/null
touch ${tmpfasta}
while IFS= read -r line
do 
	if [[ $line = \>* ]]
	then
		if [ "$prefix" = "$defaultprefix" ]; then
			prefix="${line#*>}"
			prefix="${prefix%%[^a-zA-Z0-9_]*}"
		fi
	else
		echo -n $line >> ${tmpfasta}
	fi
done < "${fasta:-/dev/stdin}"
echo "" >> ${tmpfasta}

# remove special characters
sed -i $'s/[^[:print:]\t]//g' ${tmpfasta}

# check if the prefix includes a path (if given on the command line)
# if it doesn't include an absolute path, put the output directory in the current directory
if [[ "$prefix" =~ .*"/".* ]]; then
	  outputdir=$prefix
	  prefix=$( basename $prefix )
	  [[ $outputdir != /* ]] && outputdir=$PWD/$outputdir
else
	  outputdir=$PWD/$prefix
fi

# append the prefix to the fasta file
echo -e ">${prefix}\n$(cat ${tmpfasta})" > ${tmpfasta}

if [ ! -d ${outputdir} ] && [ ! -h ${outputdir} ]
then
	mkdir -p ${outputdir}
fi

cd ${scriptdir}
source scripts/activate_conda_env.sh

echo "Precomputing protein alignments..."
python my_precompute_alignments_mmseqs.py ${tmpfasta}  \
	    uniref30_2103_db \
	    ${aligndir}/ \
    --hhsearch_binary_path /usr/bin/hhsearch \
    --env_db colabfold_envdb_202108_db \
    --pdb70 $PWD/data/pdb70/pdb70

# fix character encoding
adir=$aligndir/$prefix
if [ -d "${adir}" ]
then
	sed -i 's/\x0//g' $adir/bfd.mgnify30.metaeuk30.smag30.a3m
	sed -i 's/\x0//g' $adir/uniref.a3m
	#cat $adir/bfd.mgnify30.metaeuk30.smag30.a3m | tr -d '\000' > $adir/bfd.mgnify30.metaeuk30.smag30.a3m
	#cat $adir/uniref.a3m | tr -d '\000' > $adir/uniref.a3m
	#ex -s +"%s/\%x00//g" -cwq $adir/bfd.mgnify30.metaeuk30.smag30.a3m
	#ex -s +"%s/\%x00//g" -cwq $adir/uniref.a3m
fi

echo "Running inference on the sequence(s) using DeepMind's pretrained parameters..."
python run_pretrained_openfold.py \
	${tmpfasta} \
	/data/db/pdb_mmcif/mmcif_files/ \
	--uniref90_database_path /data/db/uniref90/uniref90.fasta \
	--mgnify_database_path /data/db/mgnify/mgy_clusters_2018_12.fa \
	--pdb70_database_path /data/db/pdb70/pdb70 \
	--uniclust30_database_path /data/db/uniclust30/uniclust30_2018_08/uniclust30_2018_08 \
	--output_dir ${outputdir} \
	--model_name model_${model} \
	--use_precomputed_alignments ${aligndir}/

source scripts/deactivate_conda_env.sh

# zip the output
cd $outputdir/..
tar zcvf $prefix.tar.gz $prefix
