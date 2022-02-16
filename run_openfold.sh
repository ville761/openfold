#!/usr/bin/bash

# this script precomputes the protein alignments and then runs OpenFold

outputdir=$PWD

function usage {
	echo Usage: $0 [-m {1-5}] input.fasta
	echo "where the option sets the AlphaFold model to be used (default: 1)"
	exit 1
}

model=1 # use the first model by default

getopts "m:" option
if [ "$option" = "m" ]; then
	model=${OPTARG}
	if [ ${OPTIND} -eq 2 ]; then
		shift
	else
		shift; shift;
	fi
fi

if [ $# -lt 1 ]; then
	usage
fi

fasta=$1
#suffix=${fasta##*.}
#
#if [[ "$suffix" != "fasta" && "$suffix" != "seq" ]] ; then
#	usage
#fi

aligndir=$PWD/alignment_dir
scriptdir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

if [ ! -d ${aligndir} ] && [ ! -h ${aligndir} ]
then
	mkdir ${aligndir}
fi

tmpfasta=${aligndir}/tmp.fasta
rm ${tmpfasta} 2> /dev/null
touch ${tmpfasta}
while IFS= read -r line
do
	if [[ $line = \>* ]]
	then
		cut -d\| -f1 <<< $line >> ${tmpfasta}
	else
		echo $line >> ${tmpfasta}
	fi
done <"$fasta"

# remove special characters
sed -i $'s/[^[:print:]\t]//g' ${tmpfasta}

names=()
while IFS= read -r line
do
	if [[ $line = \>* ]]
	then
		names+=( $( echo "$line" | cut -c 2-) ) 
	fi
done <"${tmpfasta}"

cd ${scriptdir}
source scripts/activate_conda_env.sh

echo "Precomputing protein alignments..."
python my_precompute_alignments_mmseqs.py ${tmpfasta}  \
	    uniref30_2103_db \
	    ${aligndir}/ \
    --hhsearch_binary_path /usr/bin/hhsearch \
    --env_db colabfold_envdb_202108_db \
    --pdb70 $PWD/data/pdb70/pdb70

## fix character encoding
for name in ${names[@]}
do
	adir=$aligndir/$name
	if [[ -d "${adir}" ]]
	then
		#sed -i 's|\x0/|g' $adir/bfd.mgnify30.metaeuk30.smag30.a3m
		#sed -i 's|\x0/|g' $adir/uniref.a3m
		sed -i 's/\x0//g' $adir/bfd.mgnify30.metaeuk30.smag30.a3m
		sed -i 's/\x0//g' $adir/uniref.a3m
		#cat $adir/bfd.mgnify30.metaeuk30.smag30.a3m | tr -d '\000' > $adir/bfd.mgnify30.metaeuk30.smag30.a3m
		#cat $adir/uniref.a3m | tr -d '\000' > $adir/uniref.a3m
		#ex -s +"%s/\%x00//g" -cwq $adir/bfd.mgnify30.metaeuk30.smag30.a3m
		#ex -s +"%s/\%x00//g" -cwq $adir/uniref.a3m
	fi
done

echo "Running inference on the sequence(s) using DeepMind's pretrained parameters..."
#python run_pretrained_openfold.py \
#	${tmpfasta} \
#	/data/db/uniref90/uniref90.fasta \
#	/data/db/mgnify/mgy_clusters_2018_12.fa \
#	/data/db/pdb70/pdb70 \
#	/data/db/pdb_mmcif/mmcif_files/ \
#	/data/db/uniclust30/uniclust30_2018_08/uniclust30_2018_08 \
#	--output_dir ${outputdir} \
#	--use_precomputed_alignments ${aligndir}/
# for the new version
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
