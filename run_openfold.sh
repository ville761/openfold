#!/usr/bin/bash

# this script precomputes the protein alignments and then runs OpenFold

shopt -s extglob

function usage {
	echo "Usage: $0 [-o output] [-m {1-5}] [-t mtype] [-d device] [-b db_dir] [-r] input.fasta"
	echo "where the option -m sets the AlphaFold models to be used (default: model 1). It can be a range, such as 1-5."
    echo "There are 5 models, so the model numbers should in the range 1-5."
    echo ""
    echo "The option -t sets the model type ('monomer' (default) or 'pTM', multimer models are currently not"
    echo "supported)." 
    echo ""
    echo "The option -d can be used to choose the device. It can be 'cpu' (default) or 'cuda:0', 'cuda:1' etc."
    echo ""
    echo "Use the -r option if you want to produce the pickle files containing all the results, plus"
    echo "the PAE and plddt plots. The PAE statistics are produced only for the pTM models."
    echo ""
    echo "Long option names are supported. See the code."
    echo ""
	echo "If the input file is not given, it's read from stdin. Then it's best to give the prefix."
	echo "The prefix may contain a path to an output directory (e.g.. output/path/prefix)."
}

function exit_abnormal() {
	usage
	exit 1
}

die() { echo "$*" >&2; exit 2; }

function needs_arg() {
    if [ -z "$OPTARG" ]; then
        die "No arg for --$OPT option"
    fi
}

# default options
models=1 # use the first model by default
mtype="monomer" # default AlphaFold model type
modelrange="1-5"
save_outputs="false"
device="cpu"
db_dir="/data/db"
hhsearch="/usr/bin/hhsearch"

while getopts "hrd:p:m:t:-:" OPT; do
    if [ "$OPT" = "-" ]; then
        OPT="${OPTARG%%=*}"       
        OPTARG="${OPTARG#$OPT}"       
        OPTARG="${OPTARG#=}"       
    fi
    case "$OPT" in
        o | output) needs_arg; output="$OPTARG" ;;
        p | prefix) needs_arg; output="$OPTARG" ;; # for compatibility with an older version
        m | models) needs_arg; models="$OPTARG" ;;
        t | model_type) needs_arg; mtype="$OPTARG" ;;
        d | device) needs_arg; device="$OPTARG" ;;
        b | data) needs_arg; db_dir="$OPTARG" ;;
        e | hhsearch) needs_arg; hhsearch="$OPTARG" ;;
        s | save_output) save_output="true" ;;
        h | help) usage; exit 2 ;;
    esac
done
shift $((OPTIND-1))
if [ $# -gt 0 ]; then
    fasta=$1 # if the file name is not given, the sequence is read from stdin
fi

if [ -z "$output" ]; then
    echo "no output directory given"
    exit_abnormal
fi

# model type should be either "pTM" or "monomer" (case insensitive)
if [ "${mtype,,}" != "ptm" ] && [ "${mtype,,}" != "monomer" ]; then
   exit_abnormal
fi

# check if the option -m includes a range of models or just one model number
if [[ "$models" == @([$modelrange])-@([$modelrange]) ]]; then
    model0=${models%%-+([0-9])}
    model1=${models##+([0-9])-}
elif [[ "$models" == @([$modelrange]) ]]; then
    echo option 2
    model0=$models
    model1=$models
else
    exit_abnormal
fi

# check if the output directory includes a path (if given on the command line)
# if it doesn't include an absolute path, put the output directory in the current directory
if [[ "$output" =~ .*"/".* ]]; then
	  outputdir=$output
	  output=$( basename $output )
	  [[ $outputdir != /* ]] && outputdir=$PWD/$outputdir
else
	  outputdir=$PWD/$output
fi

if [ ! -d ${outputdir} ] && [ ! -h ${outputdir} ]
then
	mkdir -p ${outputdir}
fi

aligndir=${outputdir}/msa
if [ ! -d ${aligndir} ] && [ ! -h ${aligndir} ]
then
	mkdir -p ${aligndir}
fi

# check that the sequence if given either in stdin or as a file
if [ -z "$fasta" ] && [ -t 0 ] 
then
	echo "no sequence given" 
	exit_abnormal
fi

# Read the sequence from the given file or stdin and write it to a separate file.
# Also read the prefix from the file if not given on the command line.

# to fix: the following makes sense only if the fasta file includes only one sequence 
fastadir=${aligndir}/fasta
tmpfasta=${fastadir}/tmp.fasta
mkdir -p ${fastadir}
rm ${tmpfasta} 2> /dev/null
touch ${tmpfasta}
lseq=0
while IFS= read -r line || [ -n "$line" ]
do 
    if [[ $line = \>* ]]
    then
        # take the characters on the first line until the first space (OpenFold will otherwise return an error)
        line=$(cut -d ' ' -f 1 <<< "$line")
        line=$(echo $line | tr -cd '[:alnum:]>_')
	    echo $line >> ${tmpfasta}
        prefix=${line#*>}
    else
	    echo -n $line >> ${tmpfasta}
        lseq=$((lseq + $(echo -n $line | wc -c))) 
    fi
done < "${fasta:-/dev/stdin}"
#done < /dev/stdin
echo "" >> ${tmpfasta}
echo Lseq: ${lseq}

# append the prefix to the fasta file
#echo -e ">${prefix}\n$(cat ${tmpfasta})" > ${tmpfasta}

fasta=${tmpfasta}

scriptdir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

cd ${scriptdir}
source scripts/activate_conda_env.sh

echo "Precomputing protein alignments..."
python scripts/precompute_alignments_mmseqs_api.py ${fasta}  \
	    uniref30_2103_db \
	    ${aligndir}/ \
    --hhsearch_binary_path ${hhsearch} \
    --env_db colabfold_envdb_202108_db \
    --pdb70 ${db_dir}/pdb70/pdb70

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
for (( modelnumber=$model0; modelnumber<=$model1; modelnumber++ ))
do
    if [ "${mtype,,}" = "ptm" ]; then
        model="${modelnumber}_ptm"
    else
        model="$modelnumber"
    fi
    echo Calculating model: $model
    #cmd="python run_pretrained_openfold_all_ptm_with_plots.py ${tmpfasta} ${db_dir}/pdb_mmcif/mmcif_files/"
    cmd="python run_pretrained_openfold.py ${tmpfasta} ${db_dir}/pdb_mmcif/mmcif_files/"
    cmd="$cmd --use_precomputed_alignments ${aligndir}/"
    cmd="$cmd --output_dir ${outputdir}"
    cmd="$cmd --model_name model_${model}"
    cmd="$cmd --uniref90_database_path ${db_dir}/uniref90/uniref90.fasta"
    cmd="$cmd --mgnify_database_path ${db_dir}/mgnify/mgy_clusters_2018_12.fa"
    cmd="$cmd --pdb70_database_path ${db_dir}/pdb70/pdb70"
    cmd="$cmd --uniclust30_database_path ${db_dir}/uniclust30/uniclust30_2018_08/uniclust30_2018_08"
    cmd="$cmd --model_device ${device}"
    if [ $save_output = "true" ]; then
        cmd="$cmd --save_outputs true"
    fi
    results=$($cmd)
    echo python mkplots.py --model_name $model --lseq ${lseq} ${outputdir}/${prefix}_model_${model}_output_dict.pkl
    python mkplots.py --model_name $model --lseq ${lseq} ${outputdir}/${prefix}_model_${model}_output_dict.pkl
done

source scripts/deactivate_conda_env.sh
