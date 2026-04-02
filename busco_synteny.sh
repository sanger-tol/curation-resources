#!/bin/bash

# Scaffold IDs need to be in final tol format. 
# Query and ref short IDs need to be two letter codes eg. "Ac". 
# Query and ref fasta files need to be in working dir.
# Must be run in curation_v2 for plotting to complete (python dependencies - matplotlib, seaborn)

myReference=''
myQuery=''
myDirPath=`pwd`
buscoLineage=''
ref_short_id='Rf'
query_short_id='Qu'
filesKeep=''

usage() { 
cat << EOF
Usage: $0 [REQUIRED -r reference_species.fasta, -q query_species.fasta -l busco_lineage] [OPTIONAL -u Query_ID, -f Ref_ID,-p wd_path, -k]

Note that both reference and query fasta headers must have chromsome names as 'SUPER_n'

REQUIRED:
  -r [file] Reference species Fasta file
  -q [file] Query species fasta file
  -l [text] Name of BUSCO lineage to test against
OPTIONAL:
  -f [text] reF short ID - two letter abbreviation for reference species (default = $ref_short_id)
  -u [text] qUery short ID - two letter abbreviation for query species (default = $query_short_id)
  -p [dir] Path to directory containing files (default = $myDirPath)
  -k  add this flag to keep BUSCO output files - all BUSCO files other than the full_table.tsv will be removed by default
EOF
	exit 1
}

while getopts "r:q:l:f:u:p:kh" opt
   do
     case $opt in
        r ) myReference=$OPTARG
	if [ ! -f $myReference ] ; then echo "File $myReference does not exist"; usage; fi
	;;
        q ) myQuery=$OPTARG
	if [ ! -f $myQuery ] ; then echo "File $myQuery does not exist"; usage; fi
	;;  
        l ) buscoLineage=$OPTARG
	if [ ! -l ]; then echo "No BUSCO lineage provided"; usage; fi
	;;
	f ) ref_short_id=$OPTARG;;
	p ) myDirPath=$OPTARG;;
	u ) query_short_id=$OPTARG;;
	k ) filesKeep='-k';;
	h | *) usage;;
     esac
done

# Run BUSCO  on ref assembly.
busco -i ${myReference} -o busco5_mini${ref_short_id} -m genome -l $buscoLineage -c 32

# Run BUSCO on query assembly
busco -i ${myQuery} -o busco5_mini${query_short_id} -m genome -l $buscoLineage -c 32

# Run circos plotting 
wait 
python busco_synteny_format_and_plot.py -r $myReference -q $myQuery -l $buscoLineage -ri $ref_short_id -qi $query_short_id -p $myDirPath $filesKeep
