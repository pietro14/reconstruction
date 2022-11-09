#! /bin/bash

RUNLOG='Runlog.csv'
RUNMIN=1
RUNMAX=99999

usage(){
    echo "The script runs background scripts:"
    echo "options:"
    
    echo "-h|--help) "
    echo "-r|--runlog) "
    echo "-m|--runmin) "
    echo "-M|--runmax) "
}

# options may be followed by one colon to indicate they have a required argument
if ! options=$(getopt -ou m:M:r:h -l help,runlog:,runmin:,runmax: -- "$@")
then
    echo "stocazzo"
# something went wrong, getopt will put out an error message for us
exit 1
fi
set -- $options
while [ $# -gt 0 ]
do
    case $1 in
	-h|--help) usage; exit 0;;
	-r|--runlog) RUNLOG=$2 ; shift ;;
	-m|--runmin) RUNMIN=$2 ; shift ;;
	-M|--runmax) RUNMAX=$2 ; shift ;;
#	--) shift; break;;
#	-*) usage; echo "$0: error - unrecognized option $1" 1>&2; usage >> /dev/stderr; exit 1;;
#	*) break;;
    esac
    shift
done

echo "Parsing runlog: = $RUNLOG and selecting runs in the range [$RUNMIN - $RUNMAX]"

echo "grep run_number $RUNLOG | sed 's|\"| |g'" > cmd.sh
echo "grep \"S000:\" $RUNLOG | sed 's|\"| |g' | awk '\$1>936 && \$1<9999{print  }' " >> cmd.sh 
source cmd.sh > parsed_runlog.csv
rm cmd.sh

echo "filtered runlog is in parsed_runlog.csv"

