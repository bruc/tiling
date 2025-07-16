#!/bin/sh

driver=SQLite
if [[ "$1" != "" ]]
then
    driver="$1"
fi

function execute_sql ()
{
    local sql="$1"
    local db="$2"

    if [[ "$db" == "" ]]
    then
	db=$dbname
    fi

    if [[ "$driver" == "Pg" ]]
    then
	psql -c "$1" $db
    else
	sqlite3 -separator " | " -header -cmd "$1" $db < /dev/null
    fi
}

if [[ "$driver" == "Pg" ]]
then
    .  /pg/postgresql-16.4/profile
fi

set -x

frescobi_bin=$BUILD_PREFIX/frescobi/scripts
export PATH=${frescobi_bin}:$BUILD_PREFIX/binaries/blast-2.14.0+/bin:$PATH

d_path=$(readlink -f ../data/2023-Revio-ssAAV-pAV-CMF-GFP)

declare -a fa_s
declare -a names

barcodes=rep1:rep2

templates=ssAAV

fa_s=($(readlink -f ../references/pAV-CMV-GFP_split.fasta))
      
names=(ssAAV)

n_names=${#fa_s[@]}

run=2023-Revio

if [[ "$driver" =~ ^(Pg|SQLite)$ ]]
then
    :
else
    echo The database driver can only be Pg or SQLite.
    exit 1
fi

echo At `date`: Starting $run

echo Processing $run $d_path

tiling_d="../tiling/2023-Revio-ssAAV-pAV-CMF-GFP"
if [[ ! -d $tiling_d ]]
then
    mkdir -p $tiling_d || exit 1
fi
pushd $tiling_d

# Create the database if it doesn't exist.
# Also, get the current time in the appropriate formate.

if [[ "$driver" == "Pg" ]]
then
    dbname=$run
    if psql -l | grep -s -q $dbname
    then
	dropdb $dbname
    fi
    createdb $dbname
    psql -e $dbname < $frescobi_bin/create.sql
    curtime="now()"
else
    dbname=$(readlink -f ${run}.sqlite)
    if [[ -r $dbname ]]
    then
	rm $dbname
    fi
    sqlite3 $dbname < $frescobi_bin/create_sqlite.sql
    curtime="'$(date '+%Y-%m-%d %H:%M:%S')'"
fi

# Populate the references
for ((i=0; i<n_names; i++))
do
    $BUILD_PREFIX/binaries/blast-2.14.0+/bin/makeblastdb -in ${fa_s[i]} -dbtype nucl
    execute_sql "
delete from dbs where common_db_name = '${names[i]}';
insert into dbs (common_db_name, file_name, last_update_time, data_type)
       values ('${names[i]}',
	       '${fa_s[i]}',
	       $curtime,
	       'nucleic');
"
done
# Load the sequences.
echo At `date`: Loading sequences for  $run

if [[ ! -d tmp ]]
then
    mkdir tmp || exit 1
fi
export TMPDIR=$(pwd)/tmp
err=0
ibc=1
for barcode in $(echo ${barcodes} | tr ':' ' ')
do
    fa=$(echo ${d_path}/$barcode/0-reads/*.hifi_reads.fasta)
    if [[ ! -r $fa ]]
    then
	echo Unable to find $fa
	err=1
	continue
    fi
    loadseq.pl -dbname=$dbname \
	       -driver=$driver \
	       -library=$barcode \
	       -repository=cgen \
	       -uniquify \
	       -uniquify_separator=_ \
	       -skipdupl \
	       $fa
done
if ((err))
then
    echo Exiting with error due to missing files.
    exit 1
fi

rm -rf tmp/Dbseq.*

for table in annotation_runs annotator_history hit_id_record hits hsps
do
    execute_sql "delete from $table" 
done
if [[ "$driver" == "Pg" ]]
then
    psql -c "analyze verbose" $dbname
    psql -c "drop index annotation_runs_db_seqid" $dbname
fi

echo At `date`: Annotating sequences for  $run
n_barcodes=$(echo ${barcodes} | tr ':' '\012' | wc -l)
for ((ib=1; ib<=$n_barcodes; ib++))
do
    lib=$(echo ${barcodes} | tr ':' '\012' | head -$ib | tail -1)
    db=$(echo ${templates} | tr ':' '\012' | head -$ib | tail -1)
    if [[ "$db" == "" ]]
    then
	echo Unable to find database for $lib
	exit 1
    fi
    for db in $db
    do
	annotator.pl -db=$db \
	    -maxproc=12 \
	    -sql="select seqid from raw_seqs where library='$lib'" \
	    -dbname=$dbname \
	    -driver=$driver \
	    -method=blast+blastn \
	    -block=50 \
	    -nofilter \
	    -idcutoff=0.90 \
	    -genome_size=5000 \
	    -extraopts=" -task blastn" \
	    -tmp=$(pwd)/tmp/${db}.${lib} \
	    > annotate_${db}_${lib}.log 2>&1 || exit 1
	rm -rf $(pwd)/tmp/${db}.${lib}
    done
done
if [[ "$driver" == "Pg" ]]
then
    psql -e -c 'vacuum analyze;' $dbname
fi

# Now do the tiling:
echo At `date`: Tiling for $run

scriptdir=$(readlink -f ../../..)
by_strand="-by_strand"

if [[ ! -d old ]]
then
    mkdir -p old || exit 1
fi

folder=phase1
if [[ -d $folder ]]
then
    mv --backup=numbered $folder old || exit 1
fi
mkdir $folder || exit 1
pushd $folder

for ((ib=1; ib<=$n_barcodes; ib++))
do
    lib=$(echo ${barcodes} | tr ':' '\012' | head -$ib | tail -1)
    db=$(echo ${templates} | tr ':' '\012' | head -$ib | tail -1)
    echo At `date`: Processing $lib
    if [[ "$db" == "" ]]
    then
	echo Unable to find database for $lib
	exit 1
    fi
    $scriptdir/tile.pl \
	-dbname=$dbname \
	-driver=$driver \
	-db=$db \
	-library=${lib}  \
	-gap_limit=40 \
	$extra_other \
	$extra_arg \
	$by_strand \
	-end_gap=3 \
	-details \
	2>${lib}.tile.err 
done
popd

echo At `date`: Done
      
