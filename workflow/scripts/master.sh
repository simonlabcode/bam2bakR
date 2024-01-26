#!/bin/bash

# Source the paths and variables:
cpus=$1
cBout=$2
keepcols=$3
mut_tracks=$4
directory=$5
mut_pos=$6

if [ "$mut_pos" = "True" ]; then

    pos_cutoff=$7
    mutposout=$8
    mutposfilter=$9
    high_cutoff=${10}

fi


#day=$(date +"%y%m%d")

# create base count name: "TC" => "nT"
base=$(echo $mut_tracks | awk -v OFS="," '
                               {
                                    split($1,a,",")
                               }
                               END {
                                   for (base in a) {
                                       $base = "n"substr(a[base],1,1)
                                   }
                                   print $0
                               }')



keepcols=${keepcols}","${base}","${mut_tracks}



# I feel like this should work: cat <(echo Filename:$(basename "{1}" _counts.csv.gz))


# Read all _counts.csv.gz files and save them as master-DATE.csv.gz and cB-DATE.csv.gz
parallel -j $cpus --compress --plus "base={1}; base=\$(basename \${base}); cat <(echo Filename:\${base%_counts.csv.gz}) <(pigz -d -k -c -p 1 {1} | sed 's/\r$//')" ::: $directory/*_counts* \
    | awk -v OFS="," '
            $1 ~ /Filename/ {
                split($1, sample, ":")
                next
            }
            NR == 2 {
                header = $0
                print "sample", $0
                next
            }
            $0 == header {
                next
            }
            {
                print sample[2], $0
            }' \
    | awk -v colNames="$keepcols" '
            BEGIN {
                FS=","
                ncol=split(colNames,tmp)
                for (i in tmp) {
                    names[tmp[i]]
                }
            }
            NR == 1 {
                for (i=1; i<=NF; i++) {
                    if ($i in names) {
                        f[++nf] = i
                    }
                }              
            }
            {
                out = ""
                for (i=1; i<=ncol; i++) {
                    out = (out=="" ? "" : out ",") $(f[i])
                }
                print out
            } ' \
    | awk -v OFS="," -v FS="," '
            NR == 1 {
                print $0, "n"
                next
            }
            NR == 2 {
                sample = $1
            }
            sample != $1 {

                for (row in count) {
                    print row, count[row]
                }
                sample = $1
                delete count
            }
            {
                ++count[$0]
            }
            END {
                for (row in count) {
                    print row, count[row]
                }
            } '\
	| pigz -p $cpus > "$cBout"


# Explanation of the data flow:
#
#   cat (echo [filename as identiefier of data origin] + gzip [decompress file]) [use Parallels to create a single stream of alternating filename line + csv lines]
#    |
#    --> awk [use filename as a value for new column named "sample"; keep only the very first .csv header]
#         |
#         --> tee  --> gzip [compress to master file]
#              |
#              --> awk [look which column names in header match keepcols and print only those]
#                   |
#                   --> awk [calculate occurences of reads with the same values]
#                        |
#                        --> gzip [compress output to cB file]

echo "** cB file created: cB.csv.gz"

# Clean up files
rm -f 0


# Read all _cU.csv.gz files and save them as cU-DATE.csv.gz
if [ "$mut_pos" = "True" ]; then
    parallel -j $cpus --plus "cat <(echo Filename:{1%_cU.csv.gz}) <(pigz -d -k -c -p 1 {1})" ::: ./results/counts/*_cU.csv.gz \
        | awk -v OFS="," '
                $1 ~ /Filename/ {
                    split($1, sample, ":")
                    next
                }
                NR == 2 {
                    header = $0
                    print "sample", $0
                    next
                }
                $0 == header { 
                    next
                }
                {
                    print sample[2], $0
                }' \
        | awk '{ if (NR > 1) {$1 = substr($1, 18); print } else print }' \
        | pigz -p $cpus > "$mutposout"


   pigz -d -c "$mutposout" \
	| awk -F "," \
	      -v cutoff="$pos_cutoff" \
          -v upper="$high_cutoff" \
	      'NR == 1 || $(NF-1) >= cutoff && $(NF-1) <= upper { print}' | pigz -p $cpus >  "$mutposfilter"


echo "**  site-specific mutation file created: mutpos.csv.gz"

fi
