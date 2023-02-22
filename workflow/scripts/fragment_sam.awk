#! /usr/bin/awk -f 
# Sep file into fragments 
# pass variables: fragment_size [-v fragment_size="300000"]
#                 sample [-v sample="my_sample_name"]

BEGIN { i = 1 }

#$1 ~ /^@/ { NR = 0
#           next 
#         }

# Remember read name for second last read
NR == (fragment_size * i - 1 ) { x = $1 }

# Keep adding reads to file until fragent size limit is reached
NR < (fragment_size * i ) { print >> i"_"sample".sam" }

# If the last read is a pair to previous one, add it to the same file. Othewise, start new fragment file.
NR == (fragment_size * i ) { if ($1 == x) 
                                { 
                                    print >> i"_"sample".sam"
                                    i++
                                } 
                            else 
                                {           
                                    i++
                                    print >> i"_"sample".sam"
                                }
                        
                            }
