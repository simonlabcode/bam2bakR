## Common challenges

If there are very few T-to-C mutations in the final cB.csv file (e.g., if sample-wide mutation rates in +s4U samples are < 0.003), then you may have used the incorrect value for the `strandedness` parameter in the config. One way to tell if this is the case is by looking at one of the +s4U sample counts.csv files in `results/counts/` and checking for an abundance of A-to-G mutations. If this is the case, flip the value of `strandedness` to the opposite of whatever you used.

Related to the first point, a good sanity check after running the pipeline is going into R and checking the raw mutation rates as such:

```
library(data.table)

# To unzip and read cB, also need to have R.utils package installed
cB <- fread("path/to/cB.csv.gz")

# Assess sample-wide T-to-C mutation rate in each sample
cB[,.(mutrate = sum(TC*n)/sum(nT*n), by = sample]
  # Want to see that +s4U samples has higher mutation rate than -s4U samples
```

Similarly, checking a counts.csv file for an abundance of A-to-G mutations can be done as follows:

```
library(data.table)

counts <- fread("path/to/+s4U/counts.csv.gz")

## Check if A-to-G mutation rate is higher than T-to-C mutation rate:

# A-to-G mutation rate
sum(counts$AG)/sum(counts$nA)

# T-to-C mutation rate
sum(counts$TC)/sum(counts$nT)
```