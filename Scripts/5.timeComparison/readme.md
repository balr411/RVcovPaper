# Description

In order to reproduce the time comparison for nref = 1000, use the script timeComparisonSysTime.R.
To reproduce the memory comparison for nref = 1000, use the script memoryComparisonCaller.R.

In order to reproduce the time and memory comparisons for nref = 9999, use the scripts 
timeComparisonSysTimeNRef9999.R and memoryComparisonCallerNRef9999.R respectively.

The scripts memoryComparisonReader.R and memoryComparisonReaderNRef9999.R can be 
used to read in the results from the memory comparisons.

Note that the pipeline in Scripts/3.SKAT_weightedBurden/ must be completed for the 
genes on chromosome 2 and nref = 1000 and nref = 9999 for this to work, as we use
the summary statistics and VCFs from that analysis here. 
