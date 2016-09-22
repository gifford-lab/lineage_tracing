## Processing procedure
#### Split the reads into 8 groups by index barcodes
Assign each read to the group of which the index barcode is uniquely closest to the read's index barcode. If the closest group is not unique or the closest hamming distance is greater than 2, the read was considered ambiguous and ignored.

```
$REPOHOME/splitreads/run_split_reads.sh
```

#### Split the reads of each group into different cells using Indrop pipeline
+ Run the following to preprocess the reads
	
	```
		python $REPOHOME/indrops_V2/indrops/run.py 1
	```

  Then visually determine the reads cutoff to rule out phatom droplets (see [notebook](https://github.com/gifford-lab/lineage_tracing/blob/master/results/round2_Sept_2016/notebook.ipynb))

+ Change the cutoff in `$REPOHOME/indrops_V2/indrops/run.py` for each group accordingly, and run:

	```
		python $REPOHOME/indrops_V2/indrops/run.py 2
	```

#### Extract lineage barcodes

```
	python $REPOHOME/cell-lineage GROUPNAME TASK
```
+ `GROUPNAME`: one of the group names specified in `$REPOHOME/indrops_V2/indrops/run.py`
+ `TASK`: `sys` for integrated lineage barcode analysis and `let` for barcodelet analysis

Detailed pipeline:

1. Search lineage barcodes by matching the 22bp before and 11 bp after the barcode with **2** bp of mismatch allowed.  

2. In each barcode, any bases with quality score < **20** was labeled as 'N'. Any barcode with > **13** 'N' is _**discarded**_ as we can't be sure of its origin by mapping to other reads. (each barcode is 17 bp long)

3. Within each cell, the pair-wise Levenshtein distance of all the barcodes were calculated, only considering all the positions that is not 'N' for either barcodes.

4. In each cell, cluster the barcodes so that any pair-wise distance within one cluster is **<=2**.  After clustering, calculate the "true" barcode of each cluster in following manner: for each position in the barcode, pick the most commonly seen non-'N' nucleotide. Therefore, each cell is hence represented by this list of "true" barcode, one from each cluster.

5. For each cell, further _**trim**_ this list of "true" barcode by following two criteria to increase the mapping confidence in the next steps. 
	* with **<33%** of the 17bp as 'N'
	* with **>= 2** reads 


#### Further filter lineage barcodes and perform lineage analysis
Detailed implementation in the [notebook](https://github.com/gifford-lab/lineage_tracing/blob/master/results/round2_Sept_2016/notebook.ipynb) (function `analysis2` and `family2`)

1.	For each cell, first filter out all barcodes with less than 10 reads. Among the remaining barcodes, keep the barcode with the most reads and retrain the second ranking barcode too if it has more than half of the read count of the top ranking barcode.

2. Calculate the pair-wise (cell) distance of all the cells, where if at least one barcode in cell A and one barcode in cell B has a distance smaller than 0 with any base that is 'N' ignored,  the distance between these two cell is 0. Otherwise, it is 1.

3. Perform hierarchical clustering on the cell using the distance matrix obtained above so that the closest distance between two cluster >= 1.  The closest distance of two cluster is calculated as the smallest distance between any two cell in these two cluster.

## Results

#### Integrated barcodes
Here is a summary of the statistics of the four experiments:
![Integrated barcodes analysis](https://github.com/gifford-lab/lineage_tracing/blob/master/results/round2_Sept_2016/integrated.png?raw=true)

Some tables are put under folders corresponding to each of the [group](https://github.com/gifford-lab/lineage_tracing/tree/master/results/round2_Sept_2016/fth1_time1_rep2/)
+ \*.barcode: the unique barcodes and the number of cells they show up in before filtering (see the above section)
+ \*.barcode.afterprocessing: same as above but after filtering
+ \*.barcode.afterprocessing.bi: same as above but here for cells with two barcodes left after the filtering, we treat the combination of the barcodes as one 'barcode'.
+ \*.familyanalysis: the result of family analysis. Read the header for the meaning of each column.

#### Barcodelet
![Barcodelet analysis](https://github.com/gifford-lab/lineage_tracing/blob/master/results/round2_Sept_2016/barcodelet.png?raw=true)
