topdir=/cluster/shortreads/reads/160908_lineagetracing_zeng/160902_NS500414_0227_HLLWYBGXY/Data/Intensities/BaseCalls
python split_reads.py -f $topdir/all_R1.fastq -g $topdir/all_R2.fastq  -i ../index_analysis/index_barcodes -o split_result_hard_le2 -c 2
