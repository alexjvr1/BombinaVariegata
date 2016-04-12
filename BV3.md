#Analysis of BV3

BV3 data has just arrived (12 April 2016) - i.e. 8 days after submission. 

1. FastQC
2. Demultiplex 
3. pyRAD for entire dataset (with or without concat data??)



##1. FastQC

The data quality looks good, but we only got 160Mil reads! The HiSeq 2500 can sequence 250Mil reads, so this is well below what we would expect. I will contact Cathy and ask her about it. 


![alt_txt][stats]
[stats]:https://cloud.githubusercontent.com/assets/12142475/14467673/3ede09c4-0091-11e6-8009-89601dc69b48.png


![alt_txt][base.qual]
[base.qual]:https://cloud.githubusercontent.com/assets/12142475/14467672/3edd88aa-0091-11e6-99d7-a61dd97e3f74.png

![alt_txt][Fig3]
[Fig3]:https://cloud.githubusercontent.com/assets/12142475/14467675/3edfa05e-0091-11e6-966c-bd6be1c20670.png

![alt_txt][Fig4]
[Fig4]:https://cloud.githubusercontent.com/assets/12142475/14467677/3ee11ec0-0091-11e6-87ca-4d94870cb3a2.png

![alt_txt][Fig5]
[Fig5]:https://cloud.githubusercontent.com/assets/12142475/14467674/3edf9a00-0091-11e6-9733-fb5be86e2bf4.png

![alt_txt][Fig6]
[Fig6]:https://cloud.githubusercontent.com/assets/12142475/14467676/3ee01ffc-0091-11e6-92c7-f2a0c6bf4deb.png



##2. Demultiplex the data

I will set up the barcode file and demultiplex on FGCZ server as before. 

Barcodes: 

![alt_txt][barcodes]
[barcodes]:https://cloud.githubusercontent.com/assets/12142475/14468112/28af11aa-0093-11e6-8b45-d2286752321a.png


Demultiplex: 

Started 12 April 9:51

```
/usr/local/ngseq/stow/stacks-1.28/bin/process_radtags -i gzfastq -f /srv/gstore4users/p1795/HiSeq2500_20160405_RUN263_o2409_DataDelivery/20160405.A-BV3_R1.fastq.gz  -o ./demultiplexed -y fastq -b barcodes_BV3 --disable_rad_check -r -D


166256011 total sequences;
  28165822 ambiguous barcode drops;
  0 low quality read drops;
  0 ambiguous RAD-Tag drops;
138090189 retained reads.
```


And trim data: 

```
screen -S TrimSubset -L
for i in *.fq; do  java -jar /usr/local/ngseq/src/Trimmomatic-0.33/trimmomatic-0.33.jar SE $i $i.trim ILLUMINACLIP:/usr/local/ngseq/src/Trimmomatic-0.33/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36;done
```

Demultiplex and Trim complete by 11:16, 12 April 2016. 

I will copy the data to the GDCserver and start the pyRAD run, including all of the samples (i.e. BV2 + BV3)

From the GDC side
```
rsync -av -e "ssh -l alexjvr" /srv/kenlab/alexjvr_p1795/Bombina/BV3/demultiplexed/* gdcsrv1.ethz.ch:/gdc_home4/alexjvr/Bombina/BV3/demultiplexed.BV3/
```









