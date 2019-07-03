# WWP_Soil
Woolsey Wet Prairie Soil - 16S and ITS1

The phyloseq_object_* files are from the full MiSeq run.
02_WWP_phyloseq_subsetting.R was used to subset from full objects to WWP samples.
WWP* phyloseq.RDS are the WWP phyloseq objects.

./output/ contains all intermediate files and statistical tests

./output/figs/ contains figures


I ran ITS fwd reads through ITSxpress before running them through DADA2 with the following:
for fn in ./ITS/*.fastq.gz;do echo $fn; itsxpress --fastq $fn -s --outfile ./ITS/$fn.ITS1 --region ITS1 --taxa Fungi --log ./ITS/$fn.log --threads 10; done

metadata.xlsx was the full Illumina run metadata, including other projects. Need to copy it over here from office machine.
