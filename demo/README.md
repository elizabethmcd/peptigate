We constructed demo data using the following commands:
```
curl -JLo Amblyomma_americanum_transcriptome_assembly_data.tar.gz https://zenodo.org/records/10574110/files/Amblyomma_americanum_transcriptome_assembly_data.tar.gz?download=1
tar xf Amblyomma_americanum_transcriptome_assembly_data.tar.gz
mv transcriptome_data/* .
head -n 200 orthofuser_final_clean.fa.transdecoder.pep > orfs_amino_acids.faa
head -n 200 orthofuser_final_clean.fa.dammit.fasta > contigs.fa
head -n 204 orthofuser_final_clean.fa.transdecoder.cds > orfs_nucleotides.fa
```

We also pulled short contigs (less than 75 bp) from an internal S3 bucket and added these contigs to the `contigs.fa` file (50 contigs).
These are contigs that were filtered from the *Amblyomma* transcriptome prior to transcriptome merging (mid assembly pipeline).
