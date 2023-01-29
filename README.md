# Single cell RNA-Seq data analysis pipeline with an example of reprocessing data from a study investigating susceptibility of multiple types of pancreatic islet cells to SARS-CoV-2.

Download cellranger

```
wget -O cellranger-6.1.2.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-6.1.2.tar.gz?Expires=1641366506&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci02LjEuMi50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2NDEzNjY1MDZ9fX1dfQ__&Signature=kaV8~ZabHhyDykUhbN~F78PDQfNZ64IamgsGc1nOSghFKPr0fbZ3WJk-2eWYh7IEt-KupenYP89W1zHi4lrxF~ZBbuP4NTaKEAa-G6ILJoX-VdyFnktkXFYDHgzEJ8ABq-NM6RWn20WD3a9BITNHTIWPtxjM-NaXAuR5uc5PuAEgjSDaQ2QBAQr~1q4aSM-~vJt~ia5e8acTz9RlM24EluLqfO59VCtAorP-5iJRwvLw9DjfrTlDtWfy3M2LSXp5OGmVJH1WUQReLK~0iZX2e8~vrHlAYpuxMa0Lgil6oHQ5s6vc~Dod3Aqpjb9sM~wuVo80zi4EqJ5nq0LU8SNbiQ__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
```

Download human reference genome data file

```
curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
```

Unzip cellranger

```
Tar -xf cellranger-7.1.0.tar.gz
```

Unzip human reference genome file
```
Tar -xf refdata-gex-GRCh38-2020-A.tar.gz
```

Download SRR12831415, SRR12831416,SRR12831417 and SRR12831418 files from GEO database with accession # GSE159556.
Place Hg38 reference genome folder, cellranger folder in the same folder, hereby named "scrnaseq". Within this folder, create a folder to store all downloaded fastq files.

Run cellranger count pipeline to align sequencing reads in FASTQ files to the reference transcriptome.
Set up the command to run cellranger count.

Specify an --id, which is the folder name for output files, here this directory is called "run_count". The --fastqs is the path to the directory containing the FASTQ files. Use the --sample argument to specify which samples to use. The --transcriptome argument specifies the path to the transcriptome reference package.
```
cd scrnaseq
Cellranger-7.1.0/bin/cellranger count --id=run_count --fastqs=fastqs --sample=SRR12831418 --transcriptome=refdata-gex-GRCh38-2020-A
```
When the output of the cellranger count command says, “Pipestance completed successfully!”, this means the job is done.
