# scRNA-seq pipeline (nextflow)

This pipeline uses **starsolo** (https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md) 
and **Seurat** package (https://satijalab.org/seurat/) as the core modules to process scRNA-seq
data and generate a concise html report. The workflow was implemented by [**nextflow**](https://www.nextflow.io/).

Please visit our documentation website for more details: https://starscope.netlify.app

## Compatible Data Type

- 10X Genomics Chromium (https://www.10xgenomics.com/resources/datasets)
- ThunderBio scRNA-seq data
- Dropseq

## Usage

```
nextflow run /path/to/pipeline/dir --input sampleList.csv -c custom.config [-bg] [-resume]
```

### Install Nextflow

Please refer to: https://www.nextflow.io/

- Ensure Java 8 or later is installed: `java --version`

- Installation command: `curl -s https://get.nextflow.io | bash`

- Test with: `./nextflow run hello`

### Prepare Running Environment

The conda env was demonstrated here. One could also prepare docker image(s).

- Create conda env from file

```
conda env create -f scRNAseq_env.yml
```

### Prepare Reference

To make the output "consistent" with CellRanger result, one could 
use reference dataset prepared by 10X directly. 
Please visit their download page (https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest) to 
obtain the latest version. It's worthy noting 10X conducted some custom filtration steps based on ensembl reference. 
Please visit the [building notes](https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#grch38_#{files.refdata_GRCh38.version}) for details.


Or one could generate their own reference by STAR as long as they have genome fasta file and 
gene annotation file (GTF), the command was suggested 
in [`starsolo`'s README file](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md).

```
STAR  --runMode genomeGenerate --runThreadN ... --genomeDir ./ --genomeFastaFiles /path/to/genome.fa  --sjdbGTFfile /path/to/genes.gtf
```

A `mkref` sub-workflow was added and requires `--genomeFasta`, `--genomeGTF`, `--refoutDir` options to run. A typical
using example would be:

```
nextflow run /path/to/scRNA-seq -entry mkref --genomeFasta /path/to/genome.fa --genomeGTF /path/to/genes.gtf --refoutDir ref_out_dir -bg
```

### Prepare Whitelist

As suggested in `starsolo`'s document, the white list of 10X barcode could be accessed as:

> The 10X Chromium whitelist files can be found or inside the CellRanger distribution or on [GitHub/10XGenomics](https://github.com/10XGenomics/cellranger/tree/master/lib/python/cellranger/barcodes). Please make sure that the whitelist is compatible with the specific version of the 10X chemistry: V2 or V3.

For V2:

```
cellranger-cs/3.1.0/lib/python/cellranger/barcodes/737K-august-2016.txt

https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-august-2016.txt
```

For V3:

```
cellranger-cs/3.1.0/lib/python/cellranger/barcodes/3M-february-2018.txt.gz

https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz
```

### Prepare Sample List

The workflow takes a csv format sample list as input, which is the `sampleList.csv` in the command above. The sample list file only contains three column: **sample**, **fastq_1**, **fastq_2**. User only needs to specify sample name, path of read1 and read2. The program will assume that read1 file is the one containing barcode and UMI information. Both absolute and relative path are supported. If the library has multiple runs (i.e. multiple files of R1 and R2), user could indicate the path of those files with identical sample name, and they will be `cat` into one in the analysis.

```
sample,fastq_1,fastq_2
10_pbmc_1k,read1.fq.gz,/absolute/path/to/read2.fq.gz
```

### Local Config File

Example custom config file for 10x V3 datasets, please adjust executor parameter according to your setup. For executors of nextflow, please refer to its [doc page](https://www.nextflow.io/docs/latest/executor.html). The cell barcode and UMI position were indicated by the `soloCBstart` , `soloCBlen`, `soloUMIstart` and  `soluUMIlen`.

```
params {
  genomeDir = "/path/to/STAR/reference/"
  genomeGTF = "/path/to/reference/genes.gtf"
  whitelist = "/path/to/whitelist"
  soloCBstart = 1
  soloCBlen = 16
  soloUMIstart = 17
  soloUMIlen = 12
  enable_conda = true
}

process {
  executor = "slurm" // remove this line if use local executor
  conda = "/path/to/miniconda3/envs/starscope_env
  // adjust resources here
  withLabel: process_high {
    cpus = 8
    memory = 32.GB
  }
  withLabel: process_medium {
    cpus = 4
    memory = 20.GB
  }
  withLabel: process_low {
    cpus = 4
    memory = 20.GB
  }
}
```

## Output

After finished running, user could find all "published" result files in `results` directory.
A tsv file recording execution trace could be found as `results/pipeline_info/execution_trace_{timeStamp}.txt`.
If one found all processes were finished as "COMPLETED" status, the intermediate
directory `work` could be safely removed by `rm -rf work`.

- **results**: main output directory
  - **pipeline\_info**: store plain text execution_trace file
  - **cutqc**: store cutadapt trimming report generated by customized R markdown script
  - **starsolo**: store `starsolo` output files, including barcode/UMI matrix, summary files etc.
  - **qualimap**: store `qualimap` report
  - **report**: contain main HTML report files

- **work**: intermediate directory, could be safely removed after running complete

