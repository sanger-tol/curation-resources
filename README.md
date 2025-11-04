# Rapid curation - Curation tools guide

This section will cover the curation tools and how they work.



## Support 

[Slack Channel for curation assembly help](https://join.slack.com/t/assemblycuration/shared_invite/zt-1kx2ww71y-823ruaAxswgQGypgofBaOA)


PretextView tutorial:
* [YouTube video 1](https://youtu.be/3IL2Q4f3k3I)
* [BGA24: Introduction to Manual Curation](https://www.youtube.com/watch?v=qlWZlFXy07w)


PretextView can be obtained here https://github.com/wtsi-hpag/PretextView/releases


To generate the analysis data used for curation run one of the nextflow pipelines:

https://pipelines.tol.sanger.ac.uk/treeval

Full pipeline above can include files for displaying in JBrowse2 as a companion to the Hi-C 2D contact map


https://pipelines.tol.sanger.ac.uk/curationpretext




## Training data

Pretext maps and assembly files can be found at [the Curation Resources bucket](https://curation_resources.cog.sanger.ac.uk/) 

A training environment can be accessed (requires github account) that contains all the scripts and data, feel free to have a play!

[![Open in Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/sanger-tol/rapid-curation)


Please also see documentation:

[RAPID_CURATION_TRAINING_MANUAL](-/blob/main/RAPID_CURATION_TRAINING_MANUAL.pdf)

[Interpreting_HiC_Maps_guide](-/blob/main/Interpreting_HiC_Maps_guide.pdf)

 
## Suggested workflow

Manipulate the assembly in PretextView app using Edit mode (E).
Tag unlocalised scaffolds, haplotigs, contaminats, haplotypes etc, using meta-data tag mode (M) in PretextView.
Unlocs (linked to chromosome but no specific location found) should be placed at the start/and/or/end 
of a chromosome to be included when the chromosomes are painted in the subsequent step.

Once the map has been rearranged and meta-data tags added paint the chromosomes using Scaffold Edit mode (S).

Export the map to AGP. (Some input scaffolds won't feature in the output AGP due to PretextView resolution.  
These will be recovered in the next step).

Run [pretext-to-tasm](https://github.com/sanger-tol/agp-tpf-utils.git) , pointing it at the original FASTA and the Pretext AGP. 
This will output a FASTA that mirrors the chromosomes that have been built in PretextView, and which will have the same number 
of basepairs as the original TPF. The script will produce stats (joins/breaks/hap-dup removals).

Produce a new Pretext map from the curated fasta file to check it.


## Scripts

**Scripts and documentation for Rapid curation:**

###  [pretext-to-asm](https://github.com/sanger-tol/agp-tpf-utils.git)
to install:
```
git clone https://github.com/sanger-tol/agp-tpf-utils.git
```

This script takes the AGP output from PretextView and the original FASTA to generate new assembly fastas and supporting files

Usage: 

```
pretext-to-asm [OPTIONS]

  Uses fragments in the assembly (AGP) produced by PretextView to find
  matching fragments in the assembly which was fed into Pretext and output an
  assembly made from the input assembly fragments.

  Named Chromsomes

    Upper case letters followed by zero or more digits are assumed to be
    chromosome names. e.g. 'X', 'W', 'B1'

  Known Tags

    Contaminant tagged scaffolds are saved in a separate
    'Contaminants' file.

    When there are large numbers of contaminant scaffolds in the   assembly,
    Target tags can insted be used to label the   non-contaminant
    scaffolds and reduce the amount of labelling   necessary in PretextView.
    Any un-tagged scaffolds will then be   treated as if they were tagged with
    Contaminant.   (Any contaminants occurring before the first
    Target tag in   the PretextView AGP must still be individually
    tagged with   Contaminant.)

    Haplotig taggged scaffolds are saved in a separate 'Haplotigs'
    file. The haplotig scaffolds receive names 'H_1' to 'H_n', sorted
    and numbered from longest to shortest.

    Unloc tagged scaffolds receive names 'CHR_unloc_1' to
    'CHR_unloc_n', added to the end of their chromosome and
    sorted and numbered from longest to shortest.

  Haplotypes

    Any other tags are assumed to be the name of a haplotype, and their
    assemblies are placed in separate files. Unplaced scaffolds for each
    haplotype are identified by their names beginning with the haplotype's
    name followed by an underscore. i.e. 'Hap2_' for 'Hap2'

Options:
  -a, --assembly PATH             Assembly before curation, usually a FASTA
                                  file. FASTA files will be indexed, creating
                                  a '.fai' and a '.agp' file alongside the
                                  assembly if they are missing or are older
                                  than the FASTA.  [required]
  -p, --pretext PATH              Assembly file from Pretext, which is usually
                                  an AGP.  [required]
  -o, --output FILE               Output file, usually a FASTA file. If not
                                  given, prints to STDOUT in 'STR' format. The
                                  output file type is determined from its
                                  extension. If the outuput is FASTA ('.fa'),
                                  an AGP format file ('.fa.agp') is also
                                  written. Other output files are named after
                                  the output file minus its extension.
  -c, --autosome-prefix TEXT      Prefix for naming autosomal chromosomes.
                                  [default: SUPER_]
  -f, --clobber / --no-clobber    Overwrite an existing output file.
                                  [default: clobber]
  -l, --log-level [DEBUG|INFO|WARNING|ERROR|CRITICAL]
                                  Diagnostic messages to show.  [default:
                                  INFO]
  -w, --write-log / -W, --no-write-log
                                  Write messages into a '.log' file alongside
                                  the output file  [default: write-log]
  --help                          Show this message and exit.

```

See https://github.com/sanger-tol/agp-tpf-utils# for complete documentation

### Other scripts

Telomere identification script [telo_finder.py](-/blob/main/telo_finder.py):

The telo_finder script aims to detect likely telomere motifs in assemblies.
It assumes the assemblies are high quality (hi-fi or similar) and that there are a number of chromosomes (it won't return a result with just one chromosome for example).
It uses a hard-coded set of non-redundant telomere motifs with which it attempts to distinguish true telomere sequences from background noise.


Usage: 
```
telo_finder.py [-h] [--size size] [--klo klo] [--khi khi] [--ends ends] fasta

Finds most likely telomere motif in Hi-fi or equivalent quality assembly where
telomeres are expected to occur at the ends of multiple scaffolds

positional arguments:
  fasta        assembly fasta

optional arguments:
  -h, --help   show this help message and exit
  --size size  top and tail scf bp (default: 200)
  --klo klo    min kmer (default: 4)
  --khi khi    max kmer (default: 15)
  --ends ends  ends to scan (default: 1000)
  
```

