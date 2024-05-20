# rustymicrobial
Just a personal collection of bundled Rust tools to carry out microbial bioinformatics work.

To compile the program, make sure you have rust installed so you have 'cargo' in your path.
Then do the following:

```
git clone https://github.com/patrickmunk/rustymicrobial

cd rustymicrobial

cargo build --release

target/release/rustymicrobial --help
```
The following subcomands are available with distinct functionalies:

| Subcommand | Description |
| --- | --- |
| count | counts amino acids in protein multifasta |
| countkmers | counts kmers of length k in DNA multifasta |
| countnucleotides | counts kmers in each entry in a multifasta file |
| extractflanks | extract flanks of length k around exact substrings in a reference multifasta |
| help | help message for the program |
| mpt | metagenomic predicted temperature calculated on a proteome (a protein multifasta) |
| translate | translates a multifasta of genes into their corresponding proteins with the standard genetic code |
