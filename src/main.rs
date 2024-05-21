use structopt::StructOpt;
use bio::io::fasta;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;

#[derive(StructOpt)]
#[structopt(name = "rustymicro", about = "A suite of bioinformatics tools.")]
enum RustyMicro {
    #[structopt(name = "translate")]
    Translate {
        #[structopt(short = "i", long = "input", help = "Input file")]
        input: String,
        #[structopt(short = "o", long = "output", help = "Output file")]
        output: String,
    },
    #[structopt(name = "count")]
    Count {
        #[structopt(short = "i", long = "input", help = "Input file")]
        input: String,
        #[structopt(short = "o", long = "output", help = "Output file")]
        output: String,
    },
    #[structopt(name = "countnucleotides")]
    CountNucleotides {
        #[structopt(short = "i", long = "input", help = "Input file")]
        input: String,
        #[structopt(short = "o", long = "output", help = "Output file")]
        output: String,
    },
    #[structopt(name = "mpt")]
    MPT {
        #[structopt(short = "i", long = "input", help = "Input file")]
        input: String,
        #[structopt(short = "o", long = "output", help = "Output file")]
        output: String,
    },
    #[structopt(name = "countkmers")]
    CountKmers {
        #[structopt(short = "i", long = "input", help = "Input file")]
        input: String,
        #[structopt(short = "k", long = "kmer", help = "Kmer length")]
        kmer: usize,
        #[structopt(short = "o", long = "output", help = "Output file")]
        output: String,
    },
    #[structopt(name = "extractflanks")]
    ExtractFlanks {
        #[structopt(short = "q", long = "query", help = "Query file")]
        query: String,
        #[structopt(short = "r", long = "reference", help = "Reference file")]
        reference: String,
        #[structopt(short = "k", long = "kmer", help = "Kmer length")]
        kmer: usize,
        #[structopt(short = "o", long = "output", help = "Output file")]
        output: String,
    },
    #[structopt(name = "codonextractor")]
    CodonExtractor {
        #[structopt(short = "i", long = "input", help = "Input file")]
        input: String,
        #[structopt(short = "o", long = "output", help = "Output file")]
        output: String,
    },
    // Add other subcommands here...
}

fn main() {
    let opt = RustyMicro::from_args();

    match opt {
        RustyMicro::Translate { input, output } => {
            let reader = fasta::Reader::from_file(&input).unwrap();
            let mut writer = fasta::Writer::to_file(&output).unwrap();

            // Define the complete standard genetic code
            let genetic_code: HashMap<&str, &str> = [
                ("TTT", "F"), ("TTC", "F"), 
                ("TTA", "L"), ("TTG", "L"), ("CTT", "L"), ("CTC", "L"), ("CTA", "L"), ("CTG", "L"),
                ("ATT", "I"), ("ATC", "I"), ("ATA", "I"), 
                ("ATG", "M"),
                ("GTT", "V"), ("GTC", "V"), ("GTA", "V"), ("GTG", "V"),
                ("TCT", "S"), ("TCC", "S"), ("TCA", "S"), ("TCG", "S"), ("AGT", "S"), ("AGC", "S"),
                ("CCT", "P"), ("CCC", "P"), ("CCA", "P"), ("CCG", "P"),
                ("ACT", "T"), ("ACC", "T"), ("ACA", "T"), ("ACG", "T"),
                ("GCT", "A"), ("GCC", "A"), ("GCA", "A"), ("GCG", "A"),
                ("TAT", "Y"), ("TAC", "Y"),
                ("CAT", "H"), ("CAC", "H"),
                ("CAA", "Q"), ("CAG", "Q"),
                ("AAT", "N"), ("AAC", "N"),
                ("AAA", "K"), ("AAG", "K"),
                ("GAT", "D"), ("GAC", "D"),
                ("GAA", "E"), ("GAG", "E"),
                ("TGT", "C"), ("TGC", "C"),
                ("TGG", "W"),
                ("CGT", "R"), ("CGC", "R"), ("CGA", "R"), ("CGG", "R"), ("AGA", "R"), ("AGG", "R"),
                ("GGT", "G"), ("GGC", "G"), ("GGA", "G"), ("GGG", "G"),
                // ... include all 64 codons
            ].iter().cloned().collect();

            for result in reader.records() {
                let record = result.unwrap();
                let dna_seq = record.seq();
                let mut protein = String::new();

                // Translate each codon to an amino acid
                for codon in dna_seq.chunks(3) {
                    let codon_str = std::str::from_utf8(codon).unwrap();
                    if let Some(aa) = genetic_code.get(codon_str) {
                        protein.push_str(aa);
                    }
                }

                let result = fasta::Record::with_attrs(
                    record.id(),
                    Some(&format!("translation of {}", record.desc().unwrap_or(""))),
                    protein.as_bytes(),
                );

                writer.write_record(&result).unwrap();
            }
        },
        RustyMicro::Count { input, output } => {
            let reader = fasta::Reader::from_file(&input).unwrap();
            let mut writer = File::create(&output).unwrap();

            writeln!(writer, "name\tA\tR\tN\tD\tC\tQ\tE\tG\tH\tI\tL\tK\tM\tF\tP\tS\tT\tW\tY\tV").unwrap();

            for result in reader.records() {
                let record = result.unwrap();
                let protein_seq = record.seq();
                let mut counts = HashMap::new();

                for aa in protein_seq {
                    *counts.entry(*aa as char).or_insert(0) += 1;
                }

                write!(writer, "{}", record.id()).unwrap();
                for aa in "ARNDCQEGHILKMFPSTWYV".chars() {
                    write!(writer, "\t{}", counts.get(&aa).unwrap_or(&0)).unwrap();
                }
                writeln!(writer).unwrap();
            }
        },
        RustyMicro::CountNucleotides { input, output } => {
            let reader = fasta::Reader::from_file(&input).unwrap();
            let mut writer = File::create(&output).unwrap();

            writeln!(writer, "name\tA\tT\tC\tG").unwrap();

            for result in reader.records() {
                let record = result.unwrap();
                let dna_seq = record.seq();
                let mut counts = HashMap::new();

                for nt in dna_seq {
                    *counts.entry(*nt as char).or_insert(0) += 1;
                }

                write!(writer, "{}", record.id()).unwrap();
                for nt in "ATCG".chars() {
                    write!(writer, "\t{}", counts.get(&nt).unwrap_or(&0)).unwrap();
                }
                writeln!(writer).unwrap();
            }
        },
        RustyMicro::MPT { input, output } => {
            let reader = fasta::Reader::from_file(&input).unwrap();
            let mut writer = File::create(&output).unwrap();

            let mut total_aa_count = 0;
            let mut ivywrel_count = 0;

            for result in reader.records() {
                let record = result.unwrap();
                let protein_seq = record.seq();

                for aa in protein_seq {
                    let aa = *aa as char;
                    total_aa_count += 1;
                    if "IVYWREL".contains(aa) {
                        ivywrel_count += 1;
                    }
                }
            }

            let fraction_ivywrel = ivywrel_count as f64 / total_aa_count as f64;
            let mpt = 937.0 * fraction_ivywrel - 335.0;

            writeln!(writer, "{}\t{}", input, mpt).unwrap();
        },
        RustyMicro::CountKmers { input, kmer, output } => {
            let reader = fasta::Reader::from_file(&input).unwrap();
            let mut writer = File::create(&output).unwrap();

            // Generate all possible kmers of length k.
            let mut kmers = Vec::new();
            for i in 0..4usize.pow(kmer as u32) {
                let mut kmer_str = String::new();
                for j in 0..kmer {
                    let base = match (i / 4usize.pow(j as u32)) % 4 {
                        0 => 'A',
                        1 => 'T',
                        2 => 'C',
                        _ => 'G',
                    };
                    kmer_str.push(base);
                }
                kmers.push(kmer_str);
            }

            // Write the header to the output file.
            write!(writer, "name").unwrap();
            for kmer in &kmers {
                write!(writer, "\t{}", kmer).unwrap();
            }
            writeln!(writer).unwrap();

            for result in reader.records() {
                let record = result.unwrap();
                let dna_seq = record.seq();

                // Count the occurrence of each kmer in the sequence.
                let mut counts = HashMap::new();
                for i in 0..dna_seq.len() - kmer + 1 {
                    let kmer_str = std::str::from_utf8(&dna_seq[i..i+kmer]).unwrap();
                    *counts.entry(kmer_str).or_insert(0) += 1;
                }

                // Write the counts to the output file.
                write!(writer, "{}", record.id()).unwrap();
                for kmer in &kmers {
                    let kmer_str = kmer.as_str();
                    write!(writer, "\t{}", counts.get(kmer_str).unwrap_or(&0)).unwrap();
                }
                writeln!(writer).unwrap();
            }
        },
        RustyMicro::ExtractFlanks { query, reference, kmer, output } => {
            let query_reader = fasta::Reader::from_file(&query).unwrap();
            let reference_reader = fasta::Reader::from_file(&reference).unwrap();
            let mut writer = fasta::Writer::to_file(&output).unwrap();
        
            // Read the query sequences into a Vec.
            let mut queries = Vec::new();
            for result in query_reader.records() {
                let record = result.unwrap();
                queries.push(record.seq().to_owned());
            }
        
            // For each reference sequence, check if it contains any of the query sequences.
            for result in reference_reader.records() {
                let record = result.unwrap();
                let reference_seq = record.seq();
        
                for query in &queries {
                    if let Some(start) = reference_seq.windows(query.len()).position(|window| window == query.to_vec()) {
                        // Extract the flanks of length k on either side of the query sequence.
                        let left_flank_start = start.saturating_sub(kmer);
                        let right_flank_end = std::cmp::min(start + query.len() + kmer, reference_seq.len());
        
                        let flanks = &reference_seq[left_flank_start..right_flank_end];
        
                        // Write the flanks to the output file.
                        let flank_record = fasta::Record::with_attrs(record.id(), None, flanks);
                        writer.write_record(&flank_record).unwrap();
                    }
                }
            }
        },
        RustyMicro::CodonExtractor { input, output } => {
            let reader = fasta::Reader::from_file(&input).unwrap();
            let mut writer = File::create(&output).unwrap();

            // Generate all possible codons
            let bases = vec!['A', 'T', 'C', 'G'];
            let mut all_codons = Vec::new();
            for a in &bases {
                for b in &bases {
                    for c in &bases {
                        all_codons.push(format!("{}{}{}", a, b, c));
                    }
                }
            }

            writeln!(writer, "name\t{}", all_codons.join("\t")).unwrap();

            for result in reader.records() {
                let record = result.unwrap();
                let dna_seq = record.seq();
                let mut counts = HashMap::new();

                for codon in dna_seq.chunks(3) {
                    if codon.len() == 3 {
                        let codon_str = std::str::from_utf8(codon).unwrap();
                        *counts.entry(codon_str.to_string()).or_insert(0) += 1;
                    }
                }

                write!(writer, "{}", record.id()).unwrap();
                for codon in &all_codons {
                    write!(writer, "\t{}", counts.get(codon).unwrap_or(&0)).unwrap();
                }
                writeln!(writer).unwrap();
            }
        },
        // Handle other subcommands here...
    }
}
