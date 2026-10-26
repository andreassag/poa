use std::env;
use std::process;
use poa::{
    align_to_graph, export_dot, export_gfa, export_msa, extract_consensus, fuse_alignment,
    read_fasta, write_fasta, FastaRecord, PoaGraph, PoaParams,
};

fn print_usage() {
    eprintln!(
        "Usage: poa [options] -i <input.fa>\n\n\
         Options:\n\
           -i, --input <file>           Input multi-FASTA file (required)\n\
           -c, --consensus-out <file>   Write consensus sequence to FASTA\n\
           -m, --msa-out <file>         Write multiple sequence alignment to FASTA\n\
           -g, --gfa-out <file>         Export variation graph to GFA v1\n\
           -d, --dot-out <file>         Export variation graph to Graphviz DOT\n\
           --mode <global|local>        Alignment mode (default: global)\n\
           --match <int>                Match score (default: 2)\n\
           --mismatch <int>             Mismatch penalty (default: -4)\n\
           --gap-open <int>             Gap open penalty (default: -6)\n\
           --gap-extend <int>           Gap extend penalty (default: -2)\n"
    );
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        print_usage();
        process::exit(1);
    }

    let mut in_path = String::new();
    let mut consensus_out = String::new();
    let mut msa_out = String::new();
    let mut gfa_out = String::new();
    let mut dot_out = String::new();
    let mut params = PoaParams::default();

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "-i" | "--input" => {
                i += 1;
                if i < args.len() { in_path = args[i].clone(); }
            }
            "-c" | "--consensus-out" => {
                i += 1;
                if i < args.len() { consensus_out = args[i].clone(); }
            }
            "-m" | "--msa-out" => {
                i += 1;
                if i < args.len() { msa_out = args[i].clone(); }
            }
            "-g" | "--gfa-out" => {
                i += 1;
                if i < args.len() { gfa_out = args[i].clone(); }
            }
            "-d" | "--dot-out" => {
                i += 1;
                if i < args.len() { dot_out = args[i].clone(); }
            }
            "--mode" => {
                i += 1;
                if i < args.len() && args[i] == "local" {
                    params.is_local = true;
                }
            }
            "--match" => {
                i += 1;
                if i < args.len() { params.match_score = args[i].parse().unwrap_or(2); }
            }
            "--mismatch" => {
                i += 1;
                if i < args.len() { params.mismatch_score = args[i].parse().unwrap_or(-4); }
            }
            "--gap-open" => {
                i += 1;
                if i < args.len() { params.gap_open = args[i].parse().unwrap_or(-6); }
            }
            "--gap-extend" => {
                i += 1;
                if i < args.len() { params.gap_extend = args[i].parse().unwrap_or(-2); }
            }
            _ => {}
        }
        i += 1;
    }

    if in_path.empty() {
        eprintln!("Error: -i <input.fa> is required.");
        print_usage();
        process::exit(1);
    }

    let records = read_fasta(&in_path).unwrap_or_else(|e| {
        eprintln!("Error reading FASTA: {}", e);
        process::exit(1);
    });

    if records.is_empty() {
        eprintln!("Error: No FASTA sequences found in {}", in_path);
        process::exit(1);
    }

    println!("[poa] Initializing POA graph with sequence 1 ({} bp)...", records[0].seq.len());
    let mut graph = PoaGraph::new();
    graph.init_with_sequence(&records[0].seq);

    for (idx, r) in records.iter().enumerate().skip(1) {
        println!("[poa] Aligning sequence {} ({}, {} bp) to graph...", idx + 1, r.name, r.seq.len());
        let aln = align_to_graph(&graph, &r.seq, &params);
        fuse_alignment(&mut graph, &r.seq, &aln, idx as u32);
    }

    println!("[poa] Built POA graph with {} nodes and {} sequences.", graph.num_nodes(), graph.num_sequences);

    let consensus = extract_consensus(&graph);
    println!("[poa] Extracted consensus sequence ({} bp):", consensus.len());
    let cons_str = String::from_utf8_lossy(&consensus);
    if cons_str.len() <= 100 {
        println!("  {}", cons_str);
    } else {
        println!("  {}...", &cons_str[..100]);
    }

    if !consensus_out.is_empty() {
        let cons_record = FastaRecord {
            name: "POA_Consensus".to_string(),
            seq: consensus.clone(),
        };
        write_fasta(&consensus_out, &[cons_record]).unwrap();
        println!("[poa] Saved consensus to {}", consensus_out);
    }

    if !msa_out.is_empty() {
        export_msa(&graph, &records, &consensus, &msa_out).unwrap();
        println!("[poa] Saved MSA to {}", msa_out);
    }

    if !gfa_out.is_empty() {
        export_gfa(&graph, &gfa_out).unwrap();
        println!("[poa] Exported GFA v1 to {}", gfa_out);
    }

    if !dot_out.is_empty() {
        export_dot(&graph, &dot_out).unwrap();
        println!("[poa] Exported Graphviz DOT to {}", dot_out);
    }
}
