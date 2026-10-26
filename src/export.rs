use crate::fasta::FastaRecord;
use crate::graph::PoaGraph;
use std::fs::File;
use std::io::Write;

pub fn export_gfa(graph: &PoaGraph, path: &str) -> Result<(), String> {
    let mut file = File::create(path).map_err(|e| format!("Cannot create {}: {}", path, e))?;

    writeln!(file, "H\tVN:Z:1.0").map_err(|e| e.to_string())?;

    for i in 0..graph.num_nodes() {
        let b = graph.bases[i] as char;
        // Calculate coverage as max of in/out weights or 1
        let mut cov = 1;
        for e in &graph.out_edges[i] {
            cov = cov.max(e.weight);
        }
        for e in &graph.in_edges[i] {
            cov = cov.max(e.weight);
        }
        writeln!(file, "S\ts{}\t{}\tLN:i:1\tic:i:{}", i, b, cov).map_err(|e| e.to_string())?;
    }

    for u in 0..graph.num_nodes() {
        for e in &graph.out_edges[u] {
            writeln!(file, "L\ts{}\t+\ts{}\t+\t0M\trc:i:{}", u, e.target, e.weight)
                .map_err(|e| e.to_string())?;
        }
    }

    Ok(())
}

pub fn export_dot(graph: &PoaGraph, path: &str) -> Result<(), String> {
    let mut file = File::create(path).map_err(|e| format!("Cannot create {}: {}", path, e))?;

    writeln!(file, "digraph POA {{").map_err(|e| e.to_string())?;
    writeln!(file, "  rankdir=LR;").map_err(|e| e.to_string())?;
    writeln!(file, "  node [shape=circle, fontname=\"Courier\"];").map_err(|e| e.to_string())?;

    for i in 0..graph.num_nodes() {
        let b = graph.bases[i] as char;
        let color = match b {
            'A' => "#4285F4", // Blue
            'C' => "#EA4335", // Red
            'G' => "#FBBC05", // Yellow
            'T' => "#34A853", // Green
            _ => "#9AA0A6",
        };
        writeln!(
            file,
            "  n{} [label=\"{}\", color=\"{}\", style=\"filled\", fillcolor=\"{}22\"];",
            i, b, color, color
        )
        .map_err(|e| e.to_string())?;
    }

    for u in 0..graph.num_nodes() {
        for e in &graph.out_edges[u] {
            let width = (e.weight as f32).max(1.0);
            writeln!(
                file,
                "  n{} -> n{} [label=\"{}\", penwidth={}];",
                u, e.target, e.weight, width
            )
            .map_err(|e| e.to_string())?;
        }
    }

    writeln!(file, "}}").map_err(|e| e.to_string())?;
    Ok(())
}

pub fn export_msa(
    graph: &PoaGraph,
    records: &[FastaRecord],
    consensus: &[u8],
    path: &str,
) -> Result<(), String> {
    let mut file = File::create(path).map_err(|e| format!("Cannot create {}: {}", path, e))?;

    let consensus_nodes = crate::consensus::extract_consensus_nodes(graph);
    let mut is_consensus_node = vec![false; graph.num_nodes()];
    for &u in &consensus_nodes {
        is_consensus_node[u as usize] = true;
    }

    // Export row-aligned consensus
    writeln!(file, ">Consensus").map_err(|e| e.to_string())?;
    let mut cons_row = String::new();
    for &u in &graph.topo_order {
        if is_consensus_node[u as usize] {
            cons_row.push(graph.bases[u as usize] as char);
        } else {
            cons_row.push('-');
        }
    }
    writeln!(file, "{}", cons_row).map_err(|e| e.to_string())?;

    // Determine nodes visited by each input sequence
    for (s_id, r) in records.iter().enumerate() {
        let sid = s_id as u32;
        let mut visited = vec![false; graph.num_nodes()];

        for u in 0..graph.num_nodes() {
            for e in &graph.out_edges[u] {
                if e.seq_ids.contains(&sid) {
                    visited[u] = true;
                    visited[e.target as usize] = true;
                }
            }
        }
        // Sequence 0 starts at node 0
        if sid == 0 && graph.num_nodes() > 0 {
            visited[0] = true;
        }

        let mut row = String::new();
        for &u in &graph.topo_order {
            if visited[u as usize] {
                row.push(graph.bases[u as usize] as char);
            } else {
                row.push('-');
            }
        }

        writeln!(file, ">{}", r.name).map_err(|e| e.to_string())?;
        writeln!(file, "{}", row).map_err(|e| e.to_string())?;
    }

    Ok(())
}
