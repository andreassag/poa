pub mod fasta;
pub mod graph;
pub mod align;
pub mod consensus;
pub mod export;

pub use fasta::{read_fasta, write_fasta, FastaRecord};
pub use graph::{PoaGraph, Edge};
pub use align::{align_to_graph, fuse_alignment, PoaParams, PoaAlignment};
pub use consensus::extract_consensus;
pub use export::{export_gfa, export_dot, export_msa};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_graph_initialization() {
        let mut graph = PoaGraph::new();
        let seq = b"ACGT";
        graph.init_with_sequence(seq);

        assert_eq!(graph.num_nodes(), 4);
        assert_eq!(graph.topo_order.len(), 4);
        assert_eq!(graph.bases, b"ACGT");
    }

    #[test]
    fn test_identical_sequence_fusion() {
        let mut graph = PoaGraph::new();
        let seq = b"ACGTACGT";
        graph.init_with_sequence(seq);

        let params = PoaParams::default();
        let aln = align_to_graph(&graph, seq, &params);
        fuse_alignment(&mut graph, seq, &aln, 1);

        // Identical sequence should reinforce existing nodes without adding new nodes
        assert_eq!(graph.num_nodes(), seq.len());
        assert_eq!(graph.num_sequences, 2);
    }

    #[test]
    fn test_bubble_creation() {
        let mut graph = PoaGraph::new();
        graph.init_with_sequence(b"ACGT");

        let mutated = b"ACTT"; // G -> T mismatch at pos 2
        let params = PoaParams::default();
        let aln = align_to_graph(&graph, mutated, &params);
        fuse_alignment(&mut graph, mutated, &aln, 1);

        // Should create 1 new sibling node for 'T'
        assert_eq!(graph.num_nodes(), 5);
        let cons = extract_consensus(&graph);
        assert!(!cons.is_empty());
    }

    #[test]
    fn test_consensus_majority() {
        let mut graph = PoaGraph::new();
        graph.init_with_sequence(b"AACCC");

        let s2 = b"AACCC";
        let s3 = b"AATCC"; // Single mismatch in s3

        let params = PoaParams::default();
        let aln2 = align_to_graph(&graph, s2, &params);
        fuse_alignment(&mut graph, s2, &aln2, 1);

        let aln3 = align_to_graph(&graph, s3, &params);
        fuse_alignment(&mut graph, s3, &aln3, 2);

        let cons = extract_consensus(&graph);
        // Majority is 'C' (2 vs 1)
        assert_eq!(cons, b"AACCC");
    }
}
