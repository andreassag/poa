use crate::graph::PoaGraph;

pub fn extract_consensus_nodes(graph: &PoaGraph) -> Vec<u32> {
    let n = graph.num_nodes();
    if n == 0 {
        return Vec::new();
    }

    let mut score = vec![0u64; n];
    let mut parent: Vec<Option<u32>> = vec![None; n];

    // Dynamic programming over topological order
    for &u in &graph.topo_order {
        let u_idx = u as usize;
        let mut best_score = 0;
        let mut best_pred = None;

        for e in &graph.in_edges[u_idx] {
            let p = e.target as usize;
            let edge_score = score[p] + (e.weight as u64);
            if edge_score > best_score {
                best_score = edge_score;
                best_pred = Some(e.target);
            }
        }

        score[u_idx] = best_score;
        parent[u_idx] = best_pred;
    }

    // Find end node with maximum score among terminal nodes
    let mut best_end = 0;
    let mut max_total = 0;
    for i in 0..n {
        if graph.out_edges[i].is_empty() && score[i] >= max_total {
            max_total = score[i];
            best_end = i;
        }
    }

    if max_total == 0 {
        // Fallback: overall max score node
        for i in 0..n {
            if score[i] >= max_total {
                max_total = score[i];
                best_end = i;
            }
        }
    }

    // Backtrace consensus path
    let mut nodes = Vec::new();
    let mut cur: Option<u32> = Some(best_end as u32);

    while let Some(u) = cur {
        nodes.push(u);
        cur = parent[u as usize];
    }

    nodes.reverse();
    nodes
}

pub fn extract_consensus(graph: &PoaGraph) -> Vec<u8> {
    let nodes = extract_consensus_nodes(graph);
    nodes.iter().map(|&u| graph.bases[u as usize]).collect()
}
