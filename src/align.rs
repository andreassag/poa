use crate::graph::PoaGraph;
use std::cmp::max;

#[derive(Clone, Copy, Debug)]
pub struct PoaParams {
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub is_local: bool,
}

impl Default for PoaParams {
    fn default() -> Self {
        Self {
            match_score: 2,
            mismatch_score: -4,
            gap_open: -6,
            gap_extend: -2,
            is_local: false,
        }
    }
}

pub struct PoaAlignment {
    pub score: i32,
    pub path: Vec<(Option<u32>, Option<usize>)>, // (graph_node, query_pos)
}

pub fn align_to_graph(graph: &PoaGraph, query: &[u8], params: &PoaParams) -> PoaAlignment {
    let num_nodes = graph.num_nodes();
    let m = query.len();
    if num_nodes == 0 || m == 0 {
        return PoaAlignment {
            score: 0,
            path: Vec::new(),
        };
    }

    let inf = 1_000_000_000;
    // DP matrices indexed by [node][query_pos] (0..=m)
    let mut dp_m = vec![vec![-inf; m + 1]; num_nodes];
    let mut dp_d = vec![vec![-inf; m + 1]; num_nodes];
    let mut dp_i = vec![vec![-inf; m + 1]; num_nodes];

    // Topological evaluation
    for &u in &graph.topo_order {
        let u_idx = u as usize;
        let u_base = graph.bases[u_idx];

        // Query pos 0:
        if graph.in_edges[u_idx].is_empty() {
            dp_m[u_idx][0] = if params.is_local { 0 } else { params.gap_open + params.gap_extend };
            dp_d[u_idx][0] = dp_m[u_idx][0];
        } else {
            let mut best_m0 = -inf;
            let mut best_d0 = -inf;
            for e in &graph.in_edges[u_idx] {
                let p = e.target as usize;
                best_d0 = max(best_d0, max(dp_m[p][0] + params.gap_open + params.gap_extend, dp_d[p][0] + params.gap_extend));
                best_m0 = max(best_m0, best_d0);
            }
            dp_d[u_idx][0] = best_d0;
            dp_m[u_idx][0] = if params.is_local { max(0, best_m0) } else { best_m0 };
        }

        for j in 1..=m {
            let q_base = query[j - 1];
            let cost = if u_base == q_base { params.match_score } else { params.mismatch_score };

            // Transition from predecessor nodes
            if graph.in_edges[u_idx].is_empty() {
                // Starting node: transition from query insertion
                let init_pred = if params.is_local { 0 } else { params.gap_open + (j as i32) * params.gap_extend };
                dp_m[u_idx][j] = init_pred + cost;
            } else {
                let mut best_pred = -inf;
                let mut best_del = -inf;
                for e in &graph.in_edges[u_idx] {
                    let p = e.target as usize;
                    let p_max = max(dp_m[p][j - 1], max(dp_d[p][j - 1], dp_i[p][j - 1]));
                    best_pred = max(best_pred, p_max);

                    best_del = max(best_del, max(dp_m[p][j] + params.gap_open + params.gap_extend, dp_d[p][j] + params.gap_extend));
                }
                dp_m[u_idx][j] = best_pred + cost;
                dp_d[u_idx][j] = best_del;
            }

            // Insertion in query
            dp_i[u_idx][j] = max(
                dp_m[u_idx][j - 1] + params.gap_open + params.gap_extend,
                dp_i[u_idx][j - 1] + params.gap_extend,
            );

            if params.is_local {
                dp_m[u_idx][j] = max(0, dp_m[u_idx][j]);
            }
        }
    }

    // Find best terminal node
    let mut best_score = -inf;
    let mut best_node = 0;
    let mut best_j = m;

    for u in 0..num_nodes {
        let score = max(dp_m[u][m], max(dp_d[u][m], dp_i[u][m]));
        if score > best_score {
            best_score = score;
            best_node = u;
        }
    }

    // Backtrace
    let mut path = Vec::new();
    let mut cur_u: Option<usize> = Some(best_node);
    let mut cur_j = best_j;

    let mut cur_state = if dp_m[best_node][best_j] >= dp_d[best_node][best_j] && dp_m[best_node][best_j] >= dp_i[best_node][best_j] {
        'M'
    } else if dp_d[best_node][best_j] >= dp_i[best_node][best_j] {
        'D'
    } else {
        'I'
    };

    while let Some(u) = cur_u {
        if cur_state == 'M' {
            if cur_j == 0 {
                path.push((Some(u as u32), None));
                break;
            }
            path.push((Some(u as u32), Some(cur_j - 1)));

            if graph.in_edges[u].is_empty() {
                break;
            }

            let mut pred_u = None;
            let mut next_state = 'M';
            let mut max_val = -inf;

            for e in &graph.in_edges[u] {
                let p = e.target as usize;
                let j_prev = cur_j - 1;
                if dp_m[p][j_prev] > max_val {
                    max_val = dp_m[p][j_prev];
                    pred_u = Some(p);
                    next_state = 'M';
                }
                if dp_d[p][j_prev] > max_val {
                    max_val = dp_d[p][j_prev];
                    pred_u = Some(p);
                    next_state = 'D';
                }
                if dp_i[p][j_prev] > max_val {
                    max_val = dp_i[p][j_prev];
                    pred_u = Some(p);
                    next_state = 'I';
                }
            }

            cur_u = pred_u;
            cur_j -= 1;
            cur_state = next_state;
        } else if cur_state == 'D' {
            path.push((Some(u as u32), None));

            if graph.in_edges[u].is_empty() {
                break;
            }

            let mut pred_u = None;
            let mut next_state = 'M';
            let mut max_val = -inf;

            for e in &graph.in_edges[u] {
                let p = e.target as usize;
                if dp_m[p][cur_j] + params.gap_open + params.gap_extend > max_val {
                    max_val = dp_m[p][cur_j] + params.gap_open + params.gap_extend;
                    pred_u = Some(p);
                    next_state = 'M';
                }
                if dp_d[p][cur_j] + params.gap_extend > max_val {
                    max_val = dp_d[p][cur_j] + params.gap_extend;
                    pred_u = Some(p);
                    next_state = 'D';
                }
            }

            cur_u = pred_u;
            cur_state = next_state;
        } else {
            // cur_state == 'I'
            if cur_j == 0 {
                break;
            }
            path.push((None, Some(cur_j - 1)));

            if dp_m[u][cur_j - 1] + params.gap_open + params.gap_extend >= dp_i[u][cur_j - 1] + params.gap_extend {
                cur_state = 'M';
            } else {
                cur_state = 'I';
            }
            cur_j -= 1;
        }
    }

    path.reverse();

    PoaAlignment {
        score: best_score,
        path,
    }
}

pub fn fuse_alignment(graph: &mut PoaGraph, query: &[u8], alignment: &PoaAlignment, seq_id: u32) {
    if alignment.path.is_empty() {
        return;
    }

    let mut prev_node: Option<u32> = None;

    for &(opt_node, opt_q) in &alignment.path {
        let cur_node = match (opt_node, opt_q) {
            (Some(node_id), Some(q_idx)) => {
                let q_base = query[q_idx];
                if graph.bases[node_id as usize] == q_base {
                    node_id
                } else {
                    // Check if sibling ring already contains this base
                    let mut found_sib = None;
                    let mut sib = graph.aligned_sibling[node_id as usize];
                    while sib != node_id {
                        if graph.bases[sib as usize] == q_base {
                            found_sib = Some(sib);
                            break;
                        }
                        sib = graph.aligned_sibling[sib as usize];
                    }

                    if let Some(s) = found_sib {
                        s
                    } else {
                        // Allocate new parallel sibling node
                        let new_id = graph.add_node(q_base);
                        graph.link_sibling(node_id, new_id);
                        new_id
                    }
                }
            }
            (Some(node_id), None) => node_id,
            (None, Some(q_idx)) => {
                graph.add_node(query[q_idx])
            }
            (None, None) => continue,
        };

        if let Some(p) = prev_node {
            graph.add_edge(p, cur_node, seq_id);
        }
        prev_node = Some(cur_node);
    }

    graph.num_sequences += 1;
    graph.update_topological_order();
}
