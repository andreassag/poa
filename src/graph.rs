#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Edge {
    pub target: u32,
    pub weight: u32,
    pub seq_ids: Vec<u32>,
}

#[derive(Clone, Debug, Default)]
pub struct PoaGraph {
    pub bases: Vec<u8>,
    pub in_edges: Vec<Vec<Edge>>,
    pub out_edges: Vec<Vec<Edge>>,
    pub aligned_sibling: Vec<u32>, // Circular ring of alternative alleles at same column
    pub topo_order: Vec<u32>,
    pub num_sequences: u32,
}

impl PoaGraph {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn num_nodes(&self) -> usize {
        self.bases.len()
    }

    pub fn add_node(&mut self, base: u8) -> u32 {
        let id = self.bases.len() as u32;
        self.bases.push(base);
        self.in_edges.push(Vec::new());
        self.out_edges.push(Vec::new());
        self.aligned_sibling.push(id); // Self-loop ring
        id
    }

    pub fn add_edge(&mut self, from: u32, to: u32, seq_id: u32) {
        let f_idx = from as usize;
        let t_idx = to as usize;

        // Check if out-edge already exists
        let mut found = false;
        for e in &mut self.out_edges[f_idx] {
            if e.target == to {
                e.weight += 1;
                e.seq_ids.push(seq_id);
                found = true;
                break;
            }
        }
        if !found {
            self.out_edges[f_idx].push(Edge {
                target: to,
                weight: 1,
                seq_ids: vec![seq_id],
            });
        }

        // Check if in-edge already exists
        let mut in_found = false;
        for e in &mut self.in_edges[t_idx] {
            if e.target == from {
                e.weight += 1;
                e.seq_ids.push(seq_id);
                in_found = true;
                break;
            }
        }
        if !in_found {
            self.in_edges[t_idx].push(Edge {
                target: from,
                weight: 1,
                seq_ids: vec![seq_id],
            });
        }
    }

    // Link a new node into the circular aligned_sibling ring of an existing node
    pub fn link_sibling(&mut self, existing: u32, new_node: u32) {
        let next = self.aligned_sibling[existing as usize];
        self.aligned_sibling[existing as usize] = new_node;
        self.aligned_sibling[new_node as usize] = next;
    }

    // Recompute topological sorting via Kahn's algorithm
    pub fn update_topological_order(&mut self) {
        let n = self.num_nodes();
        let mut in_degree = vec![0; n];
        for u in 0..n {
            for e in &self.out_edges[u] {
                in_degree[e.target as usize] += 1;
            }
        }

        let mut queue = Vec::new();
        for i in 0..n {
            if in_degree[i] == 0 {
                queue.push(i as u32);
            }
        }

        let mut order = Vec::with_capacity(n);
        let mut head = 0;
        while head < queue.len() {
            let u = queue[head];
            head += 1;
            order.push(u);

            for e in &self.out_edges[u as usize] {
                let v = e.target as usize;
                in_degree[v] -= 1;
                if in_degree[v] == 0 {
                    queue.push(v as u32);
                }
            }
        }

        // If cycle detected or nodes disconnected, append remaining
        if order.len() < n {
            for i in 0..n {
                if !order.contains(&(i as u32)) {
                    order.push(i as u32);
                }
            }
        }

        self.topo_order = order;
    }

    // Initialize graph with first sequence
    pub fn init_with_sequence(&mut self, seq: &[u8]) {
        self.bases.clear();
        self.in_edges.clear();
        self.out_edges.clear();
        self.aligned_sibling.clear();
        self.topo_order.clear();

        if seq.is_empty() {
            return;
        }

        let mut prev = self.add_node(seq[0]);
        for &b in &seq[1..] {
            let cur = self.add_node(b);
            self.add_edge(prev, cur, 0);
            prev = cur;
        }

        self.num_sequences = 1;
        self.update_topological_order();
    }
}
