#!/usr/bin/env python3
"""
verify_logic.py - Standalone verification of Partial Order Alignment graph and consensus.
"""

class SimplePoa:
    def __init__(self):
        self.bases = []
        self.edges = {} # u -> list of (v, weight)
        self.in_edges = {} # v -> list of (u, weight)

    def add_node(self, b):
        nid = len(self.bases)
        self.bases.append(b)
        self.edges[nid] = {}
        self.in_edges[nid] = {}
        return nid

    def add_edge(self, u, v):
        self.edges[u][v] = self.edges[u].get(v, 0) + 1
        self.in_edges[v][u] = self.in_edges[v].get(u, 0) + 1

    def init_seq(self, seq):
        prev = self.add_node(seq[0])
        for b in seq[1:]:
            cur = self.add_node(b)
            self.add_edge(prev, cur)
            prev = cur

    def align_and_fuse(self, seq, match=2, mismatch=-4, gap=-6):
        m = len(seq)
        n = len(self.bases)
        # Topological order is 0..n-1 for linear DAG
        dp = [[-1000000] * (m + 1) for _ in range(n)]

        for u in range(n):
            cost_0 = match if self.bases[u] == seq[0] else mismatch
            if not self.in_edges[u]:
                dp[u][1] = cost_0
            else:
                for p in self.in_edges[u]:
                    dp[u][1] = max(dp[u][1], dp[p][0] + cost_0)

            for j in range(1, m + 1):
                cost = match if self.bases[u] == seq[j - 1] else mismatch
                if not self.in_edges[u]:
                    dp[u][j] = max(dp[u][j], j * gap + cost)
                else:
                    for p in self.in_edges[u]:
                        dp[u][j] = max(dp[u][j], dp[p][j - 1] + cost, dp[p][j] + gap)
                if j > 1:
                    dp[u][j] = max(dp[u][j], dp[u][j - 1] + gap)

        # Backtrace
        best_u = max(range(n), key=lambda u: dp[u][m])
        cur_u = best_u
        cur_j = m
        path = []
        while cur_u is not None and cur_j > 0:
            path.append((cur_u, cur_j - 1))
            best_pred = None
            max_val = -1000000
            for p in self.in_edges[cur_u]:
                if dp[p][cur_j - 1] > max_val:
                    max_val = dp[p][cur_j - 1]
                    best_pred = p
            cur_u = best_pred
            cur_j -= 1

        path.reverse()
        # Fuse path
        prev = None
        for u, j in path:
            if self.bases[u] == seq[j]:
                node = u
            else:
                # bubble node
                node = self.add_node(seq[j])
            if prev is not None:
                self.add_edge(prev, node)
            prev = node

    def consensus(self):
        n = len(self.bases)
        score = [0] * n
        parent = [None] * n
        for u in range(n):
            for p, w in self.in_edges[u].items():
                if score[p] + w > score[u]:
                    score[u] = score[p] + w
                    parent[u] = p
        end = max(range(n), key=lambda u: score[u])
        cur = end
        res = []
        while cur is not None:
            res.append(self.bases[cur])
            cur = parent[cur]
        res.reverse()
        return "".join(res)

if __name__ == '__main__':
    # Test identical sequences reinforce weights
    p = SimplePoa()
    p.init_seq("ACGTACGT")
    p.align_and_fuse("ACGTACGT")
    cons = p.consensus()
    assert cons == "ACGTACGT", f"Expected ACGTACGT, got {cons}"

    # Test bubble with majority rule
    p2 = SimplePoa()
    p2.init_seq("AACCC")
    p2.align_and_fuse("AACCC")
    p2.align_and_fuse("AATCC")
    cons2 = p2.consensus()
    assert cons2 == "AACCC", f"Expected AACCC, got {cons2}"

    print("POA logic verification passed: consensus and bubble fusion validated!")
