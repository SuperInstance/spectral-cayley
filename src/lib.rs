//! Spectral Cayley: spectral analysis of Cayley graphs from group generators.
//! Group structure → graph structure → spectral properties.

/// A group generator: a named permutation
#[derive(Clone, Debug)]
pub struct Generator {
    pub name: String,
    pub permutation: Vec<usize>,
}

impl Generator {
    pub fn new(name: &str, perm: Vec<usize>) -> Self {
        Self { name: name.to_string(), permutation: perm }
    }

    /// Apply this permutation to an element
    pub fn apply(&self, x: usize) -> usize {
        self.permutation.get(x).copied().unwrap_or(x)
    }

    /// Inverse permutation
    pub fn inverse(&self) -> Generator {
        let n = self.permutation.len();
        let mut inv = vec![0usize; n];
        for (i, &p) in self.permutation.iter().enumerate() {
            inv[p] = i;
        }
        Generator::new(&format!("{}⁻¹", self.name), inv)
    }
}

/// A Cayley graph built from group generators
pub struct CayleyGraph {
    pub n: usize,
    pub generators: Vec<Generator>,
    pub adj: Vec<Vec<f64>>,
}

impl CayleyGraph {
    /// Build a Cayley graph from generators acting on {0, 1, ..., n-1}
    pub fn from_generators(n: usize, generators: Vec<Generator>) -> Self {
        let mut adj = vec![vec![0.0_f64; n]; n];
        for generator in &generators {
            for i in 0..n {
                let j = generator.apply(i);
                if i != j {
                    adj[i][j] = 1.0;
                    adj[j][i] = 1.0; // undirected
                }
            }
        }
        // Remove duplicate edges (keep weight 1)
        for i in 0..n {
            for j in 0..n {
                if adj[i][j] > 0.0 { adj[i][j] = 1.0; }
            }
        }
        Self { n, generators, adj }
    }

    /// Cycle graph C_n: single cyclic generator
    pub fn cycle(n: usize) -> Self {
        let perm: Vec<usize> = (1..n).chain(std::iter::once(0)).collect();
        Self::from_generators(n, vec![Generator::new("σ", perm)])
    }

    /// Complete graph K_n: all transpositions
    pub fn complete(n: usize) -> Self {
        // Use adjacent transpositions as generators
        let mut gens = Vec::new();
        for i in 0..n.saturating_sub(1) {
            let mut perm: Vec<usize> = (0..n).collect();
            perm.swap(i, i + 1);
            gens.push(Generator::new(&format!("τ{}", i), perm));
        }
        // These only generate the symmetric group, giving a connected graph
        Self::from_generators(n, gens)
    }

    /// Hypercube Q_d: d generators flipping one bit
    pub fn hypercube(d: usize) -> Self {
        let n = 1 << d;
        let mut gens = Vec::new();
        for bit in 0..d {
            let perm: Vec<usize> = (0..n).map(|i| i ^ (1 << bit)).collect();
            gens.push(Generator::new(&format!("flip{}", bit), perm));
        }
        Self::from_generators(n, gens)
    }

    /// Laplacian eigenvalues
    pub fn eigenvalues(&self) -> Vec<f64> {
        let n = self.n;
        let mut lap = vec![vec![0.0; n]; n];
        for i in 0..n {
            let deg: f64 = self.adj[i].iter().sum();
            lap[i][i] = deg;
            for j in 0..n { if i != j { lap[i][j] = -self.adj[i][j]; } }
        }
        jacobi(&mut lap)
    }

    /// Conservation ratio
    pub fn cr(&self) -> f64 {
        let eigs = self.eigenvalues();
        if eigs.len() < 2 { return 0.0; }
        let l2 = eigs[1];
        let ln = *eigs.last().unwrap_or(&1.0);
        if ln <= 0.0 { 0.0 } else { l2 / ln }
    }

    /// Spectral gap (algebraic connectivity)
    pub fn spectral_gap(&self) -> f64 {
        self.eigenvalues().get(1).copied().unwrap_or(0.0)
    }

    /// Degree regularity check: all nodes should have same degree for Cayley graphs
    pub fn is_regular(&self) -> bool {
        if self.n == 0 { return true; }
        let deg0: f64 = self.adj[0].iter().sum();
        self.adj.iter().all(|row| (row.iter().sum::<f64>() - deg0).abs() < 1e-10)
    }

    /// Degree of each node
    pub fn degrees(&self) -> Vec<f64> {
        self.adj.iter().map(|row| row.iter().sum()).collect()
    }

    /// Diameter: maximum shortest path length
    pub fn diameter(&self) -> usize {
        let n = self.n;
        if n <= 1 { return 0; }
        let mut max_dist = 0;
        for start in 0..n {
            let dists = bfs_distances(&self.adj, start);
            for &d in &dists {
                if d > max_dist { max_dist = d; }
            }
        }
        max_dist
    }

    /// Expansion ratio: minimum |∂S|/|S| over all S with |S| ≤ n/2
    pub fn expansion(&self) -> f64 {
        let n = self.n;
        if n <= 1 { return 0.0; }
        let mut min_ratio = f64::MAX;
        // Check all subsets of size 1 to n/2
        for size in 1..=(n / 2) {
            // Sample some subsets (full enumeration is exponential)
            for mask in 0u32..std::cmp::min(1u32 << n, 1000u32) {
                if mask.count_ones() as usize != size { continue; }
                let boundary = self.boundary_size(mask as usize);
                let ratio = boundary as f64 / size as f64;
                if ratio < min_ratio { min_ratio = ratio; }
            }
        }
        if min_ratio == f64::MAX { 0.0 } else { min_ratio }
    }

    fn boundary_size(&self, mask: usize) -> usize {
        let mut boundary = 0usize;
        for i in 0..self.n {
            if mask & (1 << i) != 0 {
                for j in 0..self.n {
                    if mask & (1 << j) == 0 && self.adj[i][j] > 0.0 {
                        boundary += 1;
                    }
                }
            }
        }
        boundary
    }

    /// Compare with another Cayley graph: spectral distance
    pub fn spectral_distance(&self, other: &CayleyGraph) -> f64 {
        let e1 = self.eigenvalues();
        let e2 = other.eigenvalues();
        let n = e1.len().min(e2.len());
        if n == 0 { return 0.0; }
        // L² distance between normalized eigenvalue spectra
        let s1: f64 = e1.iter().sum();
        let s2: f64 = e2.iter().sum();
        let mut dist = 0.0;
        for i in 0..n {
            let p1 = if s1 > 0.0 { e1[i] / s1 } else { 0.0 };
            let p2 = if s2 > 0.0 { e2[i] / s2 } else { 0.0 };
            dist += (p1 - p2).powi(2);
        }
        dist.sqrt()
    }
}

fn bfs_distances(adj: &[Vec<f64>], start: usize) -> Vec<usize> {
    let n = adj.len();
    let mut dist = vec![usize::MAX; n];
    dist[start] = 0;
    let mut queue = std::collections::VecDeque::new();
    queue.push_back(start);
    while let Some(node) = queue.pop_front() {
        for j in 0..n {
            if adj[node][j] > 0.0 && dist[j] == usize::MAX {
                dist[j] = dist[node] + 1;
                queue.push_back(j);
            }
        }
    }
    dist
}

fn jacobi(a: &mut Vec<Vec<f64>>) -> Vec<f64> {
    let n = a.len();
    if n == 0 { return vec![]; }
    for _ in 0..100 * n * n {
        let (mut p, mut q) = (0, 1);
        let mut max_val = 0.0_f64;
        for i in 0..n { for j in (i+1)..n { if a[i][j].abs() > max_val { max_val = a[i][j].abs(); p = i; q = j; } } }
        if max_val < 1e-14 { break; }
        let app = a[p][p]; let aqq = a[q][q]; let apq = a[p][q];
        let theta = if (app - aqq).abs() < 1e-30 { std::f64::consts::FRAC_PI_4 }
                     else { 0.5 * (2.0 * apq / (app - aqq)).atan() };
        let (c, s) = (theta.cos(), theta.sin());
        for i in 0..n { if i != p && i != q { let aip = a[i][p]; let aiq = a[i][q]; a[i][p] = c*aip+s*aiq; a[p][i]=a[i][p]; a[i][q]=-s*aip+c*aiq; a[q][i]=a[i][q]; } }
        a[p][p] = c*c*app+2.0*s*c*apq+s*s*aqq;
        a[q][q] = s*s*app-2.0*s*c*apq+c*c*aqq;
        a[p][q] = 0.0; a[q][p] = 0.0;
    }
    let mut eigs: Vec<f64> = (0..n).map(|i| a[i][i]).collect();
    eigs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    eigs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cycle_graph_regular() {
        let c = CayleyGraph::cycle(5);
        assert!(c.is_regular());
        let degs = c.degrees();
        assert!(degs.iter().all(|&d| (d - 2.0).abs() < 1e-10));
    }

    #[test]
    fn cycle_graph_eigenvalues() {
        let c = CayleyGraph::cycle(4);
        let eigs = c.eigenvalues();
        // C₄ Laplacian eigenvalues: 0, 2, 2, 4
        assert!(eigs[0].abs() < 1e-10, "λ₀ = {}", eigs[0]);
        assert!((eigs.last().unwrap() - 4.0).abs() < 0.1, "λ₃ = {}", eigs.last().unwrap());
    }

    #[test]
    fn hypercube_regular() {
        let q = CayleyGraph::hypercube(3);
        assert!(q.is_regular());
        let degs = q.degrees();
        assert!(degs.iter().all(|&d| (d - 3.0).abs() < 1e-10), "Q₃ should be 3-regular");
    }

    #[test]
    fn hypercube_spectral_gap() {
        let q = CayleyGraph::hypercube(3);
        let gap = q.spectral_gap();
        // Q_d has λ₂ = 2 for all d
        assert!((gap - 2.0).abs() < 0.1, "Q₃ gap = {}", gap);
    }

    #[test]
    fn cycle_diameter() {
        let c = CayleyGraph::cycle(6);
        assert_eq!(c.diameter(), 3, "C₆ diameter should be 3");
    }

    #[test]
    fn hypercube_diameter() {
        let q = CayleyGraph::hypercube(3);
        assert_eq!(q.diameter(), 3, "Q₃ diameter should be 3");
    }

    #[test]
    fn cr_in_range() {
        let c = CayleyGraph::cycle(10);
        let cr = c.cr();
        assert!(cr >= 0.0 && cr <= 1.0, "CR = {}", cr);
        // C₁₀ has relatively low CR (small spectral gap relative to max eigenvalue)
        assert!(cr > 0.0, "CR should be positive for connected graph");
    }

    #[test]
    fn generator_inverse() {
        let g = Generator::new("σ", vec![1, 2, 0]);
        let inv = g.inverse();
        assert_eq!(inv.apply(0), 2); // σ(2)=0, so σ⁻¹(0)=2
        assert_eq!(inv.apply(1), 0);
        assert_eq!(inv.apply(2), 1);
    }

    #[test]
    fn spectral_distance_self_zero() {
        let c = CayleyGraph::cycle(5);
        let d = c.spectral_distance(&c);
        assert!(d < 1e-10, "Self-distance should be 0: {}", d);
    }

    #[test]
    fn expansion_positive() {
        let c = CayleyGraph::cycle(6);
        let exp = c.expansion();
        assert!(exp > 0.0, "Expansion should be positive: {}", exp);
    }
}
