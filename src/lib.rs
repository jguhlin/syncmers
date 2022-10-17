/// Syncmers as defined by Dutta et al. 2022, https://www.biorxiv.org/content/10.1101/2022.01.10.475696v2.full
/// Esp Fig 1b
/// Planning to implement other methods soon
///
/// TODO: Add Iterator impl's

// use std::iter::{FilterMap, Enumerate};
// use std::slice::Windows;
use std::cmp::Ordering;

use pulp::Arch;

// TODO:Denote the reverse complement of x by Embedded Image. For a given order, the canonical form of a k-mer x, denoted by Canonical(x), is the smaller of x and Embedded Image. For example, under the lexicographic order, Canonical(CGGT) = ACCG.
// Canonical(x) = min(x, revcomp(x))

// Copied from ffforf.
pub fn complement(c: &mut u8) {
    if *c != b'N' {
        if *c & 2 != 0 {
            *c ^= 4;
        } else {
            *c ^= 21;
        }
    }
}

pub fn revcomp(sequence: &mut [u8]) {
    let arch = Arch::new();
    arch.dispatch(|| {
        sequence.reverse();
        sequence.make_ascii_uppercase();
        sequence.iter_mut().for_each(complement);
    });
}

pub fn is_revcomp_min(seq: &[u8]) -> bool {
    assert!(!seq.is_empty());
    for i in 0..seq.len() {
        let mut c = seq[seq.len() - i - 1];
        complement(&mut c);
        match seq[i].cmp(&c) {
            Ordering::Less => return true,
            Ordering::Greater => return false,
            Ordering::Equal => continue,
        }
    }

    false
}

// Best as determined by criterion benchmarks
pub fn find_syncmers(k: usize, s: usize, ts: &[usize], seq: &[u8]) -> Vec<usize> {
    assert!(seq.len() > k);
    assert!(s < k);
    assert!(ts.iter().all(|&t| t <= k - s));
    assert!(ts.len() < 5, "Only supports up to 4 syncmers. Email if you'd like more (or change the lines)");

    let cmp = match ts.len() {
        1 => |x: &[usize], y: usize| -> bool { x[0] == y },
        2 => |x: &[usize], y: usize| -> bool { x[0] == y || x[1] == y },
        3 => |x: &[usize], y: usize| -> bool { x[0] == y || x[1] == y || x[2] == y },
        4 => |x: &[usize], y: usize| -> bool { x[0] == y || x[1] == y || x[2] == y || x[3] == y },
        _ => unreachable!(),
    };

   seq.windows(k)
        .enumerate()
        .filter_map(|(i, kmer)| {
            let min_pos = kmer
                .windows(s)
                .enumerate()
                .min_by(|(_, a), (_, b)| a.cmp(b));

            if cmp(ts, min_pos.unwrap().0) {
                Some(i)
            } else {
                None
            }
           
        })
        .collect::<Vec<_>>()
}

// NOTE: "By convention, ties are broken by choosing the leftmost position"

/// 1-parameter syncmer method
/// t is 0-based (unlike in the paper)
/// NOTE: Sequence should be all upper case (or all lower case)
// TODO: Remove?
pub struct Syncmers {
    pub k: usize,
    pub s: usize,
    pub t: usize,
}

// type FilterMapIter<'a> = FilterMap<Enumerate<Windows<'a, u8>>, &'static fn ((usize, &'a [u8])) -> Option<usize>>;

/* Docs for getting the canonical strand.
let mut revcmp: Vec<u8>;
let mut rev = false;

// Get the canonical strand
let seq = if is_revcomp_min(seq) {
    rev = true;

    revcmp = seq.to_vec();
    revcomp(&mut revcmp);
    &revcmp
} else {
    seq
};

*/
impl Syncmers {
    pub fn new(k: usize, s: usize, t: usize) -> Self {

        assert!(s < k);
        assert!(t < k);
        Syncmers { k, s, t }
    }

    // TODO: Find a way to return just the iter (FilterMap Iter and it's long return type)

    pub fn find_all(&self, seq: &[u8]) -> Vec<usize> {
        assert!(seq.len() >= self.k);

       seq.windows(self.k)
            .enumerate()
            .filter_map(|(i, kmer)| {
                let min_pos = kmer
                    .windows(self.s)
                    .enumerate()
                    .min_by(|(_, a), (_, b)| a.cmp(b));

                if min_pos.unwrap().0 == self.t {
                    Some(i)
                } else {
                    None
                }
            })
            .collect::<Vec<_>>()
    }

    /*
    type FilterMapIter<'a> = FilterMap<Enumerate<Windows<'a, u8>>, &'static fn ((usize, &'a [u8])) -> Option<usize>>;

    pub fn find<'a>(&self, seq: &'a [u8]) -> FilterMapIter<'a> {
        assert!(seq.len() >= self.k);
        seq.windows(self.k)
            .enumerate()
            .filter_map(|(i, kmer)| {
                if let Some(min_pos) = kmer
                    .windows(self.s)
                    .enumerate()
                    .min_by(|(_, a), (_, b)| a.cmp(b)) {

                        if min_pos == self.t {
                            Some(i)
                        } else {
                            None        
                        }
                    } else {
                        None
                    }
            }) 

    } */
}

pub struct SyncmersIter<'a> {
    pub sequence: &'a [u8],
    pub k: usize,
    pub s: usize,
    pub t: usize,
}


impl<'a> Iterator for SyncmersIter<'a> {
    type Item = &'a [u8];

    fn next(&mut self) -> Option<Self::Item> {
        if self.sequence.len() < self.k {
            return None;
        }

        let kmer = &self.sequence[..self.k];
        self.sequence = &self.sequence[1..];

        Some(kmer)
    }
}

/// Multi-parameter syncmer method
/// t is 0-based (unlike in the paper)
pub struct ParameterizedSyncmers<'a> {
    pub k: usize,
    pub s: usize,
    pub t: &'a [usize],
}

impl<'a> ParameterizedSyncmers<'a> {
    pub fn new(k: usize, s: usize, t: &'a [usize]) -> Self {
        assert!(s < k);
        assert!(t.iter().all(|&t| t < k));
        ParameterizedSyncmers { k, s, t }
    }

    pub fn find_all(&self, seq: &[u8]) -> Vec<usize> {
        seq.windows(self.k)
            .enumerate()
            .filter_map(|(i, kmer)| {
                let min_pos = kmer
                    .windows(self.s)
                    .enumerate()
                    .min_by(|(_, a), (_, b)| a.cmp(b));

                if self.t.contains(&min_pos.unwrap().0) {
                    Some(i)
                } else {
                    None
                }
            })
            .collect::<Vec<_>>()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    pub fn test_syncmers_fig1b() {
        let sequence = b"CCAGTGTTTACGG";
        let syncmer_positions = find_syncmers(5, 2, &[2], sequence);
        println!("{:?}", syncmer_positions);
        assert!(syncmer_positions == vec![0, 7]);

        let ts: [usize; 1] = [2];

        let psyncmers = ParameterizedSyncmers::new(5, 2, &ts);
        let syncmer_positions = psyncmers.find_all(sequence);
        assert!(syncmer_positions == vec![0, 7]);
    }
}
