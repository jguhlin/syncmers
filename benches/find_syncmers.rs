use criterion::{black_box, criterion_group, criterion_main, Criterion};
use pulp::Arch;

pub struct Syncmers {
    pub k: usize,
    pub s: usize,
    pub t: usize,
    // pub downsample: f32,
}

impl Syncmers {
    pub fn new(k: usize, s: usize, t: usize) -> Self {
        assert!(s < k);
        assert!(t < k - s);
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

    pub fn find_all_simdbyhand(&self, seq: &[u8]) -> Vec<usize> {
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

    pub fn find_all_pulp(&self, seq: &[u8]) -> Vec<usize> {
        let arch = Arch::new();
        arch.dispatch(|| {
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
        })
    }
}

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

    pub fn find_all_pulp(&self, seq: &[u8]) -> Vec<usize> {
        let arch = Arch::new();
        arch.dispatch(|| {
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
        })
    }
}

// TODO: These should be separate groups...
fn criterion_benchmark(c: &mut Criterion) {
    let sequence = b"ccaaattgaaagtgagtgtctaatgtattaattagtgaaataatatcttgatatttctttaagggagaattctgccaaattgaaagtgagtgtctaatgtattaattagtgaaataatatcttgatatttctttaagggagaattctgccaaattgaaagtgagtgtctaatgtattaattagtgaaataatatcttgatatttctttaagggagaattctgccaaattgaaagtgagtgtctaatgtattaattagtgaaataatatcttgatatttctttaagggagaattctgccaaattgaaagtgagtgtctaatgtattaattagtgaaataatatcttgatatttctttaagggagaattctgccaaattgaaagtgagtgtctaatgtattaattagtgaaataatatcttgatatttctttaagggagaattctg";
    let mut sequence = sequence.to_vec();
    sequence.make_ascii_uppercase();

    let ts: [usize; 1] = [2];

    let mut group = c.benchmark_group("syncmers");
    group.throughput(criterion::Throughput::Bytes(sequence.len() as u64));

    group.bench_function("syncmers", |b| {
        b.iter(|| {
            let syncmers = Syncmers::new(5, 2, 2);
            syncmers.find_all(black_box(&sequence));
        })
    });

    group.bench_function("syncmers_pulp", |b| {
        b.iter(|| {
            let syncmers = Syncmers::new(5, 2, 2);
            syncmers.find_all_pulp(black_box(&sequence));
        })
    });

    group.bench_function("psyncmers_1param", |b| {
        b.iter(|| {
            let syncmers = ParameterizedSyncmers::new(5, 2, &ts);
            syncmers.find_all(black_box(&sequence));
        })
    });

    group.bench_function("psyncmers_1param_pulp", |b| {
        b.iter(|| {
            let syncmers = ParameterizedSyncmers::new(5, 2, &ts);
            syncmers.find_all_pulp(black_box(&sequence));
        })
    });

    let ts: [usize; 2] = [2, 3];

    group.bench_function("psyncmers_2param", |b| {
        b.iter(|| {
            let syncmers = ParameterizedSyncmers::new(5, 2, &ts);
            syncmers.find_all(black_box(&sequence));
        })
    });

    group.bench_function("psyncmers_2param_pulp", |b| {
        b.iter(|| {
            let syncmers = ParameterizedSyncmers::new(5, 2, &ts);
            syncmers.find_all_pulp(black_box(&sequence));
        })
    });

    group.finish();
}

criterion_group! {
    name=syncmers;
    config = Criterion::default().significance_level(0.05).measurement_time(std::time::Duration::from_secs(10));
    targets=criterion_benchmark
}

criterion_main!(syncmers);
