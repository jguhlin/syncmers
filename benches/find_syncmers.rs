use criterion::{black_box, criterion_group, criterion_main, Criterion};
use pulp::Arch;

pub fn find_syncmers(k: usize, s: usize, t: &[usize], seq: &[u8]) -> Vec<usize> {
    assert!(seq.len() > k);
    assert!(s < k);
    assert!(t.iter().all(|&t| t < k));
    assert!(t.len() < 5, "Only supports up to 4 syncmers. Email if you'd like more (or change the lines)");

    let mut ts: [usize; 4] = [0; 4];
    let t_len = t.len();
    ts[..t.len()].copy_from_slice(t);
    let ts = &ts[..t_len];

   seq.windows(k)
        .enumerate()
        .filter_map(|(i, kmer)| {
            let min_pos = kmer
                .windows(s)
                .enumerate()
                .min_by(|(_, a), (_, b)| a.cmp(b));

            if (t_len == 1 && min_pos.unwrap().0 == ts[0]) || ts.contains(&min_pos.unwrap().0) {
                Some(i)
            } else {
                None
            }
        })
        .collect::<Vec<_>>()
}

pub fn find_syncmers_anonfn(k: usize, s: usize, ts: &[usize], seq: &[u8]) -> Vec<usize> {
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

pub fn find_syncmers_anonfn_pulp(k: usize, s: usize, ts: &[usize], seq: &[u8]) -> Vec<usize> {
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

            let arch = Arch::new();
            if arch.dispatch(|| cmp(ts, min_pos.unwrap().0)) {
                Some(i)
            } else {
                None
            }
           
        })
        .collect::<Vec<_>>()
}

pub fn find_syncmers_anonfn_pulp_u8(k: usize, s: usize, ts: &[u8], seq: &[u8]) -> Vec<usize> {
    assert!(seq.len() > k);
    assert!(s < k);
    assert!(ts.iter().all(|&t| t as usize <= k as usize - s as usize));
    assert!(ts.len() < 5, "Only supports up to 4 syncmers. Email if you'd like more (or change the lines)");

    let cmp = match ts.len() {
        1 => |x: &[u8], y: usize| -> bool { x[0] == y as u8 },
        2 => |x: &[u8], y: usize| -> bool { x[0] == y as u8 || x[1] == y as u8 },
        3 => |x: &[u8], y: usize| -> bool { x[0] == y  as u8|| x[1] == y  as u8|| x[2] == y  as u8},
        4 => |x: &[u8], y: usize| -> bool { x[0] == y  as u8|| x[1] == y  as u8|| x[2] == y  as u8|| x[3] == y  as u8},
        _ => unreachable!(),
    };

   seq.windows(k)
        .enumerate()
        .filter_map(|(i, kmer)| {
            let min_pos = kmer
                .windows(s)
                .enumerate()
                .min_by(|(_, a), (_, b)| a.cmp(b));

            let arch = Arch::new();
            if arch.dispatch(|| cmp(ts, min_pos.unwrap().0)) {
                Some(i)
            } else {
                None
            }
           
        })
        .collect::<Vec<_>>()
}

pub fn find_syncmers_anonfn_u8(k: usize, s: usize, ts: &[u8], seq: &[u8]) -> Vec<usize> {
    assert!(seq.len() > k);
    assert!(s < k);
    assert!(ts.iter().all(|&t| t as usize <= k as usize - s as usize));
    assert!(ts.len() < 5, "Only supports up to 4 syncmers. Email if you'd like more (or change the lines)");
    assert!(!ts.is_empty());

    let cmp = match ts.len() {
        1 => |x: &[u8], y: usize| -> bool { x[0] == y as u8 },
        2 => |x: &[u8], y: usize| -> bool { x[0] == y as u8 || x[1] == y as u8 },
        3 => |x: &[u8], y: usize| -> bool { x[0] == y as u8 || x[1] == y as u8 || x[2] == y as u8 },
        4 => |x: &[u8], y: usize| -> bool { x[0] == y as u8 || x[1] == y as u8 || x[2] == y as u8 || x[3] == y as u8},
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

pub fn find_syncmers_anonfn_u8_contains(k: usize, s: usize, ts: &[u8], seq: &[u8]) -> Vec<usize> {
    assert!(seq.len() > k);
    assert!(s < k);
    assert!(ts.iter().all(|&t| t as usize <= k as usize - s as usize));
    assert!(ts.len() < 5, "Only supports up to 4 syncmers. Email if you'd like more (or change the lines)");

    let cmp = match ts.len() {
        1 => |x: &[u8], y: usize| -> bool { x[0] == y as u8 },
        2 => |x: &[u8], y: usize| -> bool { x.contains(&(y as u8)) },
        3 => |x: &[u8], y: usize| -> bool { x.contains(&(y as u8)) },
        4 => |x: &[u8], y: usize| -> bool { x.contains(&(y as u8)) },
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

// TODO: These should be separate groups...
fn criterion_benchmark(c: &mut Criterion) {
    let sequence = b"ccaaattgaaagtgagtgtctaatgtattaattagtgaaataatatcttgatatttctttaagggagaattctgccaaattgaaagtgagtgtctaatgtattaattagtgaaataatatcttgatatttctttaagggagaattctgccaaattgaaagtgagtgtctaatgtattaattagtgaaataatatcttgatatttctttaagggagaattctgccaaattgaaagtgagtgtctaatgtattaattagtgaaataatatcttgatatttctttaagggagaattctgccaaattgaaagtgagtgtctaatgtattaattagtgaaataatatcttgatatttctttaagggagaattctgccaaattgaaagtgagtgtctaatgtattaattagtgaaataatatcttgatatttctttaagggagaattctg";
    let mut sequence = sequence.to_vec();
    sequence.make_ascii_uppercase();

    let ts: [usize; 1] = [2];

    let mut group = c.benchmark_group("syncmers");
    group.throughput(criterion::Throughput::Bytes(sequence.len() as u64));

    // 270 MiB/s
    group.bench_function("syncmers_1param_anonfn", |b| {
        b.iter(|| {
            find_syncmers_anonfn(5, 2, &[2], black_box(&sequence));
        })
    });

    // 244 MiB/s
    group.bench_function("syncmers_2param_anonfn", |b| {
        b.iter(|| {
            find_syncmers_anonfn(5, 2, &[2, 3], black_box(&sequence));
        })
    });

    // 298 MiB/s
    group.bench_function("syncmers_anonfn_u8", |b| {
        b.iter(|| {
            find_syncmers_anonfn_u8(5, 2, &[2], black_box(&sequence));
        })
    });

    // 264 MiB/s
    group.bench_function("syncmers_2param_anonfn_u8", |b| {
        b.iter(|| {
            find_syncmers_anonfn_u8(5, 2, &[2, 3], black_box(&sequence));
        })
    });    

    group.bench_function("find_syncmers_anonfn_u8_contains", |b| {
        b.iter(|| find_syncmers_anonfn_u8_contains(31, 5, &[2], &sequence))
    });

    group.bench_function("find_syncmers_2param_anonfn_u8_contains", |b| {
        b.iter(|| find_syncmers_anonfn_u8_contains(31, 5, &[2, 3], &sequence))
    });

    /* 97 MiB/s
    group.bench_function("syncmers", |b| {
        b.iter(|| {
            find_syncmers(5, 2, &[2], black_box(&sequence));
        })
    }); */

    /* 280 MiB/s
    group.bench_function("syncmers_anonfn", |b| {
        b.iter(|| {
            find_syncmers_anonfn(5, 2, &[2], black_box(&sequence));
        })
    });
    */

    /* 84 MiB/s
    group.bench_function("syncmers_anonfn_pulp", |b| {
        b.iter(|| {
            find_syncmers_anonfn_pulp(5, 2, &[2], black_box(&sequence));
        })
    });
    */

    /* 86 MiB/s
    group.bench_function("syncmers_anonfn_pulp_u8", |b| {
        b.iter(|| {
            find_syncmers_anonfn_pulp_u8(5, 2, &[2], black_box(&sequence));
        })
    });
    */


    /* 95 MiB/s
    group.bench_function("psyncmers_2param", |b| {
        b.iter(|| {
            find_syncmers(5, 2, &[2, 3], black_box(&sequence));
        })
    });
    */

    

    // 78 MiB/s
    group.bench_function("psyncmers_2param_anonfn_pulp", |b| {
        b.iter(|| {
            find_syncmers_anonfn_pulp(5, 2, &[2, 3], black_box(&sequence));
        })
    });

    // 78 MiB/s
    /*
    group.bench_function("psyncmers_anonfn_pulp_u8", |b| {
        b.iter(|| {
            find_syncmers_anonfn_pulp_u8(5, 2, &[2, 3], black_box(&sequence));
        })
    }); */

    // 78 MiB/s
    /*group.bench_function("psyncmers_anonfn_u8", |b| {
        b.iter(|| {
            find_syncmers_anonfn_pulp_u8(5, 2, &[2, 3], black_box(&sequence));
        })
    });*/

    group.finish();
}

criterion_group! {
    name=syncmers;
    config = Criterion::default().significance_level(0.05).measurement_time(std::time::Duration::from_secs(10));
    targets=criterion_benchmark
}

criterion_main!(syncmers);
