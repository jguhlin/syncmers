# Syncmers Library in Rust

Syncmers as defined by Dutta et al. 2022, https://www.biorxiv.org/content/10.1101/2022.01.10.475696v2.full
Esp Fig 1b / Algorithm 1. Planning to implement other methods soon.

## Definition
Using the parameterized syncmer scheme, a syncmer is a kmer whose smallest smer is at a given target position (t).

## Extract Syncmers from &[u8]
```rust
let sequence = b"CCAGTGTTTACGG";
let syncmers = find_syncmers(5, 2, &[2], None, sequence);
assert!(syncmers == vec![b"CCAGT", b"TTACG"]);
println!("{:?}", syncmers);
```

## Extract Syncmers from &[u8], downsampling to 20%
```rust
let sequence = b"CCAGTGTTTACGG";
let syncmers = find_syncmers(5, 2, &[2], Some(0.2), sequence);
assert!(syncmers == vec![b"CCAGT", b"TTACG"]);
println!("{:?}", syncmers);
```

## Extract Syncmers from &[u8], keeping 80%
```rust
let sequence = b"CCAGTGTTTACGG";
let syncmers = find_syncmers(5, 2, &[2], Some(0.8), sequence);
assert!(syncmers == vec![b"CCAGT", b"TTACG"]);
println!("{:?}", syncmers);
```

## Find positions of Syncmers
```rust
let sequence = b"CCAGTGTTTACGG";
let syncmer_positions = find_syncmers_pos(5, 2, &[2], None, sequence);
println!("{:?}", syncmer_positions);
assert!(syncmer_positions == vec![0, 7]);
```

# TODO
Make sure X's are never the start / end of syncmers

# Changelog
0.1.4: Added downsampling support