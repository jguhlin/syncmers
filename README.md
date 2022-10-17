# Syncmers Library in Rust

Syncmers as defined by Dutta et al. 2022, https://www.biorxiv.org/content/10.1101/2022.01.10.475696v2.full
Esp Fig 1b / Algorithm 1. Planning to implement other methods soon

## Extract Syncmers from &[u8]
```rust
let sequence = b"CCAGTGTTTACGG";
let syncmers = find_syncmers(5, 2, &[2], sequence);
assert!(syncmers == vec![b"CCAGT", b"TTACG"]);
println!("{:?}", syncmers);
```

## Find positions of Syncmers
```rust
let sequence = b"CCAGTGTTTACGG";
let syncmer_positions = find_syncmers_pos(5, 2, &[2], sequence);
println!("{:?}", syncmer_positions);
assert!(syncmer_positions == vec![0, 7]);
```