# Sassy C bindings

First clone the repo and build using `cargo build --release --features c`. 
Then create your C file, see [example.c](example.c):


### Compile executable
Make sure `-L` points to the `target/release/` where the `.so` file is.
```bash
gcc -std=c11 -I. example.c \
    -L ../target/release -lsassy -lm \
    -o sassy_example
```

### Running executable
```bash
LD_LIBRARY_PATH=../target/release ./sassy_example
```

### Setup
The rust source for the library is in `src/c.rs` when the `c` feature is
enabled. `cbindgen` can be used to generate the corresponding [`sassy.h`](sassy.h) header:

```bash
cbindgen --config cbindgen.toml --crate sassy --output c/sassy.h
```
