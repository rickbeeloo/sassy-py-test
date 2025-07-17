# Sassy C bindings

First clone the repo and build using `cargo build --release --features c`. 
Then create your C file, see [example.c](example.c):


### Compile executable
Make sure `-L` points to the `target/release/` where the `.so` file is.
```bash
gcc -std=c11 -Wall -Wextra -I. example.c \
    -L ../target/release -lsassy -lm \
    -o sassy_example
```

### Running executable
```bash
LD_LIBRARY_PATH=../target/release ./sassy_example
```