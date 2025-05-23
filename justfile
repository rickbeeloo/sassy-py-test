bench bench='' *args='':
    cargo criterion --offline --plotting-backend disabled --bench bench -- "{{bench}}" {{args}}

build:
    cargo build -r --bench bench

stat bench='' *args='': build
    perf stat cargo bench --bench bench -- --profile-time 5 "{{bench}}" {{args}}

flame bench='' *args='':
    cargo flamegraph --release --bench bench -- --profile-time 5 "{{bench}}" 

record bench='' *args='': build
    perf record -g cargo bench --bench bench -- --profile-time 2 "{{bench}}" {{args}}
    perf report -n

cache bench='' *args='': build
    perf record -e cache-misses -g cargo bench --bench bench -- --profile-time 2 "{{bench}}" {{args}}
    perf report -n
    
report:
    perf report -n

cpufreq:
    sudo cpupower frequency-set --governor powersave -d 2.6GHz -u 2.6GHz
cpufreq-high:
    sudo cpupower frequency-set --governor powersave -d 5.0GHz -u 5.0GHz

heaptrack bench='' *args='': build
    heaptrack cargo bench --bench bench -- --profile-time 2 "{{bench}}" {{args}}
