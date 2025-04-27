
bench bench='' *args='':
    cargo criterion --offline --plotting-backend disabled --bench bench -- "{{bench}}" {{args}}

build:
    cargo build -r --bench bench

stat bench='' *args='': build
    perf stat cargo bench --bench bench -- --profile-time 5 "{{bench}}" {{args}}

record bench='' *args='': build
    perf record cargo bench --bench bench -- --profile-time 2 "{{bench}}" {{args}}
    perf report -n

report:
    perf report -n

cpufreq:
    sudo cpupower frequency-set --governor powersave -d 0.8GHz -u 2.6GHz
