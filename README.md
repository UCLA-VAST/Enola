# Enola Compiler
Compilation for Dynamically Field-Programmable Qubit Arrays with Efficient and Provably Near-Optimal Scheduling.
Open source under the BSD 3-Clause license.

## Logistics
- We recommend to run the compiler in a Python3 virtual environment.
- The file `enola/router/codegen.py` is based on [a file from OLSQ-DPQA](https://github.com/UCLA-VAST/DPQA/blob/main/animation.py), so you need to install `networkx` and `matplotlib`.
- To compile circuits in the QASM format, you need `qiskit` installed. We used the API in `qiskit-1.1.0` specifically.
- The code for Misra&Gries coloring algorithm in `enola/scheduler/` is from this [Misra-Gries-coloring](https://codeberg.org/alifara/Misra-Gries-coloring).

## Repo structure
- `run.py` is an example of using the compiler on a circuit specified by a qubit interaction graph. Refer to `python run.py -h` for options.
- `run_qasm.py` is an example of using the compiler for circuits in the QASM format. Refer to `python run_qasm.py -h` for options.
- `enola/` contains the source files implementing Enola.
- `graphs.json` contains all the random 3-regular graphs we used as qubit interaction graphs.
- `animation.py` contains the class `Animator` that generates animations from the full code produced by Enola based on the implementation from [Animation.py in OLSQ-DPQA](https://github.com/UCLA-VAST/DPQA/blob/main/animation.py). Refer to `python animation.py -h` for options. 
- `simulator.py` contains the class `Simulator` that calculate the circuit fidelity based on Enola instructions. Refer to `python simulator.py -h` for options. 
- `qasm_exp/` is the directory containing one generic circuit example from QASMBench.
- `results/` is the default directory for the results.
  - `results/code/` contains the code files generated from compilation results.
  - `results/animations/` contains animation generated from full code files.
  - `results/fidelity/` contains fidelity estimation based on code files.

## How to use
- Run `python run.py <S> <I>` where `<S>` is the number of vertices in the random 3-regular graph, `<I>` is the id of the graph. To try other graphs, please edit `run.py` as needed.
  - Default setting: simulated annealing with dynamic placement, maximal independent set with sorting heuristic for routing strategy, no restriction for vectex number in solving MIS, and no animation code generation.    
  - The most scalable setting: trivial layout (`--trivial_layout`), return to initial mapping after each Rydberg stage (`--r2i`), using maximal independent set for routing (`--routing_strategy=maximalis`), and limit the number of vertices to 1000 in solving MIS (`--window`). \
  For example, `python run.py 30 0 --trivial_layout --r2i --routing_strategy=maximalis --window`
- (Optional) To generate animation, run `python animation.py <F>` where `<F>` is the full code file, e.g., `results/code/rand3reg_30_0_code_full.json`.