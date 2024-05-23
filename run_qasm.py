from enola.enola import Enola
import qiskit.qasm2
import argparse
import qiskit.circuit
from qiskit import transpile

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('qasm', metavar='S', type=str, help='qasm file name')
    parser.add_argument('--arch', help='architecture dimension', type=int, default=16)
    parser.add_argument('--routing_strategy', help='routing strategy: mis, maximalis, maximalis_sort', type=str, default="maximalis")
    parser.add_argument('--trivial_layout', action='store_true', default=False)
    parser.add_argument('--r2i',  help='whether reverse to initial mapping after each Rydberg', action='store_true', default=False)
    parser.add_argument('--window',  help='restrict vertex number in mis to 1000', action='store_true', default=False)
    parser.add_argument('--full_code',  help='generate full code for animation. The size of the code file could be large.', action='store_true', default=False)
    args = parser.parse_args()

    list_gate_two_qubit = []
    with open(args.qasm, 'r') as f:
        qasm_str = f.read()
        circuit = qiskit.qasm2.loads(qasm_str)
        cz_circuit = transpile(circuit, basis_gates=['cz', 'rx', 'ry', 'rz', 'h', 't'])
        instruction = cz_circuit.data
        for ins in instruction:
            if ins.operation.num_qubits == 2:
                list_gate_two_qubit.append((ins.qubits[0]._index, ins.qubits[1]._index))
    print(list_gate_two_qubit)
    filename = args.qasm.split('/')[-1]
    filename = filename.split('.')[0]


    tmp = Enola(filename,
        dir='./results/',
        trivial_layout = args.trivial_layout,
        routing_strategy=args.routing_strategy,
        reverse_to_initial=args.r2i,
        dependency = True,
        use_window=args.window,
        full_code=args.full_code
    )

    tmp.setArchitecture([args.arch, args.arch, args.arch, args.arch])
    tmp.setProgram(list_gate_two_qubit)
    tmp.solve(save_file=True)



