from enola.router.router_mis import route_qubit_mis
from enola.router.codegen import CodeGen
import time

# from memory_profiler import profile

# @profile
def route_qubit(n_x: int, n_y: int, n_q: int, list_full_gates: list, qubit_mapping: list, routing_strategy: str, \
                reverse_to_initial: bool, l2: bool, use_window: bool):
    """
    generate rearrangement layers between two Rydberg layers
    """
    program_list = []
    final_mapping = qubit_mapping
    time_mis = 0
    time_codeGen = 0
    time_placement = 0
    for index_list_gate in range(len(list_full_gates)):
        # extract sets of movement that can be perform simultaneously
        data = []
        # t_sc = t_s
        t_s = time.time()
        data, final_mapping, time_placement_tmp = route_qubit_mis((n_x, n_y), n_q, index_list_gate, list_full_gates, list(final_mapping), routing_strategy, reverse_to_initial, l2, use_window)
        time_mis += (time.time() - t_s - time_placement_tmp)
        time_placement += time_placement_tmp
        data['n_x'] = n_x
        data['n_y'] = n_y
        data['n_r'] = n_y
        data['n_c'] = n_x
        # print("#layers: {}".format(len(data["layers"])))
        t_s = time.time()
        codegen = CodeGen(data)
        program = codegen.builder(no_transfer=False)
        tmp = program.emit_full()
        if index_list_gate == 0:
            program_list += tmp
        else:
            program_list += tmp[2:]
        time_codeGen += (time.time() - t_s)
        print("[INFO] Enola: Solve for Rydberg stage {}/{}.".format(index_list_gate+2, len(list_full_gates)))
    return program_list, time_mis, time_codeGen, time_placement
        
        # TODO: for each set of movement, collect the qubits that can be transfered simultaneously