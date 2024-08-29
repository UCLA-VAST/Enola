from enola.scheduler.gate_scheduler import gate_scheduling
from enola.placer.placer import place_qubit
from enola.router.router import route_qubit
from enola.router.codegen import global_dict
from typing import Sequence
import time
import json

class Enola:
    """class to solve QLS problem."""

    def __init__(self, name: str, dir: str = None, print_detail: bool = False, 
                 placement_opt = -1,
                 routing_opt = -1,
                 trivial_layout = False, routing_strategy = "maximalis_sorted",\
                to_verify: bool = False, reverse_to_initial:bool = False,\
                initial_mapping: list = None, dependency: bool = False, l2: bool = False, use_window: bool = False, full_code: bool = False):
        self.dir = dir
        self.n_q = 0 # number of qubits
        self.n_g = 0 # number of gates
        self.n_x = 0 # dim for arch
        self.n_y = 0 # dim for arch
        self.n_c = 0 # dim for arch
        self.n_r = 0 # dim for arch
        self.print_detail = print_detail
        self.all_commutable = False
        self.result_json = {}
        self.result_json['name'] = name
        self.result_json['layers'] = []
        self.to_verify = to_verify
        if placement_opt == 0:
            self.trivial_layout = True
            self.reverse_to_initial = True
        elif placement_opt == 1:
            self.trivial_layout = False
            self.self.reverse_to_initial = True
        elif placement_opt == 2:
            self.trivial_layout = False
            self.self.reverse_to_initial = False
        if routing_opt == 0:
            self.routing_strategy = "maximalis_sorted"
            self.use_window = True
        elif routing_opt == 1:
            self.routing_strategy = "maximalis_sorted"
            self.use_window = False
        elif routing_opt == 2:
            self.routing_strategy = "mis"
            self.use_window = False
        self.trivial_layout = trivial_layout
        self.reverse_to_initial = reverse_to_initial
        self.given_initial_mapping = initial_mapping
        self.has_dependency = dependency
        self.routing_strategy = routing_strategy
        self.use_window = use_window
        self.l2 = l2
        global_dict["full_code"] = full_code

    def setArchitecture(self, bounds: Sequence[int]):
        # bounds = [number of X, number of Y, number of C, number of R]
        self.n_x, self.n_y, self.n_c, self.n_r = bounds

    def setProgram(self, program: Sequence[Sequence[int]], nqubit: int = None):
        # assume program is a iterable of pairs of qubits in 2Q gate
        # assume that the qubit indices used are consecutively 0, 1, ...
        self.n_g = len(program)
        self.g_i = [i for i in range(self.n_g)]
        self.g_q = [(min(pair), max(pair)) for pair in program]
        self.g_s = tuple(['CRZ' for _ in range(self.n_g)])
        if not nqubit:
            for gate in program:
                self.n_q = max(gate[0], self.n_q)
                self.n_q = max(gate[1], self.n_q)
            self.n_q += 1
        else:
            self.n_q = nqubit

    def solve(self, save_file: bool = True):
        runtime_analysis = {}
        print("[INFO] Enola: Start Solving")
        if self.n_q > self.n_x * self.n_y:
            print("[Error] #qubits > #sites. There may be a problem.")
        # self.writeSettingJson()
        t_s = time.time()
        # gate scheduling with graph coloring
        print("[INFO] Enola: Run scheduling")
        if self.has_dependency:
            result_scheduling = self.asap()
        else:
            result_scheduling = gate_scheduling(self.n_q, self.g_q)
        if self.to_verify:
            self.verify_scheduling(result_scheduling)
        list_gates = []
        for gates in result_scheduling:
            tmp = [self.g_q[i] for i in gates]
            list_gates.append(tmp)
        runtime_analysis["scheduling"] = time.time()- t_s
        print("[INFO] Enola: Time for scheduling: {}s".format(runtime_analysis["scheduling"]))

        t_p = time.time()
        if self.given_initial_mapping is not None:
            qubit_mapping = self.given_initial_mapping
        else:
            # qubit placement for layout
            qubit_mapping = []
            if self.trivial_layout:
                length = self.n_x
                print(length)
                x = 0
                y = 0
                for i in range(self.n_q):
                    qubit_mapping.append((x, y))
                    x += 1
                    if x % length == 0:
                        x = 0
                        y += 1
                # return
            else:
                qubit_mapping = place_qubit((self.n_x, self.n_y), self.n_q, list_gates, self.l2)
        runtime_analysis["placement"] = time.time()- t_p
        print("[INFO] Enola: Time for placement: {}s".format(runtime_analysis["placement"]))
        if self.to_verify:
            self.verify_qubit_mapping(qubit_mapping)
        # qubit movement between layers
        program_list, time_mis, time_codeGen, time_placement = route_qubit(self.n_x, self.n_y, self.n_q, list_gates, qubit_mapping, self.routing_strategy, self.reverse_to_initial, self.l2, self.use_window)
        runtime_analysis["routing"] = time_mis
        runtime_analysis["codegen"] = time_codeGen
        runtime_analysis["placement"] = runtime_analysis["placement"] + time_placement
        runtime_analysis["total"] = time.time()- t_s
        print("[INFO] Enola: Time for routing: {}s".format(runtime_analysis["routing"]))
        print("[INFO] Enola: Toal Time: {}s".format(runtime_analysis["total"]))
        if save_file:
            if not self.dir:
                self.dir = "./results/"
            if global_dict["full_code"]:
                with open(self.dir + f"code/{self.result_json['name']}_code_full.json", 'w') as f:
                    json.dump(program_list, f)
                    for instruction in program_list:
                        instruction["state"] = {}
            with open(self.dir + f"code/{self.result_json['name']}_code.json", 'w') as f:
                json.dump(program_list, f)
            with open(self.dir + f"time/{self.result_json['name']}_time.json", 'w') as f:
                json.dump(runtime_analysis, f)
        return program_list
    
    def asap(self):
        # as soon as possible algorithm for self.g_q
        list_scheduling = []
        list_qubit_time = [0 for i in range(self.n_q)]
        for i, gate in enumerate(self.g_q):
            tq0 = list_qubit_time[gate[0]]
            tq1 = list_qubit_time[gate[1]]
            tg = max(tq0, tq1)
            if tg >= len(list_scheduling):
                list_scheduling.append([])
            list_scheduling[tg].append(i)

            tg += 1
            list_qubit_time[gate[0]] = tg
            list_qubit_time[gate[1]] = tg
        return list_scheduling
    def verify_scheduling(self, result_scheduling: list):
        time_gate_scheduled = [-1 for i in range(len(self.g_q))]
        success = True
        for i, stage in enumerate(result_scheduling):
            qubit_gate = [-1 for i in range(self.n_q)]
            for gate in stage:
                if time_gate_scheduled[gate] > -1:
                    success = False
                    print("[Error] Gate {} is already scheduled in stage {}, but is assigned to stage {} again.".format(gate, time_gate_scheduled[gate], i))
                time_gate_scheduled[gate] = i
                q0 = self.g_q[gate][0]
                q1 = self.g_q[gate][1]
                if qubit_gate[q0] > -1:
                    success = False
                    print("[Error] Qubit {} is already used in gate {}, but is involved in gate {} at the same stage ({}) again.".format(q0, qubit_gate[q0], gate, i))
                if qubit_gate[q1] > -1:
                    success = False
                    print("[Error] Qubit {} is already used in gate {}, but is involved in gate {} at the same stage ({}) again.".format(q1, qubit_gate[q1], gate, i))
        
        for i, t in enumerate(time_gate_scheduled):
            if t == -1:
                success = False
                print("[Error] Gate {} is not scheduled.".format(i))
        if success:
            print(f"[INFO] Gate Scheduling Verification: Pass.")

    def verify_qubit_mapping(self, result_mapping: list):
        success = True
        physical_qubit_idx = [[-1 for j in range(self.n_y)] for i in range(self.n_x)]
        for i, location in enumerate(result_mapping):
            if location[0] >= self.n_x or location[0] < 0:
                print("[Error] Qubit {} is mapped outside the chip ({},{}).".format(i, location[0], location[1]))
                success = False
            if location[1] >= self.n_y or location[1] < 0:
                print("[Error] Qubit {} is mapped outside the chip ({},{}).".format(i, location[0], location[1]))
                success = False
            if physical_qubit_idx[location[0]][location[1]] > -1:
                print("[Error] Qubit {} is overlapped with qubit {} at ({},{}).".format(i, physical_qubit_idx[location[0]][location[1]], location[0], location[1]))
            physical_qubit_idx[location[0]][location[1]] = i
        if len(result_mapping) != self.n_q:
            print("[Error] Not all qubits are mapped: len(result_mapping)= {}, #qubit={}.".format(len(result_mapping), self.n_q))
            success = False
        if success:
            print(f"[INFO] Qubit Placement Verification: Pass.")
    
    def verify_code(self, program_list: list):
        raise NotImplementedError
