import sys
import math
from random import randrange, shuffle, uniform, seed

seed(0)

class SAPlacerPartial:
    """class to find a qubit layout via SA."""

    def __init__(self, l2: bool = False):
        self.initialize_param()
        self.n_qubit = 0
        self.chip_dim = (0, 0)
        self.list_gate = []
        self.movement = []
        self.list_qubit_list_gate = []
        self.list_position = []
        self.initial_mapping = []
        self.target_qubits = []
        self.l2 = l2

    def initialize_param(self):
        # sa related parameter    
        self.sa_t = 100000.0
        self.sa_t1 = 4.0
        self.sa_t_frozen = 0.000001
        self.sa_p = 0.987 # the initial probability to accept up-hill so-
        self.sa_l  = 400 
        self.sa_n = 0  
        self.sa_k = 7 # use for updating t
        self.sa_c = 100 # use for updating t
        self.sa_iter_limit = 10000;

        self.sa_init_perturb_num = 100
        self.sa_uphill_avg_cnt = 0
        self.sa_uphill_sum = 0
        self.sa_delta_cost_cnt = 0
        self.sa_delta_sum = 0
        self.sa_delta = 0
        self.sa_n_trials = 1
        self.best_mapping = None
        self.best_cost = sys.maxsize
        self.current_mapping = None
        self.current_mapping_physical_to_program = None
        self.current_cost = sys.maxsize

    def run(self, chip_dim: tuple, n_qubit: int, list_gate: list, qubit_mapping: list, target_qubits: list):
        """ 
        Run SA to find a good mapping
        """
        # print("target_qubits")
        # print(target_qubits)
        print("[INFO] Enola: Start SA-based placement")
        if self.l2:
            print("[INFO] Enola: Use l2 model for WL")
        self.initialize_param()
        self.initial_mapping = qubit_mapping
        self.target_qubits = target_qubits
        self.n_qubit = n_qubit
        length = int(math.sqrt(n_qubit)) + 4
        self.chip_dim = (min(chip_dim[0], length), min(chip_dim[1], length))
        if self.chip_dim[0] * self.chip_dim[1] < self.n_qubit:
            self.chip_dim = chip_dim
        self.list_gate = list_gate
        self.list_qubit_list_gate = [[] for i in range(self.n_qubit)]
        for i, gates in enumerate(self.list_gate):
            for j, gate in enumerate(gates):
                self.list_qubit_list_gate[gate[0]].append((i, j))
                self.list_qubit_list_gate[gate[1]].append((i, j))
        # large iteration
        for trial in range(self.sa_n_trials):
            self.init_sa_solution()
            n_reject = 0
            self.sa_n = 0;
            while self.sa_t > self.sa_t_frozen:
                self.sa_n += 1
                self.sa_delta_cost_cnt = 0
                self.sa_delta_sum = 0
                n_reject = 0
                for i in range(self.sa_l):
                    self.make_movement() # make movement and calculate cost difference
                    self.sa_delta_cost_cnt += 1
                    self.sa_delta_sum += abs(self.sa_delta)
                    # print("movement")
                    # print(self.movement)
                    # print("before: self.current_mapping")
                    # print(self.current_mapping)
                    # self.recover()
                    # print("after: self.current_mapping")
                    # print(self.current_mapping)
                    # input()
                    if self.sa_delta <= 0:
                        self.current_cost += self.sa_delta
                        if self.best_cost - self.current_cost > 1e-9:
                            self.update_optimal_sol()
                    else:   # delta > 0
                        if self.accept_worse_sol():
                            self.current_cost += self.sa_delta
                        else:
                            # undo current movement
                            self.recover()
                            n_reject += 1
                self.update_temperature()
                if self.sa_n > self.sa_iter_limit:
                    break
            print("[INFO] Enola: SA-Based Placer: Iter {}, cost: {:4f}".format(trial, self.best_cost));  
        # final refinement
        # print("qubit mapping")
        # print(self.best_mapping)
        # for i in range(self.n_qubit):
        #     for x in range(self.chip_dim[0]):
        #         for y in range(self.chip_dim[1]):
        #             self.movement = (i, x, y)
        #             self.make_movement()
        #             self.current_cost += self.sa_delta
        #             if self.best_cost - self.current_cost > 1e-9:
        #                 self.update_optimal_sol()
        #             else:
        #                 self.recover()
        # print("[INFO] SA-Based Placer: Final cost: {:4f}".format(self.best_cost)) 

    def init_sa_solution(self):
        """ 
        Random generate a placement solution 
        """
        # random generate a placement
        self.list_position = []
        for x in range(self.chip_dim[0]):
            for y in range(self.chip_dim[1]):
                if not ((x,y) in self.initial_mapping):
                    self.list_position.append((x,y))
        for i in self.target_qubits:
            self.list_position.append(self.initial_mapping[i])
        # print("self.list_position")
        # print(self.list_position)
        # bias = randrange(self.chip_dim[0] * self.chip_dim[1] - self.n_qubit)
        # list_possible_position = [i for i in range(bias, bias + self.n_qubit)]
        self.current_mapping = self.initial_mapping
        self.current_mapping_physical_to_program = [[-1 for i in range(self.chip_dim[1])] for i in range(self.chip_dim[0])]
        for i in range(self.n_qubit):
            self.current_mapping_physical_to_program[self.current_mapping[i][0]][self.current_mapping[i][1]] = i
        
        # self.current_mapping = [(randrange(self.chip_dim[0]), randrange(self.chip_dim[1])) for i in self.n_qubit]
        self.current_cost = self.get_cost()
        if self.best_cost - self.current_cost > 1e-9:
            self.update_optimal_sol()

        # Initial Perturbation 
        self.init_perturb()
        # Initialize 
        self.sa_uphill_sum = 0
        self.sa_uphill_avg_cnt = 0
        self.sa_delta_sum = 0
        self.sa_delta_cost_cnt = 0
    
    def init_perturb(self):
        uphill_sum = 0
        sum_cost = 0
        uphill_cnt = 0
        for i in range(self.sa_init_perturb_num):
            self.make_movement()
            
            self.current_cost += self.sa_delta
            if self.best_cost - self.current_cost > 1e-9:
                self.update_optimal_sol()
            if self.sa_delta > 0:
                uphill_sum += self.sa_delta;
                uphill_cnt += 1
            
            sum_cost += self.current_cost
        self.sa_t1 = (uphill_sum / uphill_cnt) / ((-1)*math.log(self.sa_p))
        self.sa_t = self.sa_t1

    def accept_worse_sol(self):
        accept = (uniform(0, 1)) <= math.exp(-(self.sa_delta)/(self.sa_t))
        return accept

    def update_temperature(self):
        if self.sa_n <= self.sa_k:
            self.sa_t = (self.sa_t1 * abs(self.sa_delta_sum) / self.sa_delta_cost_cnt) / self.sa_n / self.sa_c
        else:
            self.sa_t = (self.sa_t1 * abs(self.sa_delta_sum) / self.sa_delta_cost_cnt) / self.sa_n
        
    def make_movement(self, given_movement = False):
        # get movement
        if given_movement:
            qubit_to_move = self.movement[0]
            new_x = self.movement[1]
            new_y = self.movement[2]
        else:
            qubit_idx = randrange(len(self.target_qubits))
            qubit_to_move = self.target_qubits[qubit_idx]
            position_idx = randrange(len(self.list_position))
            new_x = self.list_position[position_idx][0]
            new_y = self.list_position[position_idx][1]
        old_x = self.current_mapping[qubit_to_move][0]
        old_y = self.current_mapping[qubit_to_move][1]
        qubit_be_affected = self.current_mapping_physical_to_program[new_x][new_y]
        self.movement = (qubit_to_move, new_x, new_y, old_x, old_y)
        # calculate original cost
        ori_cost = 0
        set_affected_gate = set()
        for gate in self.list_qubit_list_gate[qubit_to_move]:
            set_affected_gate.add(gate)
        if qubit_be_affected > -1:
            assert(qubit_be_affected in self.target_qubits)
            for gate in self.list_qubit_list_gate[qubit_be_affected]:
                set_affected_gate.add(gate)
        
        for gate in set_affected_gate:
            level = gate[0]
            weight = max((1 - 0.1 * level), 0.1)
            q0 = self.list_gate[gate[0]][gate[1]][0]
            q1 = self.list_gate[gate[0]][gate[1]][1]
            if self.l2:
                p0 = self.current_mapping[q0]
                p1 = self.current_mapping[q1]
                dis = pow((p0[0]-p1[0]), 2) + pow((p0[1]-p1[1]), 2)
            else:
                dis = math.dist(self.current_mapping[q0], self.current_mapping[q1]) # Euclidean distance    
            ori_cost += (weight * dis) # !
        
        # make movement
        self.current_mapping[qubit_to_move] = (new_x, new_y)
        self.current_mapping_physical_to_program[new_x][new_y] = qubit_to_move
        self.current_mapping_physical_to_program[old_x][old_y] = qubit_be_affected
        if qubit_be_affected > -1:
            self.current_mapping[qubit_be_affected] = (old_x, old_y)
        
        # calculate new cost
        new_cost = 0
        for gate in set_affected_gate:
            level = gate[0]
            weight = max((1 - 0.1 * level), 0.1)
            q0 = self.list_gate[gate[0]][gate[1]][0]
            q1 = self.list_gate[gate[0]][gate[1]][1]
            if self.l2:
                p0 = self.current_mapping[q0]
                p1 = self.current_mapping[q1]
                dis = pow((p0[0]-p1[0]), 2) + pow((p0[1]-p1[1]), 2)
            else:
                dis = math.dist(self.current_mapping[q0], self.current_mapping[q1]) # Euclidean distance    
            new_cost += (weight * dis) # !
            
        self.sa_delta = new_cost - ori_cost

    def recover(self):
        qubit_to_move, old_x, old_y, new_x, new_y = self.movement
        qubit_be_affected = self.current_mapping_physical_to_program[new_x][new_y]
        self.current_mapping[qubit_to_move] = (new_x, new_y)
        self.current_mapping_physical_to_program[new_x][new_y] = qubit_to_move
        self.current_mapping_physical_to_program[old_x][old_y] = qubit_be_affected
        if qubit_be_affected > -1:
            self.current_mapping[qubit_be_affected] = (old_x, old_y)
    
    def get_cost(self):
        cost = 0
        for level, gates in enumerate(self.list_gate):
            weight = max((1 - 0.1 * level), 0.1)
            tmp_cost = 0
            for gate in gates:
                if self.l2:
                    p0 = self.current_mapping[gate[0]]
                    p1 = self.current_mapping[gate[1]]
                    dis = pow((p0[0]-p1[0]), 2) + pow((p0[1]-p1[1]), 2)
                else:
                    dis = math.dist(self.current_mapping[gate[0]], self.current_mapping[gate[1]]) # Euclidean distance
                tmp_cost += dis
            cost += (tmp_cost * weight)
        # todo: calculate overlapping cost. now just keep the solution valid.
        # todo: calculate constraint violation cost
        return cost

    def update_optimal_sol(self):
        self.best_mapping = self.current_mapping.copy()
        self.best_cost = self.current_cost

def place_qubit_partial(chip_dim: tuple, n_qubit: int, list_gate: list, qubit_mapping: list, target_qubits: list, l2):
    """
    generate qubit initial layout
    """
    # print("placer: list_gate")
    # print(list_gate)
    sa_placer = SAPlacerPartial(l2)
    sa_placer.run(chip_dim, n_qubit, list_gate, qubit_mapping, target_qubits)
    return sa_placer.best_mapping