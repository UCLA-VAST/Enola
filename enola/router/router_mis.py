from networkx import maximal_independent_set, Graph
import subprocess
import time
from enola.placer.placer_parital_mapping import place_qubit_partial
import math

def compatible_2D(a, b):
    # a = (start_row, end_row, start_col, end_col)
    if a[0] == b[0] and a[1] != b[1]:
        return False
    if a[1] == b[1] and a[0] != b[0]:
        return False
    if a[0] < b[0] and a[1] >= b[1]:
        return False
    if a[0] > b[0] and a[1] <= b[1]:
        return False

    if a[2] == b[2] and a[3] != b[3]:
        return False
    if a[3] == b[3] and a[2] != b[2]:
        return False
    if a[2] < b[2] and a[3] >= b[3]:
        return False
    if a[2] > b[2] and a[3] <= b[3]:
        return False

    # if a[0] < b[0] and a[1] > b[1]:
    #     return False
    # if a[0] > b[0] and a[1] < b[1]:
    #     return False

    # if a[2] < b[2] and a[3] > b[3]:
    #     return False
    # if a[2] > b[2] and a[3] < b[3]:
    #     return False
    return True

def maximalis_solve(n, edges):
    G = Graph()
    for i in range(n):
        G.add_node(i)
    for edge in edges:
        G.add_edge(edge[0], edge[1])
    result = maximal_independent_set(G, seed=0) 
    return result

def maximalis_solve_sort(n, edges):
    # assum the vertices are sorted based on qubit distance
    is_node_conflict = [False for _ in range(n)]
    node_neighbors = {i: [] for i in range(n)}
    for edge in edges:
        node_neighbors[edge[0]].append(edge[1])
        node_neighbors[edge[1]].append(edge[0])
    result = []
    for i in range(len(is_node_conflict)):
        if is_node_conflict[i]:
            continue
        else:
            result.append(i)
            for j in node_neighbors[i]:
                is_node_conflict[j] = True
    return result

def kamis_solve(n, edges, batch):
    adj = [[] for _ in range(n)]
    for edge in edges:
        adj[edge[0]].append(edge[1])
        adj[edge[1]].append(edge[0])
    for i in range(n):
        adj[i].sort()
    with open(f"mis/tmp/mis_{batch}.in", "w") as f:
        f.write(f"{n} {len(edges)}\n")
        for i in range(n):
            tmp = ""
            for j in adj[i]:
                tmp += str(j+1)
                tmp += " "
            tmp += "\n"
            f.write(tmp)
    with open(f"mis/tmp/mis_{batch}.log", "w") as f:
        subprocess.run(
            ["mis/redumis", f"mis/tmp/mis_{batch}.in", "--output", f"mis/tmp/mis_{batch}.out", "--time_limit", "3600"], stdout=f)
    with open(f"mis/tmp/mis_{batch}.out", "r") as f:
        lines = f.readlines()
    return [i for i, line in enumerate(lines) if line.startswith("1")]


def route_qubit_mis(chip_dim: tuple, n_q: int, index_list_gate:int, list_gates: list, qubit_mapping: list, \
                    routing_strategy: str, reverse_to_initial: bool, \
                    l2: bool, use_window: bool):
    # print("list_gates")
    # print(list_gates[index_list_gate])
    start_time = time.time()
    layer_gate_number = []
    layer_aod = [] # use to record aods per layer
    set_aod_qubit = set()
    layer_time = [] # use to record the solving time per layer
    initial_mapping = qubit_mapping.copy()
    layers = []
    time_placement = 0
    vector_threshold = 1000
    layers.append(
    {
        "qubits": [{
            "id": i,
            "a": 0,
            "x": qubit_mapping[i][0],
            "y": qubit_mapping[i][1],
            "c": qubit_mapping[i][0],
            "r": qubit_mapping[i][1],
        } for i in range(n_q)],
        "gates": []
    }
    )
    # ! sort remain_graph based on qubit distance if using maximal is
    remain_graph = list_gates[index_list_gate].copy()
    if not(routing_strategy == "mis" or routing_strategy == "maximalis"):
        # print("remain_graph")
        # print(remain_graph)
        # print("qubit_mapping")
        # print(qubit_mapping)
        remain_graph = sorted(remain_graph, key=lambda x: (math.dist(qubit_mapping[x[0]], qubit_mapping[x[1]])), reverse=True)
        # print("remain_graph")
        # print(remain_graph)
    batch = 0
    
    last_layer = {
                    "qubits": [{
                        "id": i,
                        "a": 0,
                        "x": qubit_mapping[i][0],
                        "y": qubit_mapping[i][1],
                        "c": qubit_mapping[i][0],
                        "r": qubit_mapping[i][1],
                    } for i in range(n_q)],
                    "gates": []
                }
    gates_dict = []
    while remain_graph:
        vectors = []
        if use_window:
            vector_length = min(vector_threshold, 2*len(remain_graph))
            vectors = [(0,0,0,0, ) for _ in range(vector_length)]
            # t_tmp = time.time()
            i = 0
            for edge in remain_graph:
                if i + 1 < vector_length:
                    (q0_row, q0_col) = qubit_mapping[edge[0]]
                    (q1_row, q1_col) = qubit_mapping[edge[1]]
                    vectors[i] = (q0_col, q1_col, q0_row, q1_row, )
                    vectors[i+1] = (q1_col, q0_col, q1_row, q0_row, )
                    i += 2
        else:
            vectors = [(0,0,0,0, ) for _ in range(2*len(remain_graph))]
            i = 0
            for edge in remain_graph:
                (q0_row, q0_col) = qubit_mapping[edge[0]]
                (q1_row, q1_col) = qubit_mapping[edge[1]]
                vectors[i] = (q0_col, q1_col, q0_row, q1_row, )
                vectors[i+1] = (q1_col, q0_col, q1_row, q0_row, )
                i += 2
            
        violations = []
        for i in range(len(vectors)):
            for j in range(i+1, len(vectors)):
                if not compatible_2D(vectors[i], vectors[j]):
                    violations.append((i, j))
        # print("time for graph construction: {}".format(time.time() - t_tmp))
        # print("number of violation: {}".format(len(violations)))
        # t_tmp = time.time()
        if routing_strategy == "mis":
            execute_gates = kamis_solve(len(vectors), violations, batch)
        elif routing_strategy == "maximalis":
            execute_gates = maximalis_solve(len(vectors), violations)
        else:
            execute_gates = maximalis_solve_sort(len(vectors), violations)
        # print("time for mis solving: {}".format(time.time() - t_tmp))
        # t_tmp = time.time()
        layer_gate_number.append(len(execute_gates))
        aod_qubits = []
        target_qubits = {}
        # tmp_moving_distance = []
        for i in execute_gates:
            if i % 2 == 0:
                aod_qubits.append(remain_graph[i//2][0])
                set_aod_qubit.add(remain_graph[i//2][0])
                target_qubits[remain_graph[i//2][0]] = remain_graph[i//2][1]
            else:
                aod_qubits.append(remain_graph[i//2][1])
                set_aod_qubit.add(remain_graph[i//2][1])
                target_qubits[remain_graph[i//2][1]] = remain_graph[i//2][0]
        layer_aod.append(aod_qubits)

        tmp = [remain_graph[g // 2] for g in execute_gates]
        # gates_dict = []
        for gate in tmp:
            if gate in list_gates[index_list_gate]:
                gates_dict.append(
                    {
                        "id": list_gates[index_list_gate].index(gate),
                        "q0": gate[0],
                        "q1": gate[1]
                    }
                )
            else:
                gates_dict.append(
                    {
                        "id": list_gates[index_list_gate].index((gate[1], gate[0])),
                        "q0": gate[1],
                        "q1": gate[0]
                    }
                )
        
        for i in range(n_q):
            if i in aod_qubits:
                qubit_mapping[i] = qubit_mapping[target_qubits[i]]
        # for i in aod_qubits:
        #     qubit_mapping[i] = qubit_mapping[target_qubits[i]]
        layers.append(
            {
                "qubits": [{
                    "id": i,
                    "a": 0,
                    "x": qubit_mapping[i][0],
                    "y": qubit_mapping[i][1],
                    "c": qubit_mapping[i][0],
                    "r": qubit_mapping[i][1],
                } for i in range(n_q)],
                # "gates": gates_dict,
                "gates": [],
            }
        )
        for qubit in layers[-2]["qubits"]:
            if qubit["id"] in layer_aod[-1]:
                qubit["a"] = 1
        # print("layers")
        # print(layers)

        tmp = [edge for g, edge in enumerate(
            remain_graph) if 2*g not in execute_gates and 2*g+1 not in execute_gates]
        remain_graph = tmp
        batch += 1
        # print("time for post processsing: {}".format(time.time() - t_tmp))
        layer_time.append(float(time.time() - start_time))

    layers[-1]["gates"] = gates_dict
    # print("data")
    # print(data)
    if index_list_gate + 1 < len(list_gates) or reverse_to_initial:
        if reverse_to_initial:
            for q in range(len(layers[-1]["qubits"])):
                layers[-1]["qubits"][q]["a"] = layers[-2]["qubits"][q]["a"]
            
            reverse_layers = []
            # print("update layer {} by layer {}".format(len(layers)-1, len(layers)-2))
            for q in range(len(layers[-1]["qubits"])):
                layers[len(layers)-1]["qubits"][q]["a"] = layers[len(layers)-2]["qubits"][q]["a"]
                layers[len(layers)-1]["qubits"][q]["c"] = layers[len(layers)-2]["qubits"][q]["c"]
                layers[len(layers)-1]["qubits"][q]["r"] = layers[len(layers)-2]["qubits"][q]["r"]
            for i in range(len(layers)-2, 0, -1):
                # reverse_layers.append(dict(layers[i]))
                reverse_layers.append({
                                    "qubits": [{
                                        "id": j,
                                        "a": 0,
                                        "x": layers[i]["qubits"][j]["x"],
                                        "y": layers[i]["qubits"][j]["y"],
                                        "c": layers[i]["qubits"][j]["c"],
                                        "r": layers[i]["qubits"][j]["r"],
                                    } for j in range(n_q)],
                                    "gates": [],
                                })
                reverse_layers[len(reverse_layers)-1]["gates"] = []
                # print("update reverse layer {} by layer {}".format(len(reverse_layers)-1, i-1))
                for q in range(len(reverse_layers[-1]["qubits"])):
                    reverse_layers[len(reverse_layers)-1]["qubits"][q]["a"] = layers[i-1]["qubits"][q]["a"]
                    reverse_layers[len(reverse_layers)-1]["qubits"][q]["c"] = layers[i-1]["qubits"][q]["c"]
                    reverse_layers[len(reverse_layers)-1]["qubits"][q]["r"] = layers[i-1]["qubits"][q]["r"]

                        
            layers = layers + reverse_layers
            layers.append(last_layer)
        else:
            # find final mapping here
            print("[INFO] Enola: Finding a mapping for Rydberg stage {}/{}".format(index_list_gate+2,len(list_gates)))
            list_aod_qubit = list(set_aod_qubit)
            list_gate_next_layers = list_gates[index_list_gate:] # consider current gate
            t_tmp = time.time()
            initial_mapping = place_qubit_partial(chip_dim, n_q, list_gate_next_layers, initial_mapping, list_aod_qubit, l2)
            time_placement += (time.time() - t_tmp) 
            if not(routing_strategy == "mis" or routing_strategy == "maximalis"):
                list_aod_qubit = sorted(list_aod_qubit, key=lambda x: (math.dist(qubit_mapping[x], initial_mapping[x])), reverse=True)
            while len(list_aod_qubit) > 0:
                vectors = []
            
                if use_window:
                    vector_length = min(vector_threshold, len(list_aod_qubit))
                    vectors = [(0,0,0,0, ) for _ in range(vector_length)]
                    # t_tmp = time.time()
                    for i, q in enumerate(list_aod_qubit):
                        if i < vector_length:
                            vectors[i] = (qubit_mapping[q][0], initial_mapping[q][0], qubit_mapping[q][1], initial_mapping[q][1])
                else:
                    vectors = [(0,0,0,0, ) for _ in range(len(list_aod_qubit))]
                    for i, q in enumerate(list_aod_qubit):
                        vectors[i] = (qubit_mapping[q][0], initial_mapping[q][0], qubit_mapping[q][1], initial_mapping[q][1])

                violations = []
                for i in range(len(vectors)):
                    for j in range(i+1, len(vectors)):
                        if not compatible_2D(vectors[i], vectors[j]):
                            violations.append((i, j))
                if routing_strategy == "mis":
                    moved_qubits = kamis_solve(len(vectors), violations, batch)
                else:
                    moved_qubits = maximalis_solve(len(vectors), violations)
                aod_qubits = []
                for i in moved_qubits:
                    aod_qubits.append(list_aod_qubit[i])
                    qubit_mapping[list_aod_qubit[i]] = initial_mapping[list_aod_qubit[i]]
                layer_aod.append(aod_qubits)

                layers.append(
                    {
                        "qubits": [{
                            "id": i,
                            "a": 0,
                            "x": qubit_mapping[i][0],
                            "y": qubit_mapping[i][1],
                            "c": qubit_mapping[i][0],
                            "r": qubit_mapping[i][1],
                        } for i in range(n_q)],
                        "gates": [],
                    }
                )
                for qubit in layers[-2]["qubits"]:
                    if qubit["id"] in layer_aod[-1]:
                        qubit["a"] = 1

                tmp = [q for q in list_aod_qubit if  list_aod_qubit.index(q) not in moved_qubits]
                list_aod_qubit = tmp
                batch += 1


    for q in range(len(layers[-1]["qubits"])):
        layers[-1]["qubits"][q]["a"] = 0

    data = {
        "runtime": float(time.time() - start_time),
        "no_transfer": False,
        "layers": layers,
        "n_q": n_q,
        "g_q": list_gates,
    }
    return data, initial_mapping, time_placement