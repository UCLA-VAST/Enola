from enola.enola import Enola
import argparse
import json

if __name__ == "__main__":
    with open('./graphs.json', 'r') as f:
        graphs = json.load(f)
    parser = argparse.ArgumentParser()
    parser.add_argument('size', metavar='S', type=int, help='#qubit in graph: 30, 40, 50, 60, 70, 80, 90, 100, 200, 500, 1000, 2000, 5000, 10000.')
    parser.add_argument('id', metavar='I', type=int, help='index of the graph: 0~9.')
    parser.add_argument('--suffix', type=str,
                        help='suffix to the file name.')
    parser.add_argument('--arch', help='architecture dimension', type=int, default=16)
    parser.add_argument('--routing_strategy', help='routing strategy: mis, maximalis, maximalis_sorted', type=str, default="maximalis_sorted")
    parser.add_argument('--trivial_layout', action='store_true', default=False)
    parser.add_argument('--r2i',  help='whether reverse to initial mapping after each Rydberg stage', action='store_true', default=False)
    parser.add_argument('--window',  help='restrict vertex number in mis to 1000', action='store_true', default=False)
    parser.add_argument('--full_code',  help='generate full code for animation. The size of the code file could be large.', action='store_true', default=False)
    args = parser.parse_args()

    filename = 'rand3reg_' + str(args.size) + '_' + str(args.id)

    if args.suffix:
        filename += '_' + args.suffix
    tmp = Enola(
        filename,
        dir='./results/',
        trivial_layout = args.trivial_layout,
        routing_strategy=args.routing_strategy,
        reverse_to_initial=args.r2i,
        use_window=args.window,
        full_code=args.full_code
    )

    tmp.setArchitecture([args.arch, args.arch, args.arch, args.arch])
    if str(args.size) in graphs.keys() and args.id in range(10):
        tmp.setProgram(graphs[str(args.size)][args.id])
    else:
        import networkx as nx
        graphs[args.size] = []
        for i in range(10):
            g = list(nx.random_regular_graph(3, args.size, i).edges)
            graphs[args.size].append(g)
        tmp.setProgram(graphs[args.size][args.id])
        with open('./graphs_new.json', 'w') as f:
            json.dump(graphs, f)
        # raise ValueError(f'No such graph {args.size}_{args.id}.')
    tmp.solve(save_file=True)



