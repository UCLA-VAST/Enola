#!/usr/bin/env python
# coding: utf-8

from typing import Sequence, Tuple

from enola.scheduler.sequence_view import SequenceView
from enola.scheduler.graph import Graph
from enola.scheduler.custom_types import *
import rustworkx as rx

def find_maximal_fan(g: Graph, x: NodeType, f: NodeType) -> FanType:
    fan_x = [f]
    fan_is_maximal = False
    x_neigh = g.adj_list[x]
    fan_x_last_node = f
    while not fan_is_maximal:
        fan_is_maximal = True
        for v in x_neigh:
            if (
                v not in fan_x
                and g.edge_is_colored(x, v)
                and g.color_is_free_at_vertex(g.get_edge_color(x, v), fan_x_last_node)
            ):
                fan_x.append(v)
                fan_x_last_node = v
                fan_is_maximal = False
    return fan_x


def find_colors_cd(
    g: Graph, x: NodeType, fan_x: FanType
) -> Tuple[ColorType, ColorType]:
    l = fan_x[-1]
    c = d = 1
    while not g.color_is_free_at_vertex(c, x):
        c += 1
    while not g.color_is_free_at_vertex(d, l):
        d += 1
    return c, d


def find_and_invert_cd_path(
    g: Graph, u: NodeType, c: ColorType, d: ColorType
) -> LengthType:
    """`u` is the `X` of the maximal fan.
    returns length of cd-path"""

    path_is_maximal = False
    seen = {u}  # TODO: can be optimized, a cd-trail may also be ok.
    # so just keep track of the previous vertex.
    while not path_is_maximal:
        path_is_maximal = True
        for v in g.adj_list[u]:
            if d == g.get_edge_color(u, v) and v not in seen:
                g.set_edge_color(u, v, c)  # invert edge color
                u = v  # set vertex for next iter
                c, d = d, c  # swap colors
                path_is_maximal = False  # to loop over the new vertex
                seen.add(v)  # TODO
                # seen |= {v}  # TODO
                break  # next iter
    return len(seen) - 1


def find_w_in_fan(g: Graph, d: ColorType, fan_x: FanType) -> Tuple[IndexType, NodeType]:
    for i, u in enumerate(fan_x):
        if g.color_is_free_at_vertex(d, u):
            return i, u
    return -1, None  # this line should be unreachable


def rotate_fan(g: Graph, x: NodeType, fan_prime_x: FanViewType):
    for u, uplus in zip(fan_prime_x, fan_prime_x[1:]):
        c = g.get_edge_color(x, uplus)
        g.set_edge_color(x, u, c)
#    g.rm_edge_color(
#        x, fan_prime_x[-1]
#    )  # Redundant since we set the edge's color after this function anyway.


def graph_coloring(n_qubit: int, list_gate: list) -> int:
    g = Graph(n_qubit, len(list_gate))

    delta = -1
    for edge in list_gate:
        u = edge[0]
        v = edge[1]
        temp_degree = g.add_edge(u, v)
        if temp_degree > delta:
            delta = temp_degree

    max_color = -1
    for X, f in list_gate:
        fan_x = find_maximal_fan(g, X, f)
        c, d = find_colors_cd(g, X, fan_x)
        cd_path_len = find_and_invert_cd_path(g, X, c, d)
        w_idx, w = (
            find_w_in_fan(g, d, fan_x) if cd_path_len else (len(fan_x) - 1, fan_x[-1])
        )
        fan_x_view = SequenceView(
            fan_x
        )  # to prevent creating new List objects when slicing them.
        rotate_fan(g, X, fan_x_view[: w_idx + 1])
        g.set_edge_color(X, w, d)
        if d > max_color:
            max_color = d

    # print(delta, max_color)
    result = [[] for i in range(max_color)]
    for i, (u, v) in enumerate(list_gate):
        result[g.get_edge_color(u, v) - 1].append(i)
        # print(u, v, g.get_edge_color(u, v))

    coloring_is_optimal = max_color <= delta + 1
    return int(not coloring_is_optimal), result  # return 0 if everything went smoothly



def graph_coloring_rustworkx(n_qubit: int, list_gate: list):
    graph = rx.PyGraph()
    graph.add_nodes_from(list(range(n_qubit)))
    for edge in list_gate:
        graph.add_edge(edge[0], edge[1], edge)
    edge_colors = rx.graph_misra_gries_edge_color(graph)
    max_color = 0
    for i in edge_colors:
        max_color = max(max_color, edge_colors[i])
    result = [[] for i in range(max_color + 1)]
    for i in range(len(list_gate)):
        result[edge_colors[i]].append(i)
    return 0, result

def gate_scheduling(n_qubit: int, list_gate: list):
    """
    solve gate scheduling problem for all-commutable gate cases by graph coloring algorithm
    """
    # status, result = graph_coloring(n_qubit, list_gate)
    status, result = graph_coloring_rustworkx(n_qubit, list_gate)
    assert (
        status == 0
    ), "The proposed edge coloring is probably wrong as the output doesn't match Vizing's theorem."
    return result
