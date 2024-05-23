# coding: utf-8

from enola.scheduler.custom_types import *


class Graph:
    def __init__(self, n: int, m: int):
        self.n = n
        self.m = m
        self.adj_mat = [ColorType()] * ((n * (n + 1)) >> 1)
        self.adj_list = [[] for _ in range(n)]

    def _idx(self, i: NodeType, j: NodeType) -> IndexType:
        if i > j:
            i, j = j, i
        #        n__i = n - i
        #        row_offset = ((self.n * (self.n + 1)) - (n__i * (n__i + 1))) >> 1
        #        row_offset = (i * (2 * self.n - i + 1)) >> 1
        #        col_offset = j - i
        #        return col_offset + row_offset
        return ((i * (2 * self.n - i + 1)) >> 1) + j - i

    def _setColor(self, u: NodeType, v: NodeType, c: ColorType):
        self.adj_mat[self._idx(u, v)] = c

    def _getColor(self, u: NodeType, v: NodeType) -> ColorType:
        return self.adj_mat[self._idx(u, v)]

    def add_edge(self, u: NodeType, v: NodeType) -> LengthType:
        v_neigh = self.adj_list[v]
        u_neigh = self.adj_list[u]
        u_neigh.append(v)
        v_neigh.append(u)
        
        len_u = len(u_neigh)
        len_v = len(v_neigh)
        if len_u > len_v:
            return len_u
        return len_v

    def get_edge_color(self, u: NodeType, v: NodeType) -> ColorType:
        return self._getColor(u, v)

    def edge_is_colored(self, u: NodeType, v: NodeType) -> bool:
        return self.get_edge_color(u, v) > 0

    def rm_edge_color(self, u: NodeType, v: NodeType):
        self._setColor(u, v, 0)

    def set_edge_color(self, u: NodeType, v: NodeType, c: ColorType):
        self._setColor(u, v, c)

    def color_is_free_at_vertex(self, c: ColorType, u: NodeType) -> bool:
        for v in self.adj_list[u]:  # Linear search :(
            if c == self.get_edge_color(u, v):
                return False
        return True
