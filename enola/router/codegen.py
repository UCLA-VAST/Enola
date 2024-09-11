from typing import Sequence, Mapping, Any
import networkx as nx
import argparse
from abc import ABC, abstractmethod
import time 


# physics constants
R_B = 6  # rydberg range
AOD_SEP = 2  # min AOD separation
RYD_SEP = 15  # sufficient distance to avoid Rydberg
SITE_SLMS = 2  # number of SLMs in a site
SLM_SEP = AOD_SEP  # separation of SLMs inside a site
SITE_WIDTH = 4  # total width of SLMs in a site
X_SITE_SEP = RYD_SEP + SITE_WIDTH  # separation of sites in X direction
Y_SITE_SEP = RYD_SEP  # separation of sites in Y direction

# padding of the figure
X_LOW_PAD = 2 * AOD_SEP
Y_LOW_PAD = 4 * AOD_SEP
X_HIGH_PAD = 2 * AOD_SEP
Y_HIGH_PAD = 4 * AOD_SEP

# constants for animation
FPS = 24  # frames per second
INIT_FRM = 24  # initial empty frames
PT_MICRON = 8  # scaling factor: points per micron
MUS_PER_FRM = 8  # microseconds per frame
T_RYDBERG = 0.36  # microseconds for Rydberg
T_ACTIVATE = 15  # microseconds for (de)activating AOD

global_dict = {
    "full_code": False
}

# class for physical entities: qubits and AOD rows/cols

class Qubit():
    def __init__(self, id: int):
        self.array = 'SLM'
        self.id = id
        self.c = -1  # AOD coloumn index
        self.r = -1  # AOD row index
        self.x = -X_LOW_PAD-1  # real X coordinates in um
        self.y = -Y_LOW_PAD-1  # real Y coordinates in um


class Row():
    def __init__(self, id: int):
        self.id = id
        self.active = False
        self.y = -Y_LOW_PAD-1  # real Y coordinates in um


class Col():
    def __init__(self, id: int):
        self.id = id
        self.active = False
        self.x = -X_LOW_PAD-1  # real X coordinates in um


class Inst(ABC):
    """abstract class of DPQA instructions.

    In general, the __init__ of specific instruction classes looks like
        def __init__(self, *):
            super().__init__(*)
            self.verify(*)
            self.operate(*)
            super().write_code(*)
    """

    def __init__(
            self,
            type: str,
            prefix = None,
            stage: int = -1,
            reduced_keys: Sequence[str] = [],
    ):
        """init method for instructions.

        Args:
            type (str): 
            prefix (str | None, optional): provide the big operation.
                this Inst belongs to. Defaults to None.
            stage (int, optional): stage the Inst belongs to. Defaults to -1.
            reduced_keys (Sequence[str], optional): data to keep in emit()
                from emit_full(). Defaults to [].
        """
        self.type = type
        self.name = prefix + ':' + type if prefix else type
        self.stage = stage
        self.reduced_keys = reduced_keys + ['type', ]
        self.duration = -1
        self.code = {'type': self.type, 'name': self.name, }

    def write_code(
            self,
            col_objs: Sequence[Col],
            row_objs: Sequence[Row],
            qubit_objs: Sequence[Qubit],
            data: Mapping[str, Any],
    ):
        """write self.code with provided info.

        Args:
            col_objs (Sequence[Col]): Col objects used.
            row_objs (Sequence[Row]): Row objects used.
            qubit_objs (Sequence[Qubit]): Qubit objects used.
            data (Mapping[str, Any]): other info for the Inst.
        """
        for k, v in data.items():
            self.code[k] = v

        # ! modified to reduce file size
        self.code['state'] = {}
        if global_dict['full_code']:
            # get the current state of DPQA
            curr = {}
            curr['qubits'] = [
                {
                    'id': q.id,
                    'x': q.x,
                    'y': q.y,
                    'array': q.array,
                    'c': q.c,
                    'r': q.r
                } for q in qubit_objs
            ]
            curr['cols'] = [
                {
                    'id': c.id,
                    'active': c.active,
                    'x': c.x
                } for c in col_objs
            ]
            curr['rows'] = [
                {
                    'id': r.id,
                    'active': r.active,
                    'y': r.y
                } for r in row_objs
            ]
            self.code['state'] = curr

    @abstractmethod
    def verify(self):
        """verification of instructions. This is abstract because we require
        each child class to provide its own verification method.
        """
        pass

    def operate(self):
        """perform operation of instructions on Col, Row, and Qubit objects."""
        pass

    def emit(self) -> Sequence[Mapping[str, Any]]:
        """emit the minimum code for executing instructions.

        Returns:
            Sequence[Mapping[str, Any]]: code in a dict.
        """
        return ({k: self.code[k] for k in self.reduced_keys}, )

    def emit_full(self) -> Sequence[Mapping[str, Any]]:
        """emit the code with full info for any purpose.

        Returns:
            Sequence[Mapping[str, Any]]: code with full info in a dict.
        """
        return (self.code, )

    def is_trivial(self) -> bool:
        return True if self.duration == 0 else False

    def remove_trivial_insts(self):
        # this is used in ComboInst. Added here for convience.
        pass


# classes for basic instructions: Init, Move, Activate, Deactivate, Rydberg
# todo: add class Raman (single-qubit gates)


class Init(Inst):
    def __init__(self,
                 col_objs: Sequence[Col],
                 row_objs: Sequence[Row],
                 qubit_objs: Sequence[Qubit],
                 slm_qubit_idx: Sequence[int] = [],
                 slm_qubit_xys: Sequence[Sequence[int]] = [],
                 aod_qubit_idx: Sequence[int] = [],
                 aod_qubit_crs: Sequence[Sequence[int]] = [],
                 aod_col_act_idx: Sequence[int] = [],
                 aod_col_xs: Sequence[int] = [],
                 aod_row_act_idx: Sequence[int] = [],
                 aod_row_ys: Sequence[int] = [],
                 data: Mapping[str, Any] = {},):
        super().__init__(
            'Init',
            reduced_keys=[
                'slm_qubit_idx', 'slm_qubit_xys', 'aod_qubit_idx',
                'aod_qubit_crs', 'aod_col_act_idx', 'aod_col_xs',
                'aod_row_act_idx', 'aod_row_ys', 'n_q',
                'x_high', 'y_high', 'c_high', 'r_high'
            ]
        )
        for k, v in data.items():
            self.code[k] = v
        self.all_slms = []
        self.verify(slm_qubit_idx,
                    slm_qubit_xys,
                    aod_qubit_idx,
                    aod_qubit_crs,
                    aod_col_act_idx,
                    aod_col_xs,
                    aod_row_act_idx,
                    aod_row_ys,)
        self.operate(col_objs,
                     row_objs,
                     qubit_objs,
                     slm_qubit_idx,
                     slm_qubit_xys,
                     aod_qubit_idx,
                     aod_qubit_crs,
                     aod_col_act_idx,
                     aod_col_xs,
                     aod_row_act_idx,
                     aod_row_ys,)
        super().write_code(
            col_objs,
            row_objs,
            qubit_objs,
            {
                'duration': INIT_FRM,
                'slm_qubit_idx': slm_qubit_idx,
                'slm_qubit_xys': slm_qubit_xys,
                'aod_qubit_idx': aod_qubit_idx,
                'aod_qubit_crs': aod_qubit_crs,
                'aod_col_act_idx': aod_col_act_idx,
                'aod_col_xs': aod_col_xs,
                'aod_row_act_idx': aod_row_act_idx,
                'aod_row_ys': aod_row_ys,
            })

    def add_slms(self, slms):
        for slm in slms:
            if slm not in self.all_slms:
                self.all_slms.append(slm)

    def verify(
            self,
            slm_qubit_idx: Sequence[int],
            slm_qubit_xys: Sequence[Sequence[int]],
            aod_qubit_idx: Sequence[int],
            aod_qubit_crs: Sequence[Sequence[int]],
            aod_col_act_idx: Sequence[int],
            aod_col_xs: Sequence[int],
            aod_row_act_idx: Sequence[int],
            aod_row_ys: Sequence[int],
    ):
        a = len(slm_qubit_idx)
        b = len(slm_qubit_xys)
        if a != b:
            raise ValueError(
                f'{self.name}: SLM qubit arguments invalid {a} idx, {b} xys.')
        for i in slm_qubit_xys:
            if len(i) != 2:
                raise ValueError(f'{self.name}: SLM qubit xys {i} invalid.')
        for i in range(len(slm_qubit_xys)):
            for j in range(i+1, len(slm_qubit_xys)):
                if slm_qubit_xys[i] == slm_qubit_xys[j]:
                    raise ValueError(
                        f'{self.name}: SLM qubits {slm_qubit_idx[i]} '
                        f'and {slm_qubit_idx[j]} xys are the same.'
                    )
        # todo: the case when not all atoms are in SLM

    def operate(
            self,
            col_objs: Sequence[Col],
            row_objs: Sequence[Row],
            qubit_objs: Sequence[Qubit],
            slm_qubit_idx: Sequence[int],
            slm_qubit_xys: Sequence[Sequence[int]],
            aod_qubit_idx: Sequence[int],
            aod_qubit_crs: Sequence[Sequence[int]],
            aod_col_act_idx: Sequence[int],
            aod_col_xs: Sequence[int],
            aod_row_act_idx: Sequence[int],
            aod_row_ys: Sequence[int],
    ):
        for i, q_id in enumerate(slm_qubit_idx):
            qubit_objs[q_id].array = 'SLM'
            (qubit_objs[q_id].x, qubit_objs[q_id].y) = slm_qubit_xys[i]
            self.all_slms.append(slm_qubit_xys[i])

        for i, col_id in enumerate(aod_col_act_idx):
            col_objs[col_id].activate = True
            col_objs[col_id].x = aod_col_xs[i]

        for i, row_id in enumerate(aod_col_act_idx):
            row_objs[row_id].activate = True
            row_objs[row_id].y = aod_row_ys[i]

        for i, q_id in enumerate(aod_qubit_idx):
            qubit_objs[q_id].array = 'AOD'
            (qubit_objs[q_id].c, qubit_objs[q_id].r) = aod_qubit_crs[i]
            qubit_objs[q_id].x = col_objs[qubit_objs[q_id].c].x
            qubit_objs[q_id].y = row_objs[qubit_objs[q_id].r].y

    def emit_full(self):
        # all the used SLMs are counted during the whole codegen process,
        # so the emit_full of Init needs to add this info
        self.code['all_slms'] = self.all_slms
        return super().emit_full()


class Move(Inst):
    def __init__(self,
                 s: int,
                 col_objs: Sequence[Col],
                 row_objs: Sequence[Row],
                 qubit_objs: Sequence[Qubit],
                 col_idx: Sequence[int] = [],
                 col_begin: Sequence[int] = [],
                 col_end: Sequence[int] = [],
                 row_idx: Sequence[int] = [],
                 row_begin: Sequence[int] = [],
                 row_end: Sequence[int] = [],
                 prefix: str = ''):
        super().__init__('Move', prefix=prefix, stage=s)
        self.verify(
            col_objs,
            row_objs,
            col_idx,
            col_begin,
            col_end,
            row_idx,
            row_begin,
            row_end,
        )
        data = self.operate(
            col_objs,
            row_objs,
            qubit_objs,
            col_idx,
            col_begin,
            col_end,
            row_idx,
            row_begin,
            row_end,)
        super().write_code(col_objs, row_objs, qubit_objs, data)

    def operate(
            self,
            col_objs: Sequence[Col],
            row_objs: Sequence[Row],
            qubit_objs: Sequence[Qubit],
            col_idx: Sequence[int],
            col_begin: Sequence[int],
            col_end: Sequence[int],
            row_idx: Sequence[int],
            row_begin: Sequence[int],
            row_end: Sequence[int],
    ) -> Mapping[str, Any]:
        data = {}
        # calculate the max  move distance of columns
        data['cols'] = []
        max_distance = 0
        for i in range(len(col_idx)):
            distance = abs(col_end[i]-col_begin[i])
            if distance > 0:
                data['cols'].append(
                    {'id': col_idx[i],
                     'shift': col_end[i]-col_begin[i],
                     'begin': col_begin[i],
                     'end': col_end[i]})
                col_objs[col_idx[i]].x = col_end[i]
                max_distance = max(max_distance, distance)

        # calculate the max  move distance of rows
        data['rows'] = []
        for i in range(len(row_idx)):
            distance = abs(row_end[i]-row_begin[i])
            if distance > 0:
                data['rows'].append(
                    {'id': row_idx[i],
                     'shift': row_end[i]-row_begin[i],
                     'begin': row_begin[i],
                     'end': row_end[i]})
                row_objs[row_idx[i]].y = row_end[i]
                max_distance = max(max_distance, distance)

        # movement time per Bluvstein et al. units are us and um.
        self.duration = 200*((max_distance/110)**(1/2))
        data['duration'] = self.duration

        for qubit_obj in qubit_objs:
            if qubit_obj.array == 'AOD':
                qubit_obj.x = col_objs[qubit_obj.c].x
                qubit_obj.y = row_objs[qubit_obj.r].y

        return data

    def verify(
            self,
            col_objs: Sequence[Col],
            row_objs: Sequence[Row],
            col_idx: Sequence[int],
            col_begin: Sequence[int],
            col_end: Sequence[int],
            row_idx: Sequence[int],
            row_begin: Sequence[int],
            row_end: Sequence[int],
    ):
        a = len(col_idx)
        b = len(col_begin)
        c = len(col_end)
        if not (a == b and a == c):
            raise ValueError(
                f'{self.name}: col arguments invalid'
                f' {a} idx, {b} begin, {c} end.'
            )
        a = len(row_idx)
        b = len(row_begin)
        c = len(row_end)
        if not (a == b and a == c):
            raise ValueError(
                f'{self.name}: row arguments invalid'
                f' {a} idx, {b} begin, {c} end.'
            )

        activated_col_idx = []
        activated_col_xs = []
        for col_obj in col_objs:
            if col_obj.active:
                if (activated_col_idx
                        and col_obj.x < activated_col_xs[-1] + AOD_SEP):
                    raise ValueError(
                        f'{self.name}: col beginning position invalid'
                        f' col {col_obj.id} at x={col_obj.x} while '
                        f'col {activated_col_idx[-1]} at'
                        f' x={activated_col_xs[-1]}.'
                    )
                activated_col_idx.append(col_obj.id)
                activated_col_xs.append(col_obj.x)
        for i, moving_col_id in enumerate(col_idx):
            if moving_col_id not in activated_col_idx:
                raise ValueError(
                    f'{self.name}: col {moving_col_id} to move'
                    f' is not activated.'
                )
            j = activated_col_idx.index(moving_col_id)
            if col_begin[i] != activated_col_xs[j]:
                raise ValueError(
                    f'{self.name}: col {moving_col_id} beginning x not agree.')
            activated_col_xs[j] = col_end[i]
        for i in range(1, len(activated_col_xs)):
            if activated_col_xs[i - 1] + AOD_SEP > activated_col_xs[i]:
                raise ValueError(
                    f'{self.name}: col ending position invalid'
                    f' col {activated_col_idx[i-1]} at '
                    f'x={activated_col_xs[i-1]} while '
                    f'col {activated_col_idx[i]} at x={activated_col_xs[i]}.')

        activated_row_idx = []
        activated_row_ys = []
        for row_obj in row_objs:
            if row_obj.active:
                if (activated_row_idx
                        and row_obj.y < activated_row_ys[-1] + AOD_SEP):
                    raise ValueError(
                        f'{self.name}: row beginning position invalid '
                        f'row {row_obj.id} at y={row_obj.y} while '
                        f'row {activated_row_idx[-1]} at '
                        f'y={activated_row_ys[-1]}.'
                    )
                activated_row_idx.append(row_obj.id)
                activated_row_ys.append(row_obj.y)
        for i, moving_row_id in enumerate(row_idx):
            if moving_row_id not in activated_row_idx:
                raise ValueError(
                    f'{self.name}: row {moving_row_id} to move '
                    f'is not activated.'
                )
            j = activated_row_idx.index(moving_row_id)
            if row_begin[i] != activated_row_ys[j]:
                raise ValueError(
                    f'{self.name}: row {moving_row_id} beginning y not agree.')
            activated_row_ys[j] = row_end[i]
        for i in range(1, len(activated_row_ys)):
            if activated_row_ys[i - 1] + AOD_SEP > activated_row_ys[i]:
                raise ValueError(
                    f'{self.name}: row ending position invalid '
                    f'row {activated_row_idx[i-1]} at '
                    f'y={activated_row_ys[i-1]} while '
                    f'row {activated_row_idx[i]} at y={activated_row_ys[i]}.')

    def emit(self):
        code = {'type': self.type}
        code['cols'] = [{k: v for k, v in tmp.items() if k in [
            'id', 'shift']} for tmp in self.code['cols']]
        code['rows'] = [{k: v for k, v in tmp.items() if k in [
            'id', 'shift']} for tmp in self.code['rows']]
        return (code,)


class Activate(Inst):
    def __init__(self,
                 s: int,
                 col_objs: Sequence[Col],
                 row_objs: Sequence[Row],
                 qubit_objs: Sequence[Qubit],
                 col_idx: Sequence[int] = [],
                 col_xs: Sequence[int] = [],
                 row_idx: Sequence[int] = [],
                 row_ys: Sequence[int] = [],
                 pickup_qs: Sequence[int] = [],
                 prefix: str = '',):
        super().__init__(
            'Activate',
            prefix=prefix,
            stage=s,
            reduced_keys=['col_idx', 'col_xs', 'row_idx', 'row_ys']
        )
        self.verify(
            col_objs,
            row_objs,
            qubit_objs,
            col_idx,
            col_xs,
            row_idx,
            row_ys,
            pickup_qs,
        )
        self.operate(
            col_objs,
            row_objs,
            qubit_objs,
            col_idx,
            col_xs,
            row_idx,
            row_ys,
            pickup_qs,
        )
        super().write_code(
            col_objs,
            row_objs,
            qubit_objs,
            {
                'col_idx': col_idx,
                'col_xs': col_xs,
                'row_idx': row_idx,
                'row_ys': row_ys,
                'pickup_qs': pickup_qs,
                'duration': T_ACTIVATE
            }
        )

    def operate(
            self,
            col_objs: Sequence[Col],
            row_objs: Sequence[Row],
            qubit_objs: Sequence[Qubit],
            col_idx: Sequence[int],
            col_xs: Sequence[int],
            row_idx: Sequence[int],
            row_ys: Sequence[int],
            pickup_qs: Sequence[int],
    ):
        for i in range(len(col_idx)):
            col_objs[col_idx[i]].active = True
            col_objs[col_idx[i]].x = col_xs[i]
        for i in range(len(row_idx)):
            row_objs[row_idx[i]].active = True
            row_objs[row_idx[i]].y = row_ys[i]
        for q_id in pickup_qs:
            qubit_objs[q_id].array = 'AOD'

    def verify(
            self,
            col_objs: Sequence[Col],
            row_objs: Sequence[Row],
            qubit_objs: Sequence[Qubit],
            col_idx: Sequence[int],
            col_xs: Sequence[int],
            row_idx: Sequence[int],
            row_ys: Sequence[int],
            pickup_qs: Sequence[int],
    ):
        a = len(col_idx)
        b = len(col_xs)
        if a != b:
            raise ValueError(
                f'{self.name}: col arguments invalid {a} idx, {b} xs')
        a = len(row_idx)
        b = len(row_ys)
        if a != b:
            raise ValueError(
                f'f{self.name}: row arguments invalid {a} idx, {b} ys.')

        for i in range(len(col_idx)):
            if col_objs[col_idx[i]].active:
                raise ValueError(
                    f'{self.name}: col {col_idx[i]} already activated.')
            for j in range(col_idx[i]):
                if (col_objs[j].active
                        and col_objs[j].x > col_xs[i] - AOD_SEP):
                    raise ValueError(
                        f'{self.name}: col {j} at x={col_objs[j].x} is '
                        f'too left for col {col_idx[i]}'
                        f' to activate at x={col_xs[i]}.'
                    )
            for j in range(col_idx[i] + 1, len(col_objs)):
                if (col_objs[j].active
                        and col_objs[j].x - AOD_SEP < col_xs[i]):
                    raise ValueError(
                        f'{self.name}: col {j} at x={col_objs[j].x} is '
                        f'too right for col {col_idx[i]} '
                        f'to activate at x={col_xs[i]}.'
                    )
        for i in range(len(row_idx)):
            if row_objs[row_idx[i]].active:
                raise ValueError(
                    f'{self.name}: row {row_idx[i]} already activated.')
            for j in range(row_idx[i]):
                if (row_objs[j].active
                        and row_objs[j].y > row_ys[i] - AOD_SEP):
                    raise ValueError(
                        f'{self.name}: row {j} at y={row_objs[j].y} is '
                        f'too high for row {row_idx[i]} '
                        f'to activate at y={row_ys[i]}.'
                    )
            for j in range(row_idx[i] + 1, len(row_objs)):
                if (row_objs[j].active
                        and row_objs[j].y-AOD_SEP < row_ys[i]):
                    raise ValueError(
                        f'{self.name}: row {j} at y={col_objs[j].y} is '
                        f'too low for row {row_idx[i]} '
                        f'to activate at y={row_ys[i]}.'
                    )

        active_xys = []  # the traps that are newly activted by this Activate
        active_xs = [col.x for col in col_objs if col.active]
        active_ys = [row.y for row in row_objs if row.active]
        for x in active_xs:
            for y in row_ys:
                active_xys.append((x, y))
        for y in active_ys:
            for x in col_xs:
                active_xys.append((x, y))
        for x in col_xs:
            for y in row_ys:
                active_xys.append((x, y))
        
        # print("active_xys")
        # print(active_xys)
        # print("active_xs")
        # print(active_xs)
        # print("active_ys")
        # print(active_ys)
        # print("pickup_qs")
        # print(pickup_qs)

        for q_id in range(len(qubit_objs)):
            if q_id in pickup_qs:
                if (qubit_objs[q_id].x, qubit_objs[q_id].y) not in active_xys:
                    raise ValueError(
                        f'{self.name}: q {q_id} not picked up '
                        f'by col {qubit_objs[q_id].c} '
                        f'row {qubit_objs[q_id].r} at '
                        f'x={qubit_objs[q_id].x} y={qubit_objs[q_id].y}.'
                    )
            else:
                if (qubit_objs[q_id].x, qubit_objs[q_id].y) in active_xys:
                    raise ValueError(
                        f'{self.name}: q {q_id} wrongfully picked up by '
                        f'col {qubit_objs[q_id].c} row {qubit_objs[q_id].r}'
                        f' at x={qubit_objs[q_id].x} y={qubit_objs[q_id].y}.')


class Deactivate(Inst):
    def __init__(self,
                 s: int,
                 col_objs: Sequence[Col],
                 row_objs: Sequence[Row],
                 qubit_objs: Sequence[Qubit],
                 col_idx: Sequence[int] = [],
                 col_xs: Sequence[int] = [],
                 row_idx: Sequence[int] = [],
                 row_ys: Sequence[int] = [],
                 dropoff_qs: Sequence[int] = [],
                 prefix: str = '',):
        super().__init__(
            'Deactivate',
            prefix=prefix,
            stage=s,
            reduced_keys=['col_idx', 'row_idx']
        )
        self.verify(
            col_objs,
            row_objs,
            qubit_objs,
            col_idx,
            col_xs,
            row_idx,
            row_ys,
            dropoff_qs
        )
        self.operate(
            col_objs,
            row_objs,
            qubit_objs,
            col_idx,
            row_idx,
            dropoff_qs
        )
        super().write_code(
            col_objs,
            row_objs,
            qubit_objs,
            {
                'col_idx': col_idx,
                'col_xs': col_xs,
                'row_idx': row_idx,
                'row_ys': row_ys,
                'dropoff_qs': dropoff_qs,
                'duration': T_ACTIVATE
            }
        )

    def operate(
        self,
        col_objs: Sequence[Col],
        row_objs: Sequence[Row],
        qubit_objs: Sequence[Qubit],
        col_idx: Sequence[int],
        row_idx: Sequence[int],
        dropoff_qs: Sequence[int],
    ):
        for i in range(len(col_idx)):
            col_objs[col_idx[i]].active = False
        for i in range(len(row_idx)):
            row_objs[row_idx[i]].active = False
        for q_id in dropoff_qs:
            qubit_objs[q_id].array = 'SLM'

    def verify(
            self,
            col_objs: Sequence[Col],
            row_objs: Sequence[Row],
            qubit_objs: Sequence[Qubit],
            col_idx: Sequence[int],
            col_xs: Sequence[int],
            row_idx: Sequence[int],
            row_ys: Sequence[int],
            dropoff_qs: Sequence[int],
    ):
        a = len(col_idx)
        b = len(col_xs)
        if a != b:
            raise ValueError(
                f'{self.name}: col arguments invalid {a} idx, {b} xs')
        a = len(row_idx)
        b = len(row_ys)
        if a != b:
            raise ValueError(
                f'{self.name}: row arguments invalid {a} idx, {b} ys.')

        for i in range(len(col_idx)):
            if not col_objs[col_idx[i]].active:
                raise ValueError(
                    f'{self.name}: col {col_idx[i]} already dectivated.')
            for j in range(col_idx[i]):
                if (col_objs[j].active
                        and col_objs[j].x > col_xs[i] - AOD_SEP):
                    raise ValueError(
                        f'{self.name}: col {j} at x={col_objs[j].x} is '
                        f'too left for col {col_idx[i]} '
                        f'to deactivate at x={col_xs[i]}.')
            for j in range(col_idx[i]+1, len(col_objs)):
                if (col_objs[j].active
                        and col_objs[j].x - AOD_SEP < col_xs[i]):
                    raise ValueError(
                        f'{self.name}: col {j} at x={col_objs[j].x} is '
                        f'too right for col {col_idx[i]} '
                        f'to deactivate at x={col_xs[i]}.')
        for i in range(len(row_idx)):
            if not row_objs[row_idx[i]].active:
                raise ValueError(
                    f'{self.name}: row {row_idx[i]} already deactivated.')
            for j in range(row_idx[i]):
                if (row_objs[j].active
                        and row_objs[j].y > row_ys[i] - AOD_SEP):
                    raise ValueError(
                        f'{self.name}: row {j} at y={row_objs[j].y} is '
                        f'too high for row {row_idx[i]} '
                        f'to deactivate at y={row_ys[i]}.'
                    )
            for j in range(row_idx[i]+1, len(row_objs)):
                if (row_objs[j].active
                        and row_objs[j].y-AOD_SEP < row_ys[i]):
                    raise ValueError(
                        f'{self.name}: row {j} at y={col_objs[j].y} is '
                        f'too low for row {row_idx[i]} '
                        f'to deactivate at y={row_ys[i]}.'
                    )

        deactive_xys = []
        active_xs = [col.x for col in col_objs if col.active]
        for x in active_xs:
            for y in row_ys:
                deactive_xys.append((x, y))

        for q_id in range(len(qubit_objs)):
            if q_id in dropoff_qs:
                if (qubit_objs[q_id].x, qubit_objs[q_id].y) not in deactive_xys:
                    raise ValueError(
                        f'{self.name}: q {q_id} not dropped off from '
                        f'col {qubit_objs[q_id].c} row {qubit_objs[q_id].r} '
                        f'at x={qubit_objs[q_id].x} y={qubit_objs[q_id].y}.'
                    )
            elif qubit_objs[q_id].array == 'AOD':
                if (qubit_objs[q_id].x, qubit_objs[q_id].y) in deactive_xys:
                    raise ValueError(
                        f'{self.name}: q {q_id} wrongfully dropped off from '
                        f'col {qubit_objs[q_id].c} row {qubit_objs[q_id].r} '
                        f'at x={qubit_objs[q_id].x} y={qubit_objs[q_id].y}.'
                    )


class Rydberg(Inst):
    def __init__(self,
                 s: int,
                 col_objs: Sequence[Col],
                 row_objs: Sequence[Row],
                 qubit_objs: Sequence[Qubit],
                 gates: Sequence[Mapping[str, int]]):
        super().__init__('Rydberg', prefix=f'Rydberg_{s}', stage=s)
        self.verify(gates, qubit_objs)
        super().write_code(
            col_objs,
            row_objs,
            qubit_objs,
            {'gates': gates, 'duration': T_RYDBERG, }
        )

    def verify(
            self,
            gates: Sequence[Mapping[str, int]],
            qubit_objs: Sequence[Qubit]):
        # for g in gates:
        #     q0 = {'id': qubit_objs[g['q0']].id,
        #           'x': qubit_objs[g['q0']].x,
        #           'y': qubit_objs[g['q0']].y}
        #     q1 = {'id': qubit_objs[g['q1']].id,
        #           'x': qubit_objs[g['q1']].x,
        #           'y': qubit_objs[g['q1']].y}
        #     if (q0['x']-q1['x'])**2 + (q0['y']-q1['y'])**2 > R_B**2:
        #         raise ValueError(
        #             f"{self.name}: q{q0['id']} at x={q0['x']} y={q0['y']} "
        #             f"and q{q1['id']} at x={q1['x']} y={q1['y']} "
        #             f"are farther away than Rydberg range."
        #         )
        return


# class for big ops: ReloadRow, Reload, OffloadRow, Offload, SwapPair, Swap
# internally, these are lists of basic operations.

# todo: is there verification needed on the ComboInst lebvel?
class ComboInst:
    pass


class ComboInst():
    """class for combined instructions which is a sequence of combined
    instructions or DPQA instructions.
    """

    def __init__(
            self,
            type: str,
            prefix = None,
            suffix = None,
            stage: int = -1,
    ):
        """init method for combined instructions.

        Args:
            type (str): 
            prefix (str | None, optional): Defaults to None.
            suffix (str | None, optional): Defaults to None.
            stage (int, optional): Defaults to -1.
        """
        self.type = type
        self.name = (f'{prefix}:' if prefix else '') + \
            type + (f'_{suffix}' if suffix else '')
        self.stage = stage
        self.duration = -1
        self.insts = []

    def emit(self) -> Sequence[Mapping[str, Any]]:
        """combine the code of each Inst inside this ComboInst and return."""
        code = []
        for inst in self.insts:
            code += inst.emit()
        return code

    def emit_full(self) -> Sequence[Mapping[str, Any]]:
        code = []
        for inst in self.insts:
            code += inst.emit_full()
        return code

    def append_inst(self, inst):
        self.insts.append(inst)

    def prepend_inst(self, inst):
        self.insts.insert(0, inst)

    def is_trivial(self) -> bool:
        for inst in self.insts:
            if not inst.is_trivial():
                return False
        return True

    def remove_trivial_insts(self):
        nontrivial_insts = []
        for inst in self.insts:
            inst.remove_trivial_insts()
            if not inst.is_trivial():
                nontrivial_insts.append(inst)
        self.insts = nontrivial_insts


class ReloadRow(ComboInst):
    def __init__(self, s: int, r: int, prefix: str = ''):
        super().__init__('ReloadRow', prefix=prefix, suffix=str(r), stage=s)
        self.row_id = r
        self.moving_cols_id = []
        self.moving_cols_begin = []
        self.moving_cols_end = []

    def add_col_shift(self, id: int, begin: int, end: int):
        self.moving_cols_id.append(id)
        self.moving_cols_begin.append(begin)
        self.moving_cols_end.append(end)

    def generate_col_shift(
            self,
            col_objs: Sequence[Col],
            row_objs: Sequence[Row],
            qubit_objs: Sequence[Qubit],
    ):
        self.insts.append(
            Move(s=self.stage,
                 col_objs=col_objs,
                 row_objs=row_objs,
                 qubit_objs=qubit_objs,
                 col_idx=self.moving_cols_id,
                 col_begin=self.moving_cols_begin,
                 col_end=self.moving_cols_end,
                 row_idx=[],
                 row_begin=[],
                 row_end=[],
                 prefix=self.name + ':ColShift')
        )

    def generate_row_activate(
            self,
            col_objs: Sequence[Col],
            row_objs: Sequence[Row],
            qubit_objs: Sequence[Qubit],
            cols: Sequence[int],
            xs: Sequence[int],
            y: int,
            pickup_qs: Sequence[int],
    ):
        self.insts.append(
            Activate(s=self.stage,
                     col_objs=col_objs,
                     row_objs=row_objs,
                     qubit_objs=qubit_objs,
                     col_idx=cols,
                     col_xs=xs,
                     row_idx=[self.row_id, ],
                     row_ys=[y, ],
                     pickup_qs=pickup_qs,
                     prefix=self.name)
        )

    def generate_parking(
            self,
            col_objs: Sequence[Col],
            row_objs: Sequence[Row],
            qubit_objs: Sequence[Qubit],
            shift_down: int,
            col_idx: Sequence[int] = [],
            col_begin: Sequence[int] = [],
            col_end: Sequence[int] = []
    ):
        self.insts.append(
            Move(s=self.stage,
                 col_objs=col_objs,
                 row_objs=row_objs,
                 qubit_objs=qubit_objs,
                 col_idx=col_idx,
                 col_begin=col_begin,
                 col_end=col_end,
                 row_idx=[self.row_id, ],
                 row_begin=[row_objs[self.row_id].y],
                 row_end=[row_objs[self.row_id].y - shift_down],
                 prefix=self.name + ':Parking')
        )


class Reload(ComboInst):
    def __init__(self, s: int):
        super().__init__('Reload', suffix=str(s), stage=s)

    def add_row_reload(self, r: int):
        self.insts.append(ReloadRow(self.stage, r, prefix=self.name))
        return self.insts[-1]


class OffloadRow(ComboInst):
    def __init__(self, s: int, r: int, prefix: str = ''):
        super().__init__('OffloadRow', prefix=prefix, suffix=str(r), stage=s)
        self.row_id = r
        self.moving_cols_id = []
        self.moving_cols_begin = []
        self.moving_cols_end = []

    def add_col_shift(self, id: int, begin: int, end: int):
        self.moving_cols_id.append(id)
        self.moving_cols_begin.append(begin)
        self.moving_cols_end.append(end)

    def generate_col_shift(
            self,
            col_objs: Sequence[Col],
            row_objs: Sequence[Row],
            qubit_objs: Sequence[Qubit],
    ):
        self.insts.append(
            Move(s=self.stage,
                 col_objs=col_objs,
                 row_objs=row_objs,
                 qubit_objs=qubit_objs,
                 col_idx=self.moving_cols_id,
                 col_begin=self.moving_cols_begin,
                 col_end=self.moving_cols_end,
                 row_idx=[],
                 row_begin=[],
                 row_end=[],
                 prefix=self.name + ':ColShift')
        )

    def generate_row_shift(
            self,
            col_objs: Sequence[Col],
            row_objs: Sequence[Row],
            qubit_objs: Sequence[Qubit],
            site_y: int,
    ):
        self.insts.append(
            Move(s=self.stage,
                 col_objs=col_objs,
                 row_objs=row_objs,
                 qubit_objs=qubit_objs,
                 col_idx=[],
                 col_begin=[],
                 col_end=[],
                 row_idx=[self.row_id, ],
                 row_begin=[row_objs[self.row_id].y],
                 row_end=[site_y*Y_SITE_SEP, ],
                 prefix=self.name + ':RowDownShift')
        )

    def generate_row_deactivate(
            self,
            col_objs: Sequence[Col],
            row_objs: Sequence[Row],
            qubit_objs: Sequence[Qubit],
            dropoff_qs: Sequence[int],
    ):
        self.insts.append(
            Deactivate(s=self.stage,
                       col_objs=col_objs,
                       row_objs=row_objs,
                       qubit_objs=qubit_objs,
                       col_idx=[],
                       col_xs=[],
                       row_idx=[self.row_id, ],
                       row_ys=[row_objs[self.row_id].y, ],
                       dropoff_qs=dropoff_qs,
                       prefix=self.name)
        )


class Offload(ComboInst):
    def __init__(self, s: int):
        super().__init__('Offload', suffix=str(s), stage=s)

    def add_row_offload(self, r: int):
        self.insts.append(OffloadRow(self.stage, r, prefix=self.name))
        return self.insts[-1]

    def all_cols_deactivate(
            self,
            col_objs: Sequence[Col],
            row_objs: Sequence[Row],
            qubit_objs: Sequence[Qubit],
    ):
        self.insts.append(
            Deactivate(
                s=self.stage,
                col_objs=col_objs,
                row_objs=row_objs,
                qubit_objs=qubit_objs,
                col_idx=[c.id for c in col_objs if c.active],
                col_xs=[c.x for c in col_objs if c.active],
                row_idx=[],
                row_ys=[],
                dropoff_qs=[],
                prefix=self.name
            )
        )


class SwapPair(ComboInst):
    """swap a pair of atoms A and B at the same interaction site
       O stands for an empty trap

            A  B

        step1: activate row 0 and col 0 on A to pick it up
            |
         ---A--B---
            |

        step2: move A to below B
               |
            O  B
               |
           ----A-----
               |

        step3: activate row 1 to pick up B
               |
         ---O--B---
               |
         ------A---
               |

        step4: move A and B to the place A used to be
            |
        ----B--O--
            |
        ----A-----
            |

        step5: deactivate row 1 to drop off B
            |
            B  O
            |
        ----A-----
            |

        step6: move A to the place B used to be
               |
         ---B--A---
               |

        step7: deactivate row 0 and col 0 to drop off A

            B  A
    """

    def __init__(self,
                 s: int,
                 col_objs: Sequence[Col],
                 row_objs: Sequence[Row],
                 qubit_objs: Sequence[Qubit],
                 left_q_id: int,
                 right_q_id: int,
                 prefix: str = ''
                 ):
        super().__init__(
            'SwapPair',
            prefix=prefix,
            suffix=f'({left_q_id},{right_q_id})',
            stage=s
        )

        # precondition
        left_x = qubit_objs[left_q_id].x
        right_x = qubit_objs[right_q_id].x
        y = qubit_objs[left_q_id].y
        qubit_objs[left_q_id].c = 0
        qubit_objs[left_q_id].r = 0
        qubit_objs[right_q_id].c = 0
        qubit_objs[right_q_id].r = 1

        self.insts = [
            Activate(s,
                     col_objs=col_objs,
                     row_objs=row_objs,
                     qubit_objs=qubit_objs,
                     col_idx=[0, ],
                     col_xs=[left_x, ],
                     row_idx=[0, ],
                     row_ys=[y, ],
                     pickup_qs=[left_q_id, ],
                     prefix=self.name + f':PickUp_q{left_q_id}'),
            Move(s,
                 col_objs=col_objs,
                 row_objs=row_objs,
                 qubit_objs=qubit_objs,
                 col_idx=[0, ],
                 col_begin=[left_x, ],
                 col_end=[right_x, ],
                 row_idx=[0, ],
                 row_begin=[y, ],
                 row_end=[y-AOD_SEP, ],
                 prefix=self.name + f':tmp<-q@{left_x}'),
            Activate(s,
                     col_objs=col_objs,
                     row_objs=row_objs,
                     qubit_objs=qubit_objs,
                     col_idx=[],
                     col_xs=[],
                     row_idx=[1, ],
                     row_ys=[y, ],
                     pickup_qs=[right_q_id, ],
                     prefix=self.name + f':PickUp_{right_q_id}'),
            Move(s,
                 col_objs=col_objs,
                 row_objs=row_objs,
                 qubit_objs=qubit_objs,
                 col_idx=[0, ],
                 col_begin=[right_x, ],
                 col_end=[left_x, ],
                 row_idx=[],
                 row_begin=[],
                 row_end=[],
                 prefix=self.name + f':q@{left_x}<-q@{right_x}'),
            Deactivate(s,
                       col_objs=col_objs,
                       row_objs=row_objs,
                       qubit_objs=qubit_objs,
                       col_idx=[],
                       col_xs=[],
                       row_idx=[1, ],
                       row_ys=[y, ],
                       dropoff_qs=[right_q_id, ],
                       prefix=self.name + f':DropOff_q{right_q_id}'),
            Move(s,
                 col_objs=col_objs,
                 row_objs=row_objs,
                 qubit_objs=qubit_objs,
                 col_idx=[0, ],
                 col_begin=[left_x, ],
                 col_end=[right_x, ],
                 row_idx=[0, ],
                 row_begin=[y-AOD_SEP, ],
                 row_end=[y, ],
                 prefix=self.name + f':q@{right_x}<-tmp'),
            Deactivate(s,
                       col_objs=col_objs,
                       row_objs=row_objs,
                       qubit_objs=qubit_objs,
                       col_idx=[0, ],
                       col_xs=[right_x, ],
                       row_idx=[0, ],
                       row_ys=[y, ],
                       dropoff_qs=[left_q_id, ],
                       prefix=self.name + f':DropOff_q{left_q_id}'),
        ]


class Swap(ComboInst):
    def __init__(self, s: int):
        super().__init__('Swap', suffix=str(s), stage=s)

    def add_swap_pair(
            self,
            col_objs: Sequence[Col],
            row_objs: Sequence[Row],
            qubit_objs: Sequence[Qubit],
            q_id: int,
            qq_id: int
    ):
        left_q_id = q_id
        right_q_id = qq_id
        if qubit_objs[q_id].x > qubit_objs[qq_id].x:
            left_q_id = qq_id
            right_q_id = q_id
        self.insts.append(
            SwapPair(self.stage, col_objs, row_objs,
                     qubit_objs, left_q_id, right_q_id, self.name)
        )


class CodeGen():
    """Generate code files: json containing a list of dict, each one 
    corresponding to a DPQA instruction defined above. 
    """

    def __init__(
            self,
            data,
            no_transfer: bool = False,
    ):
        self.read_compiled(data)

    def read_compiled(self, data):
        self.n_q = data['n_q']
        self.x_high = int(data['n_x'])
        self.y_high = int(data['n_y'])
        self.c_high = int(data['n_c'])
        self.r_high = int(data['n_r'])
        self.layers = data['layers']
        self.n_t = len(self.layers)
        # print(self.n_t)

        """change of convention. In solve() and the SMT model, a/c/r_s govern
        the movement from stage s to stage s+1, i.e.,
            ----------- x/y_0 
            | a/c/r_0 |
            ----------- x/y_1
            | a/c/r_1 |
            ----------- x/y_2
            | a/c/r_2 |

        However, the initial stage is very special in codegen and animation,
        so we generate code in this way:
           Rydberg_0  (special)  <------- x/y_0
           
           Swap_1 (optional)  <---. 
           Reload_1  <----------\  \ 
           BigMove_1  <--------- a/c/r_1
           Offload_1  <---------/ 
           Rydberg_1  <----------- x/y_1

               ...

        Thus, the movement between stage s and s+1 should be govened by a/c/r
        with subscript s+1, e.g., BigMove_1 uses a/c/r_1 above.
        So we need to shift the a/c/r variable values here.
        """
        for i in range(self.n_t-1, 0, -1):
            for q in range(self.n_q):
                self.layers[i]['qubits'][q]['a'] =\
                    self.layers[i - 1]['qubits'][q]['a']
                self.layers[i]['qubits'][q]['c'] =\
                    self.layers[i - 1]['qubits'][q]['c']
                self.layers[i]['qubits'][q]['r'] =\
                    self.layers[i - 1]['qubits'][q]['r']

        # infer some info of the AODs
        self.aod_from_compiled()

    def aod_from_compiled(self):
        for s, layer in enumerate(self.layers):
            if s == 0:
                continue
            layer['row'] = [{'id': i, 'qs': []}
                            for i in range(self.r_high-0)]
            layer['col'] = [{'id': i, 'qs': []}
                            for i in range(self.c_high-0)]
            prev_layer = self.layers[s-1]

            # figure out in the movement from stage s-1 to s:
            # - before the movement, where is each row 'y_begin'
            # - after the movement, where is each row 'y_end'
            # - what qubits are in this row 'qs'
            # similar for each AOD column
            for i, q in enumerate(layer['qubits']):
                if layer['qubits'][i]['a']:
                    layer['row'][q['r']]['y_begin'] =\
                        prev_layer['qubits'][i]['y']
                    layer['row'][q['r']]['y_end'] = q['y']
                    layer['row'][q['r']]['qs'].append(q['id'])
                    layer['col'][q['c']]['x_begin'] =\
                        prev_layer['qubits'][i]['x']
                    layer['col'][q['c']]['x_end'] = q['x']
                    layer['col'][q['c']]['qs'].append(q['id'])

            # figure out in the movement from stage s-1 to s:
            # - before the movement, which columns have site coord = X
            # - for all these cols, what is the relevant order from left 'offset'
            # - after the movement, which columns have site coord = X
            # - for all these cols, what is the relevant order from left 'offset'
            # similar for the AOD rows
            for case in ['_begin', '_end']:
                x_cols = []
                for x in range(self.x_high):
                    cols_at_x = []
                    for c in range(self.c_high):
                        if (layer['col'][c]['qs']
                                and layer['col'][c]['x' + case] == x):
                            cols_at_x.append(c)
                    for i, c in enumerate(cols_at_x):
                        layer['col'][c]['offset' + case] = i
                    x_cols.append(cols_at_x)
                layer['x_cols' + case] = x_cols
                y_rows = []
                for y in range(self.y_high):
                    rows_at_y = []
                    for r in range(self.r_high):
                        if (layer['row'][r]['qs']
                                and layer['row'][r]['y' + case] == y):
                            rows_at_y.append(r)
                    for i, r in enumerate(rows_at_y):
                        layer['row'][r]['offset' + case] = i
                    y_rows.append(rows_at_y)
                layer['y_rows' + case] = y_rows

    def builder(self, no_transfer: bool):
        qubits = [Qubit(i) for i in range(self.n_q)]
        rows = [Row(i) for i in range(self.r_high)]
        cols = [Col(i) for i in range(self.c_high)]
        program = ComboInst('Program')


        # read to comment in read_compiled() for structure of this method.
        init = self.builder_init(cols, rows, qubits, program)  # has Rydberg_0
        # time_builder_swap = 0
        time_builder_reload = 0
        time_builder_move = 0
        time_builder_offload = 0
        time_builder_rydberg = 0
        for s in range(1, len(self.layers)):
            # print("s", s)
            # t_s = time.time()
            # self.builder_swap(s, cols, rows, qubits, program)
            # time_builder_swap += (time.time() - t_s)
            if (not no_transfer) or s == 1:
                # if we know there is not atom transfer, we can simply skip the
                # reload and offload procedures. However, we keep the first one
                # just for convenience of generating animations.
                t_s = time.time()
                self.builder_reload(s, cols, rows, qubits, program)
                time_builder_reload += (time.time() - t_s)
            # t_s = time.time()
            # print("builder_move")
            self.builder_move(s, cols, rows, qubits, program)
            # time_builder_move += (time.time() - t_s)
            if not no_transfer:
                # print("builder_offload")
                # t_s = time.time()
                self.builder_offload(s, cols, rows, qubits, program)
                # time_builder_offload += (time.time() - t_s)
            
            if len(self.layers[s]['gates']) > 0:
                # print(self.layers[s]['gates'])
                # t_s = time.time()
                self.builder_rydberg(s, cols, rows, qubits, program, init)
                # time_builder_rydberg += (time.time() - t_s)
            # data = program.emit_full()
            # import json
            # with open("test_s{}.json".format(s), 'w') as f:  
            #     json.dump(program.emit_full(), f)
        
        # print("time-builder_swap: {}".format(time_builder_swap))
        # print("time-builder_move: {}".format(time_builder_move))
        # print("time-builder_reload: {}".format(time_builder_reload))
        # print("time-builder_offload: {}".format(time_builder_offload))
        # print("time-builder_rydberg: {}".format(time_builder_rydberg))
        # t_s = time.time()
        program.remove_trivial_insts()    
        # print("time-remove_trivial_insts: {}".format(time.time() - t_s))
        # t_s = time.time()
        program.prepend_inst(init)
        # print("time-prepend_inst: {}".format(time.time() - t_s))
        
        return program

    def builder_init(
            self,
            cols: Sequence[Col],
            rows: Sequence[Row],
            qubits: Sequence[Qubit],
            program: ComboInst,
    ) -> Inst:
        slm_qubit_idx = list(range(self.n_q))  # put all qubits in SLM
        slm_qubit_xys = [(
            X_SITE_SEP * self.layers[0]['qubits'][i]['x'],
            Y_SITE_SEP * self.layers[0]['qubits'][i]['y']
        ) for i in range(self.n_q)]  # put all qubits in the left trap

        # when there are more than one qubit in a site at the beginning,
        # need to put one of them in the right trap.
        for g in self.layers[0]['gates']:
            a0 = self.layers[1]['qubits'][g['q0']]['a']
            a1 = self.layers[1]['qubits'][g['q1']]['a']
            x_left = X_SITE_SEP * self.layers[0]['qubits'][g['q0']]['x']
            x_right = x_left + SITE_WIDTH
            y = Y_SITE_SEP * self.layers[0]['qubits'][g['q0']]['y']

            # if both atoms are in AOD, use their column indices to decide
            # which one to put in the left trap and which one to the right
            # if they have the same col index, the order does not matter
            # in Reload, we will pick them up in different rows.
            if (a0 == 1
                and a1 == 1
                and self.layers[1]['qubits'][g['q0']]['c'] >
                    self.layers[1]['qubits'][g['q1']]['c']):
                slm_qubit_xys[g['q0']] = (x_right, y)
            else:
                slm_qubit_xys[g['q1']] = (x_right, y)

        init = Init(cols, rows, qubits,
                    slm_qubit_idx=slm_qubit_idx,
                    slm_qubit_xys=slm_qubit_xys,
                    data={
                        'n_q': self.n_q,
                        'x_high': self.x_high,
                        'y_high': self.y_high,
                        'c_high': self.c_high,
                        'r_high': self.r_high,
                    })

        self.builder_rydberg(0, cols, rows, qubits, program, init)

        return init

    def builder_rydberg(
            self,
            s: int,
            cols: Sequence[Col],
            rows: Sequence[Row],
            qubits: Sequence[Qubit],
            program: ComboInst,
            init: Inst
    ):
        program.append_inst(
            Rydberg(s, cols, rows, qubits, self.layers[s]['gates'])
        )
        init.add_slms([(q.x, q.y) for q in qubits])

    def builder_swap(
            self,
            s: int,
            cols: Sequence[Col],
            rows: Sequence[Row],
            qubits: Sequence[Qubit],
            program: ComboInst,
    ):
        swap_obj = Swap(s)
        prev_layer = self.layers[s-1]
        this_layer = self.layers[s]
        for q0_id in range(self.n_q):
            for q1_id in range(q0_id+1, self.n_q):
                q0_a = this_layer['qubits'][q0_id]['a']
                q1_a = this_layer['qubits'][q1_id]['a']
                if q0_a == 1 and q1_a == 1:
                    q0_x = prev_layer['qubits'][q0_id]['x']
                    q1_x = prev_layer['qubits'][q1_id]['x']
                    q0_y = prev_layer['qubits'][q0_id]['y']
                    q1_y = prev_layer['qubits'][q1_id]['y']
                    q0_c = this_layer['qubits'][q0_id]['c']
                    q1_c = this_layer['qubits'][q1_id]['c']
                    q0_r = this_layer['qubits'][q0_id]['r']
                    q1_r = this_layer['qubits'][q1_id]['r']
                    # if two qubits are at the same site and
                    # both being picked up in the same row
                    if (q0_x, q0_y, q0_r) == (q1_x, q1_y, q1_r):
                        # if their position and col indeces are in reverse order
                        if (
                            (q0_c > q1_c
                             and qubits[q0_id].x < qubits[q1_id].x
                             )
                            or (q0_c < q1_c
                                and qubits[q0_id].x > qubits[q1_id].x
                                )
                        ):
                            swap_obj.add_swap_pair(cols, rows, qubits,
                                                   q0_id, q1_id)
        program.append_inst(swap_obj)

    def builder_move(
            self,
            s: int,
            cols: Sequence[Col],
            rows: Sequence[Row],
            qubits: Sequence[Qubit],
            program: ComboInst,
    ):
        """to avoid collision, the big moves will end in a slightly adjusted
        location: 1 um to +X direction and AOD_SEP(2) um to +Y direction
        Suppose the O's are the two traps in an interaction site
                O--------O
        The atoms will finish the big move in
                O--------O
             -----A-----
                  -----B-----
        (+X is ->, +Y is down, 1um is --, and AOD_SEP is one line height)
        There can be other cases, e.g., A and B are in the same column
                O--------O
             -----A-----
             -----B-----
        Or, A and B are in the same row
                O--------O
             -----A---B------
        """

        col_idx = []
        col_begin = []
        col_end = []
        for col_id in range(self.c_high):
            if cols[col_id].active:
                col_idx.append(col_id)
                col_begin.append(cols[col_id].x)
                site_x = self.layers[s]['col'][col_id]['x_end']
                offset = self.layers[s]['col'][col_id]['offset_end']
                col_end.append(1 + site_x * X_SITE_SEP + AOD_SEP * offset)

        row_idx = []
        row_begin = []
        row_end = []
        for row_id in range(self.r_high):
            if rows[row_id].active:
                row_idx.append(row_id)
                row_begin.append(rows[row_id].y)
                site_y = self.layers[s]['row'][row_id]['y_end']
                offset = self.layers[s]['row'][row_id]['offset_end']
                row_end.append(site_y * Y_SITE_SEP + AOD_SEP * (1 + offset))

        program.append_inst(
            Move(s=s,
                 col_objs=cols,
                 row_objs=rows,
                 qubit_objs=qubits,
                 col_idx=col_idx,
                 col_begin=col_begin,
                 col_end=col_end,
                 row_idx=row_idx,
                 row_begin=row_begin,
                 row_end=row_end,
                 prefix=f'BigMove_{s}')
        )

    def builder_reload(
            self,
            s: int,
            cols: Sequence[Col],
            rows: Sequence[Row],
            qubits: Sequence[Qubit],
            program: ComboInst,
    ):
        layer = self.layers[s]
        prev_layer = self.layers[s-1]
        reload_obj = Reload(s)
        # reload row by row
        dict_site_qs = {}
        for q_id in range(len(qubits)):
            # find out the qubits with site_x and in row_id
            if (layer['qubits'][q_id]['a'] == 1):
                site_x = prev_layer['qubits'][q_id]['x']
                row_id = layer['qubits'][q_id]['r']
                
                qubits[q_id].r = row_id
                qubits[q_id].c = layer['qubits'][q_id]['c']

                if row_id in dict_site_qs.keys():
                    if site_x in dict_site_qs[row_id]:
                        dict_site_qs[row_id][site_x].append(q_id)
                    else:
                        dict_site_qs[row_id][site_x] = [q_id]
                else:
                    dict_site_qs[row_id] = {site_x: [q_id]}

        for row_id in sorted(dict_site_qs):
            reloadRow_obj = reload_obj.add_row_reload(row_id)
            pickup_qs = []
            cols_to_active = []
            x_to_activate = []
            # consider the movements in a row of sites
            for site_x in sorted(dict_site_qs[row_id]):
                site_qs = dict_site_qs[row_id][site_x]
                # shift the 1 or 2 cols that are picking up qubits
                if len(site_qs) == 2:
                    # which qubit is on the left and which on the right
                    [q_id_left, q_id_right] = site_qs
                    if layer['qubits'][q_id_left]['c'] >\
                            layer['qubits'][q_id_right]['c']:
                        tmp = q_id_left
                        q_id_left = q_id_right
                        q_id_right = tmp

                    # current location of Cols
                    col_id_left = layer['qubits'][q_id_left]['c']
                    col_id_right = layer['qubits'][q_id_right]['c']
                    lower_offset =\
                        layer['col'][col_id_left]['offset_begin']
                    upper_offset =\
                        layer['col'][col_id_right]['offset_begin']

                    # target locations of Cols
                    lower_x = qubits[q_id_left].x
                    upper_x = qubits[q_id_right].x

                    # process the Col on the left
                    if not cols[col_id_left].active:
                        cols_to_active.append(col_id_left)
                        x_to_activate.append(lower_x)
                    else:
                        reloadRow_obj.add_col_shift(
                            id=col_id_left,
                            begin=cols[col_id_left].x,
                            end=lower_x)

                    # process the Col on the right
                    if not cols[col_id_right].active:
                        cols_to_active.append(col_id_right)
                        x_to_activate.append(upper_x)
                    else:
                        reloadRow_obj.add_col_shift(
                            id=col_id_right,
                            begin=cols[col_id_right].x,
                            end=upper_x)

                elif len(site_qs) == 1:
                    # for convience later on, we still keep the *_left and
                    # *_right vars even if there is one qubit to pick up
                    q_id = site_qs[0]
                    col_id_left = layer['qubits'][q_id]['c']
                    col_id_right = col_id_left
                    lower_offset =\
                        layer['col'][col_id_left]['offset_begin']
                    upper_offset = lower_offset
                    lower_x = qubits[q_id].x
                    upper_x = lower_x
                    if not cols[col_id_left].active:
                        cols_to_active.append(col_id_left)
                        x_to_activate.append(lower_x)
                    else:  # col already active, shift it to align with q_id
                        reloadRow_obj.add_col_shift(
                            id=col_id_left,
                            begin=cols[col_id_left].x,
                            end=lower_x)

                elif len(site_qs) > 2:
                    raise ValueError(
                        f"builder reload {s} row {row_id} site {site_x}:"
                        f" more than 2 qubits"
                    )
                else:
                    continue

                # shift other Cols that are already activated. Those Cols
                # may not be picking up any qubit in this row, but due to
                # AOD order constraints, we need to shift them to have the
                # correct order when we shift the Cols that indeed are
                # picking up some qubit in this row.
                for col_id in layer['x_cols_begin'][site_x]:
                    # the Cols at site_x in the beginning between stage s-1
                    # and s (since we are processing Reload_s now)

                    if (cols[col_id].active
                        and col_id != col_id_left
                            and col_id != col_id_right):
                        # if this Col is not the one involved in loading

                        # if there is a col on the right of the cols for loading
                        if layer['col'][col_id]['offset_begin'] >\
                                upper_offset:
                            if layer['col'][col_id]['offset_begin'] - upper_offset == 1:
                                reloadRow_obj.add_col_shift(
                                id=col_id,
                                begin=cols[col_id].x,
                                end=upper_x + AOD_SEP)
                            elif layer['col'][col_id]['offset_begin'] - upper_offset == 2:
                                reloadRow_obj.add_col_shift(
                                id=col_id,
                                begin=cols[col_id].x,
                                end=upper_x + AOD_SEP * 3)
                            

                        # if there is a col on the left of the cols for loading
                        elif layer['col'][col_id]['offset_begin'] <\
                                lower_offset:
                            if layer['col'][col_id]['offset_begin'] - lower_offset == -1:
                                reloadRow_obj.add_col_shift(
                                id=col_id,
                                begin=cols[col_id].x,
                                end=lower_x - AOD_SEP)
                            elif layer['col'][col_id]['offset_begin'] - lower_offset == -2:
                                reloadRow_obj.add_col_shift(
                                id=col_id,
                                begin=cols[col_id].x,
                                end=lower_x - 3 * AOD_SEP)
                        # if there is a col in the middle of the cols for loading
                        else:
                            reloadRow_obj.add_col_shift(
                                id=col_id,
                                begin=cols[col_id].x,
                                end=lower_x + AOD_SEP *
                                (layer['col'][col_id]['offset_begin'] -
                                    lower_offset))

                pickup_qs += site_qs

            # print("pickup_qs")
            # print(pickup_qs)
            # input()
            # apply all the col shifts added previously
            reloadRow_obj.generate_col_shift(cols, rows, qubits)
            reloadRow_obj.generate_row_activate(
                cols,
                rows,
                qubits,
                cols_to_active,
                x_to_activate,
                layer['row'][row_id]['y_begin']*Y_SITE_SEP,
                pickup_qs
            )

            # shift down the finished row because later on, some other row
            # may need to adjust the Cols again and if we keep this row
            # as is, some qubits may collid into each other since at each
            # site, the y of two SLM traps are the same.
            num_rows = len(layer['y_rows_begin']
                            [layer['row'][row_id]['y_begin']])
            shift_down = (
                num_rows - layer['row'][row_id]['offset_begin'])*AOD_SEP
            # Also shift the Cols by 1 to avoid collisions.
            col_idx = [col_id for col_id,
                        col in enumerate(cols) if col.active]
            col_begin = [cols[col_id].x for col_id in col_idx]
            tmp_shift = {0: -AOD_SEP, 1: AOD_SEP, 2: 3*AOD_SEP}
            col_end = [layer['col'][col_id]['x_begin'] * X_SITE_SEP +
                        tmp_shift[layer['col'][col_id]['offset_begin']]
                        for col_id in col_idx]
            reloadRow_obj.generate_parking(
                cols,
                rows,
                qubits,
                shift_down,
                col_idx,
                col_begin,
                col_end
            )

        program.append_inst(reload_obj)

    def builder_offload(
            self,
            s: int,
            cols: Sequence[Col],
            rows: Sequence[Row],
            qubits: Sequence[Qubit],
            program: ComboInst,
    ):
        offload_obj = Offload(s)
        layer = self.layers[s]
        # the row-by-row processing is quite similar to Reload
        for row_id in range(self.r_high):
            if rows[row_id].active:
                dropoff_qs = []
                offloadRow_obj = offload_obj.add_row_offload(row_id)
                dict_site_q_slm = {}
                dict_site_q_aod = {}
                for q_id, q in enumerate(layer['qubits']):
                    if q['y'] == layer['row'][row_id]['y_end']:
                        site_x = q['x'] 
                        if (qubits[q_id].array == 'AOD'
                                and q['r'] == row_id):
                            dropoff_qs.append(q_id)
                            if site_x in dict_site_q_aod.keys():
                                dict_site_q_aod[site_x].append(q_id)
                            else:
                                dict_site_q_slm[site_x] = []
                                dict_site_q_aod[site_x] = [q_id]
                        if qubits[q_id].array == 'SLM':
                            if site_x in dict_site_q_slm.keys():
                                dict_site_q_slm[site_x].append(q_id)
                            else:
                                dict_site_q_aod[site_x] = []
                                dict_site_q_slm[site_x] = [q_id]
                for site_x in sorted(dict_site_q_slm):
                    site_q_slm = dict_site_q_slm[site_x]
                    site_q_aod = dict_site_q_aod[site_x]
                # for site_x in range(self.x_high):
                #     site_q_slm = []
                #     site_q_aod = []
                #     for q_id, q in enumerate(layer['qubits']):
                #         if (
                #             q['x'], q['y']
                #         ) == (
                #             site_x, layer['row'][row_id]['y_end']
                #         ):
                #             if (qubits[q_id].array == 'AOD'
                #                     and q['r'] == row_id):
                #                 dropoff_qs.append(q_id)
                #                 site_q_aod.append(q_id)
                #             if qubits[q_id].array == 'SLM':
                #                 site_q_slm.append(q_id)
                    if len(site_q_aod) == 2:
                        [q_id_left, q_id_right] = site_q_aod
                        if layer['qubits'][q_id_left]['c'] >\
                                layer['qubits'][q_id_right]['c']:
                            tmp = q_id_left
                            q_id_left = q_id_right
                            q_id_right = tmp
                        col_id_left = layer['qubits'][q_id_left]['c']
                        col_id_right = layer['qubits'][q_id_right]['c']
                        lower_offset = layer['col'][col_id_left]['offset_end']
                        upper_offset = layer['col'][col_id_right]['offset_end']
                        lower_x = X_SITE_SEP*site_x
                        upper_x = X_SITE_SEP*site_x + SITE_WIDTH

                        offloadRow_obj.add_col_shift(
                            id=col_id_left,
                            begin=qubits[q_id_left].x,
                            end=lower_x
                        )
                        offloadRow_obj.add_col_shift(
                            id=col_id_right,
                            begin=qubits[q_id_right].x,
                            end=upper_x
                        )
                    elif len(site_q_aod) == 1:
                        q_id = site_q_aod[0]
                        col_id_left = layer['qubits'][q_id]['c']
                        col_id_right = col_id_left
                        lower_offset = layer['col'][col_id_left]['offset_end']
                        upper_offset = lower_offset
                        lower_x = X_SITE_SEP*site_x
                        if site_q_slm:
                            lower_x = 2*X_SITE_SEP*site_x + \
                                SITE_WIDTH - qubits[site_q_slm[0]].x
                        upper_x = lower_x
                        offloadRow_obj.add_col_shift(
                            id=col_id_left,
                            begin=qubits[q_id].x,
                            end=lower_x
                        )
                    elif len(site_q_aod) > 2:
                        raise ValueError(
                            f"builder offload {s} row {row_id} site {site_x}:"
                            f" more than 2 qubits"
                        )
                    else:
                        continue

                    for col_id in layer['x_cols_end'][site_x]:
                        if (cols[col_id].active
                            and col_id != col_id_left
                                and col_id != col_id_right):
                            if layer['col'][col_id]['offset_end'] >\
                                    upper_offset:
                                offloadRow_obj.add_col_shift(
                                    id=col_id,
                                    begin=cols[col_id].x,
                                    end=upper_x + AOD_SEP *
                                    (layer['col'][col_id]['offset_end'] -
                                     upper_offset) + 1
                                )
                            elif layer['col'][col_id]['offset_end'] <\
                                    lower_offset:
                                offloadRow_obj.add_col_shift(
                                    id=col_id,
                                    begin=cols[col_id].x,
                                    end=lower_x + AOD_SEP *
                                    (layer['col'][col_id]['offset_end'] -
                                     lower_offset) - 1
                                )
                            else:
                                offloadRow_obj.add_col_shift(
                                    id=col_id,
                                    begin=cols[col_id].x,
                                    end=lower_x + AOD_SEP *
                                    (layer['col'][col_id]['offset_end'] -
                                     lower_offset)
                                )

                # align the Cols to the correct locatios
                offloadRow_obj.generate_col_shift(cols, rows, qubits)
                # the rows are at the parked location finishing a BigMove,
                # so we need to shift them to align with the SLM traps.
                offloadRow_obj.generate_row_shift(
                    cols,
                    rows,
                    qubits,
                    layer['row'][row_id]['y_end']
                )
                offloadRow_obj.generate_row_deactivate(
                    cols, rows, qubits, dropoff_qs)
        offload_obj.all_cols_deactivate(cols, rows, qubits)
        program.append_inst(offload_obj)


def set_hardware_paramters(param: dict):
    R_B = param["R_B"]  # rydberg range
    AOD_SEP = param["AOD_SEP"]  # min AOD separation
    RYD_SEP = param["RYD_SEP"]  # sufficient distance to avoid Rydberg
    SITE_SLMS = param["SITE_SLMS"]  # number of SLMs in a site
    SITE_WIDTH = param["SITE_WIDTH"]  # total width of SLMs in a site
    SLM_SEP = AOD_SEP  # separation of SLMs inside a site
    X_SITE_SEP = RYD_SEP + SITE_WIDTH  # separation of sites in X direction
    Y_SITE_SEP = RYD_SEP  # separation of sites in Y direction

    # padding of the figure
    X_LOW_PAD = 2 * AOD_SEP
    Y_LOW_PAD = 4 * AOD_SEP
    X_HIGH_PAD = 2 * AOD_SEP
    Y_HIGH_PAD = 4 * AOD_SEP

    T_RYDBERG = param["T_RYDBERG"]  # microseconds for Rydberg
    T_ACTIVATE = param["T_ACTIVATE"]  # microseconds for (de)activating AOD
