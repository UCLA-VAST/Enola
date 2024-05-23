# coding: utf-8

from typing import Sequence, Tuple
from enola.scheduler.sequence_view import SequenceView

NodeType = int
ColorType = int
IndexType = int
LengthType = int
EdgeType = Tuple[NodeType, NodeType]
FanType = Sequence[NodeType]
FanViewType = SequenceView
## Stuff to add in future: MaximalFanType, FreeColorType, IncidentColorType


## Who cares ...
# SENTINEL_COLOR = (
# 0
# )  # this should be inside ColorType, if we want to use types correctly.
