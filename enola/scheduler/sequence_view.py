# coding: utf-8

from collections.abc import Sequence


# from https://more-itertools.readthedocs.io/en/stable/_modules/more_itertools/more.html#SequenceView
class SequenceView(Sequence):
    """Return a read-only view of the sequence object *target*.

    :class:`SequenceView` objects are analogous to Python's built-in
    "dictionary view" types. They provide a dynamic view of a sequence's items,
    meaning that when the sequence updates, so does the view.

        >>> seq = ['0', '1', '2']
        >>> view = SequenceView(seq)
        >>> view
        SequenceView(['0', '1', '2'])
        >>> seq.append('3')
        >>> view
        SequenceView(['0', '1', '2', '3'])

    Sequence views support indexing, slicing, and length queries. They act
    like the underlying sequence, except they don't allow assignment:

        >>> view[1]
        '1'
        >>> view[1:-1]
        ['1', '2']
        >>> len(view)
        4

    Sequence views are useful as an alternative to copying, as they don't
    require (much) extra storage.

    """

    def __init__(self, target):
        if not isinstance(target, Sequence):
            raise TypeError
        self._target = target

    def __getitem__(self, index):
        return self._target[index]

    def __len__(self):
        return len(self._target)

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__, repr(self._target))
