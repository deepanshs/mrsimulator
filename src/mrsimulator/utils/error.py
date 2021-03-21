# -*- coding: utf-8 -*-
class NamedMethodError(Exception):
    pass


class ImmutableEventError(Exception):
    """Exception raised for Events are changed in a Named Method object.

    Attributes:
        cls_name: class name.
        message: explanation of the error.
    """

    def __init__(self, cls_name, message=None):
        if message is None:
            message = (
                f"Event objects are immutable for {cls_name} class. Remove all event "
                "objects or use `Method1D`/`Method2D` class."
            )
        super().__init__(message)
