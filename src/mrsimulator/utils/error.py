# -*- coding: utf-8 -*-
class NamedMethodError(Exception):
    """Exception raised when NamedClassMethods are renamed.

    Attributes:
        name: renamed name.
        cls_name: class name.
        message: explanation of the error.
    """

    def __init__(self, name, cls_name, message=None):
        message = message or (
            f"`name={name} != classname={cls_name}`. Use the class with the same name "
            "as the attribute name or use `Method1D` or `Method2D` class."
        )
        super().__init__(message)


class ImmutableEventError(Exception):
    """Exception raised for Events are changed in a Named Method object.

    Attributes:
        cls_name: class name.
        message: explanation of the error.
    """

    def __init__(self, cls_name, message=None):
        message = message or (
            f"Event objects are immutable for {cls_name} class. Remove all event "
            "objects or use `Method1D`/`Method2D` class."
        )
        super().__init__(message)
