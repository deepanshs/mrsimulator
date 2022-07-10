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
            "as the attribute name or use the generic `Method` class."
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
            "objects or use the generic `Method` class."
        )
        super().__init__(message)


class MixedSpectralDimensionTypeError(Exception):
    """Exception raised when SpectralDimension objects and dictionaries are used in the
    same spectral_dimension list

    Attributes:
        message: explanation of the error.
    """

    def __init__(self, message=None):
        message = message or (
            "Both dict and SpectralDimension objects found in spectral dimension list. "
            "Use either dict representation or SpectralDimension objects."
        )
        super().__init__(message)


class MissingSpectralDimensionError(Exception):
    """Exception raised when a generic Method object is missing the spectral_dimension
    argument or is passed no spectral dimensions.

    Attributes:
        message: explanation of the error.
    """

    def __init__(self, message=None):
        message = message or (
            "Method requires at least one SpectralDimension, none found."
        )
        super().__init__(message)


class MissingSpectralEventError(Exception):
    """Exception raised for missing SpectralEvent from SpectralDimension.

    Attributes:
        message: explanation of the error.
    """

    def __init__(self, message=None):
        message = message or (
            "SpectralDimension requires at least one SpectralEvent, none found."
        )
        super().__init__(message)


class FileConversionError(Exception):
    """Exception raised when after another error has been raised when attempting
    to parse an older file/dict to a compatible structure.

    Attributes:
        message: explanation of the error.
    """

    def __init__(self, message=None):
        message = message or (
            "Unable to convert the requested mrsim file/dict to a compatible "
            "structure. See the documentation at "
            "https://mrsimulator.readthedocs.io/en/stable/ to find out more."
        )
        super().__init__(message)
