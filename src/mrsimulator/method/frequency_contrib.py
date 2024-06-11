from enum import Enum

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"

FREQ_LIST_ALL = [
    "Shielding1_0",
    "Shielding1_2",
    "Quad1_2",
    "Quad2_0",
    "Quad2_2",
    "Quad2_4",
    "J1_0",
    "J1_2",
    "D1_2",
    "Quad_Shielding_cross_0",
    "Quad_Shielding_cross_2",
    "Quad_Shielding_cross_4",
    "Quad_J_cross_0",
    "Quad_J_cross_2",
    "Quad_J_cross_4",
    "Quad_Dipolar_cross_0",
    "Quad_Dipolar_cross_2",
    "Quad_Dipolar_cross_4",
]

FREQ_ENUM_SHORTCUT = {
    "Shielding": {"Shielding1_0", "Shielding1_2"},
    "Isotropic": {"Shielding1_0", "J1_0"},
    "Quad": {"Quad1_2", "Quad2_0", "Quad2_2", "Quad2_4"},
    "J": {"J1_0", "J1_2"},
    "D": {"D1_2"},
    "cross": {
        "Quad_Shielding_cross_0",
        "Quad_Shielding_cross_2",
        "Quad_Shielding_cross_4",
        "Quad_J_cross_0",
        "Quad_J_cross_2",
        "Quad_J_cross_4",
        "Quad_Dipolar_cross_0",
        "Quad_Dipolar_cross_2",
        "Quad_Dipolar_cross_4",
    },
    "Quad_Shielding_cross": {
        "Quad_Shielding_cross_0",
        "Quad_Shielding_cross_2",
        "Quad_Shielding_cross_4",
    },
    "Quad_J_cross": {"Quad_J_cross_0", "Quad_J_cross_2", "Quad_J_cross_4"},
    "Quad_D_cross": {
        "Quad_Dipolar_cross_0",
        "Quad_Dipolar_cross_2",
        "Quad_Dipolar_cross_4",
    },
    "First_order": {
        "Shielding1_0",
        "Shielding1_2",
        "Quad1_2",
        "J1_0",
        "J1_2",
        "D1_2",
    },
    "Second_order": {
        "Quad2_0",
        "Quad2_2",
        "Quad2_4",
    },
    "Zeroth_rank": {
        "Shielding1_0",
        "Quad2_0",
        "J1_0",
        "Quad_Shielding_cross_0",
        "Quad_J_cross_0",
        "Quad_Dipolar_cross_0",
    },
    "Second_rank": {
        "Shielding1_2",
        "Quad1_2",
        "Quad2_2",
        "Quad_Shielding_cross_2",
        "Quad_J_cross_2",
        "Quad_Dipolar_cross_2",
    },
    "Fourth_rank": {
        "Quad2_4",
        "Quad_Shielding_cross_4",
        "Quad_J_cross_4",
        "Quad_Dipolar_cross_4",
    },
}


class FrequencyEnum(str, Enum):
    """Enumeration for selecting specific frequency contributions. The enumerations
    are:

    Attributes
    ----------

    Shielding1_0:
        Selects first-order and zeroth-rank nuclear shielding frequency contributions.

    Shielding1_2:
        Selects first-order and second-rank nuclear shielding frequency contributions.

    Quad1_2:
        Selects first-order and second-rank quadrupolar frequency contributions.

    Quad2_0:
        Selects second-order and zeroth-rank quadrupolar frequency contributions.

    Quad2_2:
        Selects second-order and second-rank quadrupolar frequency contributions.

    Quad2_4:
        Selects second-order and fourth-rank quadrupolar frequency contributions.

    J1_0:
        Selects first-order and zeroth-rank weak J-coupling frequency contributions.

    J1_2:
        Selects first-order and second-rank weak J-coupling frequency contributions.

    D1_2:
        Selects first-order and second-rank weak dipole frequency contributions.

    Quad_Shielding_cross_0:
        Selects zeroth-rank quad-shielding cross interaction.

    Quad_Shielding_cross_2:
        Selects second-rank quad-shielding cross-interaction.

    Quad_Shielding_cross_4:
        Selects fourth-rank quad-shielding cross-interaction.

    Quad_J_cross_0:
        Selects zeroth-rank quad-J-coupling cross-interaction.

    Quad_J_cross_2:
        Selects second-rank quad-J-coupling cross-interaction.

    Quad_J_cross_4:
        Selects fourth-rank quad-J-coupling cross-interaction.

    Quad_Dipolar_cross_0:
        Selects zeroth-rank quad-dipolar coupling cross-interaction.

    Quad_Dipolar_cross_2:
        Selects second-rank quad-dipolar coupling cross-interaction.

    Quad_Dipolar_cross_4:
        Selects fourth-rank quad-dipolar coupling cross-interaction.


    There are also shortcuts for including/excluding sets of contributions together.
    Frequency contributions can be excluded by including an exclamation mark in front of
    the string, for example ``"!Shielding"`` excludes all shielding interactions. The
    allowed shortcuts are:

    Shortcuts
    ---------

    ``"Shielding"``:
        Selects all shielding interactions
    ``"Isotropic"``:
        Selects first-order zeroth-rank shielding and first-order zeroth-rank J coupling
        interactions
    ``"Quad"``:
        Selects all quadrupolar interactions
    ``"J"``:
        Selects all J coupling interactions
    ``"D"``:
        Selects all dipolar interactions
    ``"cross"``:
        Selects all cross-term interactions
    ``"Quad_Shielding_cross"``:
        Selects all quadrupolar-shielding cross-terms
    ``"Quad_J_cross"``:
        Selects all quadrupolar-J-coupling cross-terms
    ``"Quad_D_cross"``:
        Selects all quadrupolar-dipolar-coupling cross-terms
    ``"First_order"``:
        Selects all first-order interactions
    ``"Second_order"``:
        Selects all second-order interactions
    ``"Zeroth_rank"``:
        Selects all zeroth-rank interactions
    ``"Second_rank"``:
        Selects all second-rank interactions
    ``"Fourth_rank"``:
        Selects all fourth-rank interactions
    """

    # Shielding1: str = "Shielding1"
    Shielding1_0: str = FREQ_LIST_ALL[0]
    Shielding1_2: str = FREQ_LIST_ALL[1]

    # All Shielding
    # Shielding: set = {FREQ_LIST_ALL[0], FREQ_LIST_ALL[1]}

    # Quad1: str = "Quad1"
    Quad1_2: str = FREQ_LIST_ALL[2]

    # Quad2: str = "Quad2"
    Quad2_0: str = FREQ_LIST_ALL[3]
    Quad2_2: str = FREQ_LIST_ALL[4]
    Quad2_4: str = FREQ_LIST_ALL[5]

    # J1: str = "J1"
    J1_0: str = FREQ_LIST_ALL[6]
    J1_2: str = FREQ_LIST_ALL[7]

    # D1: str = "D1"
    D1_2: str = FREQ_LIST_ALL[8]

    # Quad_Shielding: str = "Quad_Shielding"
    Quad_Shielding_cross_0: str = FREQ_LIST_ALL[9]
    Quad_Shielding_cross_2: str = FREQ_LIST_ALL[10]
    Quad_Shielding_cross_4: str = FREQ_LIST_ALL[11]

    Quad_J_cross_0: str = FREQ_LIST_ALL[12]
    Quad_J_cross_2: str = FREQ_LIST_ALL[13]
    Quad_J_cross_4: str = FREQ_LIST_ALL[14]

    Quad_Dipolar_cross_0: str = FREQ_LIST_ALL[15]
    Quad_Dipolar_cross_2: str = FREQ_LIST_ALL[16]
    Quad_Dipolar_cross_4: str = FREQ_LIST_ALL[17]

    class Config:
        extra = "forbid"

    def json(self, **kwargs) -> str:
        """Parse the class object to a JSON compliant python dictionary object."""
        return self.value

    def index(self) -> int:
        """Get the index of enumeration relative to FREQ_LIST_ALL."""
        return FREQ_LIST_ALL.index(self.value)


default_freq_contrib = [FrequencyEnum(item) for item in FREQ_LIST_ALL]
