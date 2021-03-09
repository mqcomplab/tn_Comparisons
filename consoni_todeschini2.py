from math import log
from base_tn import BaseTNComparisons


class ConsoniTodeschini2(BaseTNComparisons):
    """Class to calculate the Consoni-Todeschini 2 index.

    n=2 formula:
        (log(1 + p) - log(1 + b + c))/log(1 + p)

    Attributes
    ----------
    fingerprints : {np.ndarray, int}
        Numpy array with the fingerprints that will be compared.
        The fingerprints must be also given as Numpy arrays.
        If an int is given it assumes that one is comparing
        n random fingerprints of infinite length.
    c_threshold : {None, 'dissimilar', int}
        Coincidence threshold.
    characters : list
        List with the possible characters.
    n_characters : int
        Number of different characters.

    Properties
    ----------
    n_fingerprints : int
        Number of fingerprints that will be compared.

    Methods
    -------
    __init__(self, fingerprints, c_threshold=None, w_factor="fraction")
        Initialize the object.
    assign_characters(characters)
        Assign characters.
    assign_fingerprints(fingerprints)
        Assign fingerprints.
    assign_c_threshold(c_threshold)
        Assign coincidence threshold.
    matches()
        Calculate the matches between the fingerprints.
    set_general_coincidences()
        Calculate the general coincidence matrix.
    set_reduced_coincidences()
        Calculate the reduced coincidence vector.
    set_matches()
        Calculate the matches between the fingerprints.
    set_d_vector()
        Calculate the d vector.
    set_w_factor(w_factor)
        Calculate weight factors.
    set_weighted_matches()
        Calculate weighted matches.
    set_counters()
        Calculate the counters.
    set_p()
        Calculate p.
    set_weighted_p()
        Calculate weighted p.
    ct2_1sim_wdis()
        Calculate the index with 1-sim-counters and with weighted denominator.
    ct2_1sim_dis()
         Calculate the index with 1-sim-counters and with unweighted denominator.
    """

    def __init__(self, fingerprints, characters=[0, 1], c_threshold=None, w_factor="fraction"):
        """Initialize the object.

        Parameters
        ----------
        fingerprints : np.ndrarray
            Numpy array with the fingerprints that will be compared.
            The fingerprints must be also given as Numpy arrays.
        c_threshold : {None, 'dissimilar', int}
            Coincidence threshold.
        w_factor : {"fraction", "power_n"}
            Type of weight function that will be used.
        """
        super().__init__(fingerprints, characters, c_threshold, w_factor)
        self.ct2_1sim_wdis()
        self.ct2_1sim_dis()

    def ct2_1sim_wdis(self):
        """Calculate the index with 1-sim-counters and with weighted denominator.

        Note
        ----
        (log(1 + w_p) - log(1 + w_b + w_c))/log(1 + w_p)
        """
        numerator = log(1 + self.w_p) - log(1 + self.total_w_dis)
        denominator = log(1 + self.w_p)
        self.CT2_1sim_wdis = numerator/denominator

    def ct2_1sim_dis(self):
        """Calculate the index with 1-sim-counters and with unweighted denominator.

        Note
        ----
        (log(1 + w_p) - log(1 + w_b + w_c))/log(1 + p)
        """
        numerator = log(1 + self.w_p) - log(1 + self.total_w_dis)
        denominator = log(1 + self.p)
        self.CT2_1sim_dis = numerator/denominator
