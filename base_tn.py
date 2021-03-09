import numpy as np
from math import ceil


class BaseTNComparisons(object):
    """Base class to compare arbitrary numbers of fingerprints.

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
    """

    def __init__(self, fingerprints, characters=[0, 1], c_threshold=None, w_factor="fraction"):
        """Initialize the object.

        Parameters
        ----------
        fingerprints : {np.ndarray, int}
            Numpy array with the fingerprints that will be compared.
            The fingerprints must be also given as Numpy arrays.
            If an int is given it assumes that one is comparing
            n random fingerprints of infinite length.
        c_threshold : {None, 'dissimilar', int}
            Coincidence threshold.
        w_factor : {"fraction", "power_n"}
            Type of weight function that will be used.
        """
        self.assign_characters(characters)
        self.assign_fingerprints(fingerprints)
        self.assign_c_threshold(c_threshold)
        self.set_general_coincidences()
        self.set_reduced_coincidences()
        self.set_matches()
        self.set_d_vector()
        self.set_w_factor(w_factor)
        self.set_weighted_matches()
        self.set_counters()
        self.set_p()
        self.set_weighted_p()

    @property
    def n_fingerprints(self):
        """Return number of fingerprints.

        Returns
        -------
        n_fingerprints : int
            Number of fingerprints that will be compared.

        Note: If fingerprints is an int this is taken as the number of fingerprints
              that will be compared.
        """
        if isinstance(self.fingerprints, int):
            return self.fingerprints
        else:
            return len(self.fingerprints)

    def assign_characters(self, characters):
        """Assign characters.

        Parameters
        ----------
        characters : list
            List with the possible characters.
        """
        self.characters = characters
        self.n_characters = len(self.characters)

    def assign_fingerprints(self, fingerprints):
        """Assign fingerprints.

        Parameters
        ----------
        fingerprints : {np.ndarray, int}
            Numpy array with the fingerprints that will be compared.
            The fingerprints must be also given as Numpy arrays.
            If an int is given it assumes that one is comparing
            n random fingerprints of infinite length.

        Raises
        ------
        TypeError
            If fingerprints is not a numpy array.
            If the elements of fingerprints are not numpy arrays.
        ValueError
            If fingerprints is not a positive integer.
            If less than two fingerprints are provided.
            If not all the fingerprints have the same length.
        """
        if isinstance(fingerprints, int):
            if fingerprints <= 0:
                raise ValueError("If fingerprints is given as an integer,"
                                 "it should be positive integer")
            self.fingerprints = fingerprints
        else:
            if not isinstance(fingerprints, np.ndarray):
                raise TypeError("Fingerprints must be a numpy array or an int.")
            if not all(isinstance(fingerprint, np.ndarray) for fingerprint in fingerprints):
                raise TypeError("The elements of fingerprints must be a numpy array.")
            if len(fingerprints) < 2:
                raise ValueError("A minimum of 2 fingerprints must be provided.")
            if not all([len(fingerprint) == len(fingerprints[0]) for fingerprint in fingerprints]):
                raise ValueError("All the fingerprints must have the same length.")
            self.fingerprints = fingerprints

    def assign_c_threshold(self, c_threshold):
        """Assign coincidence threshold.

        Parameters
        ----------
        c_threshold : {None, 'dissimilar', int}
            Coincidence threshold.
            None : Default, c_threshold = n_fingerprints % 2
            'dissimilar' : c_threshold = ceil(n_fingerprints / 2)
            int : Integer number < n_fingerprints

        Raises
        ------
        TypeError
            If c_threshold is not None, 'dissimilar', or an integer.
        ValueError
            If c_threshold is an integer equal or greater than n_fingerprints
        """
        if not c_threshold:
            self.c_threshold = (self.n_characters - self.n_fingerprints
                                % self.n_characters) % self.n_characters
        if isinstance(c_threshold, str):
            if c_threshold != 'dissimilar':
                raise TypeError("c_threshold must be None, 'dissimilar', or an integer.")
            else:
                self.c_threshold = ceil(self.n_fingerprints / self.n_characters)
        if isinstance(c_threshold, int):
            if c_threshold >= self.n_fingerprints:
                raise ValueError("c_threshold cannot be equal or greater than n_fingerprints.")
            self.c_threshold = c_threshold

    def set_general_coincidences(self):
        """Calculate the general coincidence matrix."""
        general_coincidences = []
        for char in self.characters:
            general_coincidences.append(np.count_nonzero(self.fingerprints == char, axis=0))
        self.general_coincidences = general_coincidences

    def set_reduced_coincidences(self):
        """Calculate the reduced coincidence vector."""
        self.reduced_coincidences = np.amax(self.general_coincidences, axis=0)

    def set_matches(self):
        """Calculate the matches between the fingerprints.

        match[k] indicates how many times there are n-k characters matching.
        """
        max_k = self.n_fingerprints
        if self.n_fingerprints % self.n_characters == 0:
            min_k = int(self.n_fingerprints / self.n_characters)
        else:
            min_k = self.n_fingerprints // self.n_characters + 1
        matches = [0] * (max_k - min_k + 1)
        for k in range(max_k, min_k - 1, -1):
            matches[self.n_fingerprints - k] = np.count_nonzero(self.reduced_coincidences == k)
        self.matches = np.array(matches)        

    def set_d_vector(self):
        """Calculate the d vector.

        Notes
        -----
        The entries of this vector are the numbers t*Kb - n,
        which measure the degree of coincidence between the given fingerprints.
        """
        self.d_vector = np.array([self.n_characters * (self.n_fingerprints - i) -
                                  self.n_fingerprints for i in range(len(self.matches))])

    def set_w_factor(self, w_factor):
        """Calculate weight factors.

        Parameters
        ----------
        w_factor : {"fraction", "power_n"}
            Type of weight function that will be used.
            'fraction' : similarity = d[k]/n
                         dissimilarity = 1 - (d[k] - n_fingerprints % 2)/n_fingerprints
            'power_n' : similarity = n**-(n_fingerprints - d[k])
                        dissimilarity = n**-(d[k] - n_fingerprints % 2)
            other values : similarity = dissimilarity = 1
        """
        if w_factor == "power_n":
            power = int(w_factor.split("_")[-1])

            def f_s(d):
                return power**-(self.n_fingerprints - d)

            def f_d(d):
                return power**-(d - self.n_fingerprints % 2)
        elif w_factor == "fraction":
            def f_s(d):
                return d/(self.n_fingerprints * (self.n_characters - 1))

            def f_d(d):
                return 1 - (d - (self.n_characters - self.n_fingerprints
                                 % self.n_characters) % self.n_characters) /\
                           (self.n_fingerprints * (self.n_characters - 1))
        else:
            def f_s(d):
                return 1

            def f_d(d):
                return 1
        weights = len(self.matches) * [0]
        for k in range(len(self.matches)):
            if self.d_vector[k] > self.c_threshold:
                weights[k] = f_s(self.d_vector[k])
            else:
                weights[k] = f_d(self.d_vector[k])
        self.weights = np.array(weights)

    def set_weighted_matches(self):
        """Calculate weighted matches."""
        self.weighted_matches = self.matches * self.weights

    def set_counters(self):
        """Calculate the counters."""
        total_sim = 0
        total_dis = 0
        total_w_sim = 0
        total_w_dis = 0
        for k in range(len(self.matches)):
            if self.d_vector[k] > self.c_threshold:
                total_sim += self.matches[k]
                total_w_sim += self.weighted_matches[k]
            else:
                total_dis += self.matches[k]
                total_w_dis += self.weighted_matches[k]
        self.total_sim = total_sim
        self.total_dis = total_dis
        self.total_w_sim = total_w_sim
        self.total_w_dis = total_w_dis

    def set_p(self):
        """Calculate p."""
        self.p = self.total_sim + self.total_dis

    def set_weighted_p(self):
        """Calculate weighted p."""
        self.w_p = self.total_w_sim + self.total_w_dis
