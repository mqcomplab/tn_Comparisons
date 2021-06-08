import numpy as np
from itertools import combinations
from scipy.special import binom
#from indices.indices_info import Indices
from indices_info import Indices
from sokal_michner import SokalMichner
import random
import glob
import time


def _file_name_gen(key):
    """Generate the name of the file containing a given index.

    Arguments
    ---------
    key : str
        Key of the Indices dict that corresponds to a given index.

    Returns
    -------
    name : str
        Name of the file containing the class corresponding to the given index.
    """
    parts = key.split("-")
    if len(parts) == 1:
        name = key.lower()
    else:
        name = ""
        for part in parts[:-1]:
            name += part.lower()+"_"
        name += parts[-1].lower()
    return name


def generate_bitstring(size):
    """Generate a random fingerprint.

    Arguments
    ---------
    size : int
        Size (or length) of the fingerprint.
    
    Returns
    -------
    fingerprint : np.ndarray
        Fingerprint as a numpy array.
    """
    if not isinstance(size, int):
        raise TypeError("size can only be an integer.")
    return np.random.choice([0, 1], size=size)


def gen_fingerprints(fp_total, fp_size):
    """Generates random fingerprints.

    Arguments
    ---------
    fp_total : int
        Number of fingerprints.
    fp_size : int
        Size (or length) of the fingerprints.

    Returns
    -------
    total_fingerprints : np.array
        Numpy array containing the fingerprints.
    """
    total_fingerprints = []
    for i in range(fp_total):
        total_fingerprints.append(generate_bitstring(fp_size))
    total_fingerprints = np.array(total_fingerprints)
    return total_fingerprints


def calc_indices(indices=Indices, fp_total=2,
                 total_fingerprints=np.array([np.array([1]), np.array([1])]),
                 n=2, c_threshold=None, w_factor="fraction"):
    """Calculate the indices and generates the output.

    Arguments
    ---------
    indices : dict
        Dictionary with the indices that will be calculated.
    fp_total : int
        Total number of fingerprints.
    total_fingerprints : np.ndarray
        Numpy array containing the fingerprints that will be compared.
    n : int
        Number of fingerprints that will be compared simultaneously.
    c_threshold : {None, 'dissimilar', int}
        Coincidence threshold.
    w_factor : {"fraction", "power_n"}
        Type of weight function that will be used.

    Raises
    ------
    TypeError
        If n is not an integer.
    ValueError
        If the number of fingerprints in total_fingerprints is not equal to fp_total.
        If n is less than 2.
        If n is greater than f_total.

    Returns
    -------
    Results : dict
        Dictionary with the results of the comparisons.
    """
    if not isinstance(n, int):
        raise TypeError("n must be an integer.")
    if len(total_fingerprints) != fp_total:
        raise ValueError("The number of fingerprints in total_fingerprints must be equal"
                         "to fp_total.")
    if n < 2:
        raise ValueError("n cannot be less than 2.")
    if n > fp_total:
        raise ValueError("n cannot be greater than fp_total.")

    # Dictionary that will contain the results of all the comparisons.
    # Its structure is: Results[index] = (class_name, [index_values])
    Results = {}
    
    for s_index in indices:
        for variant in Indices[s_index][2]:
            Results[indices[s_index][1] + "_" + variant] = (Indices[s_index][0], [])

    # Sets of n numbers that indicate which fingerprints will be compared at a given time.
    index_list = list(combinations(range(fp_total), n))
    
    # Populating the Results dict with the results of the comparisons.
    for inds in index_list:
        fingerprints = total_fingerprints[list(inds)]
        for s_index in sorted(Results):
            h = "fingerprints=np.array([np.array([1, 0]), np.array([0, 1])]), "
            h += "c_threshold = 1, "
            h += "w_factor=None"
            h = "(" + h + ")"
            exec("index = " + Results[s_index][0] + h, None, globals())
            index.__init__(fingerprints=fingerprints, c_threshold=c_threshold, w_factor=w_factor)
            exec("result = index." + s_index, None, globals())
            Results[s_index][1].append(result)
    return Results


def _indices_values(results, fp_total=2, n=2, methods=[]):
    """Generate output with the results of the comparisons.

    Arguments
    ---------
    results : dict
        Dictionary with the results of the comparisons.
    fp_total : int
        Total number of fingerprints.
    n : int
        Number of fingerprints that will be compared simultaneously.
    methods : list
        List that contains the methods used to calculate the indices.

    Returns
    -------
    s : str
        String with the results of the comparisons of the given indices.
    """
    s = ""
    for method in methods:
        s += method + "      "
    s += "\n"
    for i in range(int(binom(fp_total, n))):
        s += "{:<13d}".format(i + 1)
        for method in methods:
            l = len(method)
            s += "{:^{}.6f}     ".format(results[method][1][i], l + 1)
        s += "\n"
    s += "\n             "
    for method in methods:
        s += method + "      "
    return s


def _max(results, methods=[]):
    """Generate the maxima of the results of the comparisons.

    Arguments
    ---------
    results : dict
        Dictionary with the results of the comparisons.
    methods : list
        List that contains the methods used to calculate the indices.

    Returns
    -------
    s : str
        String with the maxima of the comparisons of the given indices.
    """
    s = ""
    s += "\nMax          "
    for method in methods:
        max = np.amax(np.array(results[method][1]))
        l = len(method)
        s += "{:^{}.6f}     ".format(max, l + 1)
    return s


def _min(results, methods=[]):
    """Generate the minima of the results of the comparisons.

    Arguments
    ---------
    results : dict
        Dictionary with the results of the comparisons.
    methods : list
        List that contains the methods used to calculate the indices.

    Returns
    -------
    s : str
        String with the minima of the comparisons of the given indices.
    """
    s = ""
    s += "\nMin          "
    for method in methods:
        min = np.amin(np.array(results[method][1]))
        l = len(method)
        s += "{:^{}.6f}     ".format(min, l + 1)
    return s


def _abs_max(results, methods=[]):
    """Generate the abs maxima of the results of the comparisons.

    Arguments
    ---------
    results : dict
        Dictionary with the results of the comparisons.
    methods : list
        List that contains the methods used to calculate the indices.

    Returns
    -------
    s : str
        String with the abs maxima of the comparisons of the given indices.
    """
    s = ""
    s += "\nAbsMax       "
    for method in methods:
        absmax = np.amax(np.abs(np.array(results[method][1])))
        l = len(method)
        s += "{:^{}.6f}     ".format(absmax, l + 1)
    return s


def _abs_min(results, methods=[]):
    """Generate the abs minima of the results of the comparisons.

    Arguments
    ---------
    results : dict
        Dictionary with the results of the comparisons.
    methods : list
        List that contains the methods used to calculate the indices.

    Returns
    -------
    s : str
        String with the abs minima of the comparisons of the given indices.
    """
    s = ""
    s += "\nAbsMin       "
    for method in methods:
        absmin = np.amin(np.abs(np.array(results[method][1])))
        l = len(method)
        s += "{:^{}.6f}     ".format(absmin, l + 1)
    return s


def _average(results, methods=[]):
    """Generate the averages of the results of the comparisons.

    Arguments
    ---------
    results : dict
        Dictionary with the results of the comparisons.
    methods : list
        List that contains the methods used to calculate the indices.

    Returns
    -------
    s : str
        String with the averages of the comparisons of the given indices.
    """
    s = ""
    s += "\nAverage      "
    for method in methods:
        average = np.average(np.array(results[method][1]))
        l = len(method)
        s += "{:^{}.6f}     ".format(average, l + 1)
    return s


def _abs_average(results, methods=[]):
    """Generate the abs averages of the results of the comparisons.

    Arguments
    ---------
    results : dict
        Dictionary with the results of the comparisons.
    methods : list
        List that contains the methods used to calculate the indices.

    Returns
    -------
    String with the abs averages of the comparisons of the given indices.
    """
    s = ""
    s += "\nAbsAverage   "
    for method in methods:
        absaverage = np.average(np.abs(np.array(results[method][1])))
        l = len(method)
        s += "{:^{}.6f}     ".format(absaverage, l + 1)
    return s


def indices_output(results, fp_total=2, fp_size=1, n=2):
    """Generate output file with the results of the comparisons.

    Arguments
    ---------
    results : dict
        Dictionary with the results of the comparisons.
    fp_total : int
        Total number of fingerprints.
    fp_size : int
        Size of the fingerprints.
    n : int
        Number of fingerprints that will be compared simultaneously.
    """
    # Generic header for the output file.
    s = "Similarity analysis\n\n"
    s += "Fingerprint size (m)\n" + str(fp_size) + "\n\n"
    s += "Total number of fingerprints\n" + str(fp_total) + "\n\n"
    s += "Fingerprints compared simultaneously (n)\n" + str(n) + "\n\n"
    s += "#            " 
    
    # Weighted indices.
    w = []
    # Unweighted indices.
    no_w = []
    for s_index in sorted(results):
        if "w" in s_index:
            w.append(s_index)
        else:
            no_w.append(s_index)
    r_w = s
    r_no_w = s

    r_w += _indices_values(results=results, fp_total=fp_total, n=n, methods=w)
    r_no_w += _indices_values(results=results, fp_total=fp_total, n=n, methods=no_w)
    
    r_w += _max(results=results, methods=w)
    r_no_w += _max(results=results, methods=no_w)

    r_w += _abs_max(results=results, methods=w)
    r_no_w += _abs_max(results=results, methods=no_w)

    r_w += _min(results=results, methods=w)
    r_no_w += _min(results=results, methods=no_w)

    r_w += _abs_min(results=results, methods=w)
    r_no_w += _abs_min(results=results, methods=no_w)

    r_w += _average(results=results, methods=w)
    r_no_w += _average(results=results, methods=no_w)

    r_w += _abs_average(results=results, methods=w)
    r_no_w += _abs_average(results=results, methods=no_w)

    with open("wFP"+str(fp_total)+"m"+str(fp_size)+"n"+str(n)+".sim", "w") as outfile:
        outfile.write(r_w)

    with open("nwFP"+str(fp_total)+"m"+str(fp_size)+"n"+str(n)+".sim", "w") as outfile:
        outfile.write(r_no_w)


def read_fps(file):
    with open(file, "r") as infile:
        raw_lines = infile.readlines()
    #fp_total = len(raw_lines)
    total_fingerprints = []
    for line in raw_lines:
        s_fp = line.strip().split()[-1]
        np_fp = list(s_fp)
        total_fingerprints.append(np_fp)
    fp_total = len(total_fingerprints)
    fp_size = len(total_fingerprints[0])
    total_fingerprints = np.array(total_fingerprints)
    return fp_total, fp_size, total_fingerprints


def time_output(t_summary, comp, fp_type, t_out_name):
    if t_out_name:
        pass
    else:
        names = []
        for key in t_summary:
            if fp_type == "MACCS":
                name = key.split("_")[0] + "_" + key.split("_")[2]
            else:
                name = key.split("_")[0] + "_" + key.split("_")[2] + "_" + key.split("_")[3]
            if name not in names:
                names.append(name)
        if len(names) == 1:
            t_out_name = names[0] + "_{}.time".format(comp)
        else:
            raise TypeError("Please, provide an unambiguous name for the output file.")
    s = "{:20}".format("fp_total")
    s += "{:20}".format("n")
    s += "{:20}\n".format("time (s)")
    for key in reversed(sorted(t_summary)):
        s += "{:<20}{:<20}{:<20.6f}\n".format(t_summary[key]["fp_total"],
                                              t_summary[key]["n"], t_summary[key]["time"])
    with open(t_out_name, "w") as outfile:
        outfile.write(s)


def get_new_indices_n(total_fingerprints, selected_n, select_from_n):
    sim_values = []
    for i in select_from_n:
        new_indices = np.append(selected_n, i)
        new_fingerprints = total_fingerprints[new_indices]
        sim_index = SokalMichner(new_fingerprints, characters=characters)
        sim_values.append(sim_index.SM_1sim_dis)
    min_sim = min(sim_values)
    min_list = [j for j, v in enumerate(sim_values) if v == min_sim]
    return select_from_n[min_list]


def get_new_indices_b_max(total_fingerprints, selected_b, select_from_b):
    all_comps = []
    for i in select_from_b:
        comps = []
        for j in selected_b:
            new_indices = [i, j]
            new_fingerprints = total_fingerprints[new_indices]
            sim_index = SokalMichner(new_fingerprints, characters=characters)
            comps.append(sim_index.SM_1sim_dis)
        all_comps.append(comps)
    sim_values = [max(comps) for comps in all_comps]
    min_sim = min(sim_values)
    min_list = [j for j, v in enumerate(sim_values) if v == min_sim]
    return select_from_b[min_list]

def get_new_indices_b_sum(total_fingerprints, selected_b, select_from_b):
    all_comps = []
    for i in select_from_b:
        comps = []
        for j in selected_b:
            new_indices = [i, j]
            new_fingerprints = total_fingerprints[new_indices]
            sim_index = SokalMichner(new_fingerprints, characters=characters)
            comps.append(sim_index.SM_1sim_dis)
        all_comps.append(comps)
    sim_values = [sum(comps) for comps in all_comps]
    min_sim = min(sim_values)
    min_list = [j for j, v in enumerate(sim_values) if v == min_sim]
    return select_from_b[min_list]


def diver_output(name, select, n_averages, b_max_averages, b_sum_averages, n_std, b_max_std, b_sum_std):
    s = "Diversity Analysis\n\n"
    s += "Selecting {} molecules from the set {}\n\n".format(select, name.split(".")[0])
    s += "n-ary similarity  "
    for r in n_averages:
        s += "{:^{}.6f}".format(r, 11)
    s += "\n"
    s += "n-ary std         "
    for r in n_std:
        s += "{:^{}.6f}".format(r, 11)
    s += "\n\n"
    s += "binary MaxMin     "
    for r in b_max_averages:
        s += "{:^{}.6f}".format(r, 11)
    s += "\n"
    s += "MaxMin std        "
    for r in b_max_std:
        s += "{:^{}.6f}".format(r, 11)
    s += "\n\n"
    s += "binary MaxSum     "
    for r in b_sum_averages:
        s += "{:^{}.6f}".format(r, 11)
    s += "\n"
    s += "MaxSum std        "
    for r in b_sum_std:
        s += "{:^{}.6f}".format(r, 11)
    s += "\n"
    
    #ratio = []
    #for b, n in zip(b_results, n_results):
    #    if np.isnan(b/n):
    #        ratio.append("-")
    #    else:
    #        ratio.append(b/n)
    #s += "S(bin)/S(n)       "
    #for r in ratio:
    #    s += "{:^{}.6}".format(r, 11)
    #s += "\n"
    name = "DA_" + str(select) + "_" + name
    with open(name, "w") as outfile:
        outfile.write(s)    


if __name__ == "__main__":
    # Imports the classes corresponding to the similarity indices.
    for key in Indices:
        # s = "from indices." + file_name_gen(key) + " import " + Indices[key][0]
        s = "from " + _file_name_gen(key) + " import " + Indices[key][0]
        exec(s)

    # Sample run with randomly generated fingerprints.

    # Coincidence threshold.
    c_threshold = None

    # Weight factor.
    w_factor = "fraction"
    
    # T = 20
    characters = ['G', 'A', 'V', 'L', 'I', 'M', 'F', 'Y', 'W', 'H', 'D', 'E', 'K', 'R', 'N', 'Q', 'S', 'T', 'C', 'P']
    
    # T = 8
    #characters = ['H', 'A', 'N', 'P', 'G', 'D', 'B', 'R']
    
    # T = 4
    #characters = ['A', 'T', 'G', 'C']
    
    t = len(characters)

    # List of files with the fingerprints.
    fp_file_list = ["CYP2_n3296_t20.txt"]
    
    # Repetitions
    reps = 3
    
    # Selected molecules
    for select in [100]:
        # Similarity index
        sim_index = "Sokal-Michner"
        
        
        if fp_file_list:
            for fp_file in fp_file_list:
                total_n = []
                total_b_sum = []
                total_b_max = []
                start = time.time()
                for rep in range(reps):
                    fp_total, fp_size, total_fingerprints = read_fps(fp_file)
                    total_indices = np.array(range(fp_total))
                    seed = random.randint(0, fp_total - 1)
                    selected_n = np.array([seed])
                    selected_b_max = np.array([seed])
                    selected_b_sum = np.array([seed])
                    n_results = []
                    b_max_results = []
                    b_sum_results = []
                    while len(selected_n) < select:
                        #print(len(selected_n))
                        #print(total_indices)
                        #print(selected_n)
                        select_from_n = np.delete(total_indices, selected_n)
                        select_from_b_max = np.delete(total_indices, selected_b_max)
                        select_from_b_sum = np.delete(total_indices, selected_b_sum)
                        new_list_n = get_new_indices_n(total_fingerprints, selected_n, select_from_n)
                        new_list_b_max = get_new_indices_b_max(total_fingerprints, selected_b_max, select_from_b_max)
                        new_list_b_sum = get_new_indices_b_sum(total_fingerprints, selected_b_sum, select_from_b_sum)
                        new_index_n = random.choice(new_list_n)
                        new_index_b_max = random.choice(new_list_b_max)
                        new_index_b_sum = random.choice(new_list_b_sum)
                        selected_n = np.append(selected_n, new_index_n)
                        selected_b_max = np.append(selected_b_max, new_index_b_max)
                        selected_b_sum = np.append(selected_b_sum, new_index_b_sum)
                        new_fingerprints_n = total_fingerprints[selected_n]
                        sim_index_n = SokalMichner(new_fingerprints_n, characters=characters)
                        new_fingerprints_b_max = total_fingerprints[selected_b_max]
                        sim_index_b_max = SokalMichner(new_fingerprints_b_max, characters=characters)
                        new_fingerprints_b_sum = total_fingerprints[selected_b_sum]
                        sim_index_b_sum = SokalMichner(new_fingerprints_b_sum, characters=characters)
                        n_results.append(sim_index_n.SM_1sim_dis)
                        b_max_results.append(sim_index_b_max.SM_1sim_dis)
                        b_sum_results.append(sim_index_b_sum.SM_1sim_dis)
                        #print(sim_index_b.SM_1sim_dis/sim_index_n.SM_1sim_dis)
                    #print(time.time() - start)
                    new_fingerprints_n = total_fingerprints[selected_n]
                    sim_index_n = SokalMichner(new_fingerprints_n, characters=characters)
                    new_fingerprints_b_max = total_fingerprints[selected_b_max]
                    sim_index_b_max = SokalMichner(new_fingerprints_b_max, characters=characters)
                    new_fingerprints_b_sum = total_fingerprints[selected_b_sum]
                    sim_index_b_sum = SokalMichner(new_fingerprints_b_sum, characters=characters)
                    n_results.append(sim_index_n.SM_1sim_dis)
                    b_max_results.append(sim_index_b_max.SM_1sim_dis)
                    b_sum_results.append(sim_index_b_sum.SM_1sim_dis)
                    total_n.append(n_results)
                    total_b_max.append(b_max_results)
                    total_b_sum.append(b_sum_results)
                total_b_max = np.array(total_b_max)
                total_b_sum = np.array(total_b_sum)
                total_n = np.array(total_n)
                b_max_averages = total_b_max.mean(axis=0)
                b_sum_averages = total_b_sum.mean(axis=0)
                n_averages = total_n.mean(axis=0)
                b_max_std = total_b_max.std(axis=0)
                b_sum_std = total_b_sum.std(axis=0)
                n_std = total_n.std(axis=0)
                diver_output(fp_file, select, n_averages, b_max_averages, b_sum_averages, n_std, b_max_std, b_sum_std)
        else:
            pass