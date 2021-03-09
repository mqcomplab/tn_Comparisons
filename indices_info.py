"""Contains miscellaneous information about the implemented indices
[index-name] : (ClassName, Abbreviation, (ways of calculating the index,))
"""

Indices = {"Austin-Colwell": ("AustinColwell", "AC", ("1sim_wdis", "1sim_dis")),
           "Consoni-Todeschini1": ("ConsoniTodeschini1", "CT1", ("1sim_wdis", "1sim_dis")),
           "Consoni-Todeschini2": ("ConsoniTodeschini2", "CT2", ("1sim_wdis", "1sim_dis")),
           "Sokal-Michner": ("SokalMichner", "SM", ("1sim_wdis", "1sim_dis")),
           "Sokal-Sneath2": ("SokalSneath2", "SS2", ("1sim_wdis", "1sim_dis")),
           "Rogers-Tanimoto": ("RogersTanimoto", "RT", ("1sim_wdis", "1sim_dis"))}
