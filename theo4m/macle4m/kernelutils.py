"""

Miscellaneous kernel utils for ML models.

"""

import numpy as np
import sys
from munkres import Munkres

##########################################################################
#
# Originally devised for the EXAFS-guided MC geometry 
# optimization. See more on 26/07/2020-(1).
#
##########################################################################

#
# builds an environment covariance matrix for a 
# given specie; see Eq. (8) of [10.1039/c6cp00415f]
# so that the Munkres algorithm can be used to minimimze
def __min_cov_matrix__(SOAPs_A, SOAPs_B):
    rC_ij = [[np.dot(SOAPs_A[i], SOAPs_B[j])
                for j in range(len(SOAPs_B))] for i in range(len(SOAPs_A))]
    return rC_ij

#
# builds an environment covariance matrix for a 
# given specie; see Eq. (8) of [10.1039/c6cp00415f]
# so that the Munkres algorithm can be used to maximize
def __max_cov_matrix__(SOAPs_A, SOAPs_B):
    rC_ij = [[1.0 - np.dot(SOAPs_A[i], SOAPs_B[j])
                for j in range(len(SOAPs_B))] for i in range(len(SOAPs_A))]
    return rC_ij

##########################################################################
#
# Originally devised for the EXAFS-guided MC geometry 
# optimization. See more e-mail from Prof. Gabor on 27/07/2020.
#
##########################################################################

#
# see more at: 
# https://doi.org/10.1021/ci700274r
# https://doi.org/10.1039/c6cp00415f
# http://software.clapper.org/munkres/
# wtd ==> what to do? find the minimum or the maximum?
def min_max(atoms_A, atoms_B, S_A, S_B, wtd = "MIN", sort = False):

    # building an environment covariance matrix for 
    # each type of atom; see Eq. (8) of [10.1039/c6cp00415f]
    datar = []
    if wtd == "MIN":
        C_ij = __min_cov_matrix__(S_A, S_B)
    else:
        C_ij = __max_cov_matrix__(S_A, S_B)
    # Now checking the (min) best match matrix
    # Or the profit matrix within the context of 
    # the higher profit example at:
    # http://software.clapper.org/munkres/
    m = Munkres()
    min_bm = m.compute(C_ij)
    if wtd == "MIN":
        for row, column in min_bm:
            datar.append(C_ij[row][column])
    else:
        for row, column in min_bm:
            datar.append(1.0 - C_ij[row][column])
    if sort: datar.sort()
    return datar

##########################################################################
#
# Originally devised for the first tests on using the CUR 
# algorithm as implemented in ASAP (see more on 03/08/2020-(3)).
#
##########################################################################

#
# builds a structural covariance matrix using 
# the average structural kernel as given by 
# Eq. (9) of [10.1039/c6cp00415f] 
# 
# see more at: 
# https://doi.org/10.1021/ci700274r
# https://doi.org/10.1039/c6cp00415f
# http://software.clapper.org/munkres/
def average(structures, soapsl):
    desc = []
    for i in range(len(structures)):
        desc.append(soapsl.get_average_soap(0,i))
    covr = np.dot(np.asmatrix(desc), np.asmatrix(desc).T)
    return covr

#
# builds a structural covariance matrix using 
# the best-match structural kernel as given by 
# Eq. (10) of [10.1039/c6cp00415f] 
# 
# see more at: 
# https://doi.org/10.1021/ci700274r
# https://doi.org/10.1039/c6cp00415f
# http://software.clapper.org/munkres/
def best_match(structures, soapsl):
    soaps_strc = []
    for ist in range(len(structures)):
        # TODO: C necessary when Z != species_Z !!!
        S, C = soapsl.get_per_atom_soaps(0, ist)
        soaps_strc.append(S)
    covr = [[0.0 for j in range(len(structures))] for i in range(len(structures))]
    for i in range(len(structures)):
        for j in range(i, len(structures)): # it's a symmetric matrix
            datai = min_max(structures[i], structures[j], 
                    soaps_strc[i], soaps_strc[j], wtd = "MAX")
            covr[i][j] = sum(datai) / len(datai)
    for j in range(len(structures)):
        for i in range(j+1,len(structures)):
            covr[i][j] = covr[j][i]
    return covr

##########################################################################
#
# Originally devised for the first tests on using a boxplot-based kernel.
# See more on 28/07/2020-(10).
#
# IMPORTANT note: to compare two different structurues the kernel 
# must be normalized according to Eq. (3) of [DOI: 10.1039/c6cp00415f].
#
##########################################################################

# 
def __calc_boxplot_score__(structures, soaps_strc, i, j):
    datai = min_max(structures[i], structures[j], 
            soaps_strc[i], soaps_strc[j], wtd = "MIN", sort = True)
    #return 1.0 - datai[len(datai)-1]
    # removing duplicates
    datai = list(dict.fromkeys(datai))
    # interquartile range using numpy.percentile; see more at
    # https://www.geeksforgeeks.org/interquartile-range-and-quartile-deviation-using-numpy-and-scipy/
    # the interpolation method used make the data more compatible
    # with the box plots generated with the tool
    # http://shiny.chemgrid.org/boxplotr/
    # first quartile (Q1)
    Q1 = np.percentile(datai, 25.0, interpolation = "lower")
    # the median
    Q2 = np.percentile(datai, 50.0, interpolation = "lower")
    # third quartile (Q3)
    Q3 = np.percentile(datai, 75.0, interpolation = "lower")
    # interquaritle range (IQR)
    IQR = Q3 - Q1
    # Now using the IQR to detect outliers
    # https://www.geeksforgeeks.org/interquartile-range-to-detect-outliers-in-data/?ref=rp
    # find the lower and upper limits as Q1 – 1.5 IQR and Q3 + 1.5 IQR, respectively
    # https://www.nature.com/articles/nmeth.2807
    UL = Q3 + 1.5 * IQR
    LL = Q1 - 1.5 * IQR
    # and then the outliers
    UOL = [] # upper
    for iul in range(len(datai)-1, -1, -1):
        if datai[iul] > UL: # it's an outlier
            UOL.append(datai[iul])
        else:
            break
    LOL = [] # lower
    for ill in range(0, len(datai)):
        if datai[ill] < LL: # it's an outlier
            LOL.append(datai[ill])
        else:
            break
    # setting the score
    scr = 1.0 - datai[iul]
    scr += 1.0 - Q3
    scr += 1.0 - Q2
    scr += 1.0 - Q1
    scr += 1.0 - datai[ill]
    for io in range(len(LOL)): scr += 1.0 - LOL[io]
    return scr

#
# builds a structural covariance matrix using 
# the boxplot-based structural kernel
#
# see more at: 
# http://software.clapper.org/munkres/
def boxplot(structures, soapsl):
    soaps_strc = []
    for ist in range(len(structures)):
        # TODO: C necessary when Z != species_Z !!!
        S, C = soapsl.get_per_atom_soaps(0, ist)
        soaps_strc.append(S)
    covr = [[0.0 for j in range(len(structures))] for i in range(len(structures))]
    # diagonal
    for i in range(len(structures)):
        covr[i][i] = __calc_boxplot_score__(structures, soaps_strc, i, i)
    # off-diagonal (up)
    for i in range(len(structures)):
        for j in range(i + 1, len(structures)): # it's a symmetric matrix
            covr[i][j] = __calc_boxplot_score__(structures, soaps_strc, i, j)
            # normalizing
            #covr[i][j] = covr[i][j] / np.sqrt(covr[i][i]*covr[j][j])
        # normalizing
        #covr[i][i] = covr[i][i] / np.sqrt(covr[i][i]*covr[i][i])
    # off-diagonal (down)
    for j in range(len(structures)):
        for i in range(j+1,len(structures)):
            covr[i][j] = covr[j][i]
    return covr

