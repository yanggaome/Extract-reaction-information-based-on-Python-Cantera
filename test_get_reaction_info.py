# Author: Yang Gao
# Email: yanggao.me@gmail.com
#################################################
# Extract reaction information from Cantera functions
# including: reaction rate parameters, third-body species and efficiencies, fall-off parameters
# reac.RA, reac.RB, reac.RE : forward temperature ABE factors
# reac.ITHB: number of enhanced third-body species, type: int
# reac.NKTB: indices of enhanced third-body species, type: list
# reac.AIK: efficiencies of enhanced third-body species, type: list
# reac.Fall: first 3-> low pressure limit ABE, 4th-7th: alpha, T3, T1, T2 (T2 is only for 7 parameters Troe)
# reaction type flag:
# reac.isReversible = False
# reac.isThirdbody = False
# reac.isFalloff = False
# reac.isChemical = False
# reac.isPLOG = False
# reac.isSimple = False
# reac.isLindemann = False
# reac.isTroe = False
# reac.isTroe6 = False
# reac.isTroe7 = False
#################################################
import cantera as ct
import numpy as np
import sys

# import function
from get_reaction_info import get_reaction_info

# initialize the chemistry, save in g
g = ct.Solution('h2-plog.cti')

# stoichiometric coefficients
# positive numbers
Sf = g.reactant_stoich_coeffs()

# number of species and reactions
KK = g.n_species
II = g.n_reactions

# open the file to write
file = open('get_reaction_info.txt','w')
    
# iterate reactions
# index starts from 0
# a[n] is the n+1th element in array a
# a[n,m] is n+1th to mth elements in array a
for i in range(0, II):
    r1 = g.reaction(i)
    # call the function to extract reaction information
    reac = get_reaction_info(g, i)

    coef_f = Sf[:,i]

    # np.nonzero returns a tuple, not an array
    # forward
    sp_index_f = np.nonzero(coef_f)
    sp_index_f = sp_index_f[0]
    r_coef_f = coef_f[sp_index_f]
    N_sp_f = len(r_coef_f)

    # compute reaction order
    # to convert the unit of pre-exponent factor A from kmol to to mole based
    reac_order = np.sum(r_coef_f)
    Convert_A = 1E3 ** (reac_order - 1)

    # get reaction type
    R_type = r1.reaction_type
    
    file.write ('C R%g: %s Reaction type: %g\n'% (i+1, r1.equation, R_type))

    if R_type == 1:
        A_high = reac.RA*Convert_A
        B_high = reac.RB
        E_high = reac.RE
        file.write('C forward temperature dependent reaction rate\n')
        file.write('      RF = %+.15E*T**%+.15E*EXP(%+.15E/T)\n' % (A_high, B_high, -E_high))

    # simple third-body, not fall-off
    elif R_type == 2:
        file.write('C have third-body\n')        
        A_high = reac.RA*1E3*Convert_A
        B_high = reac.RB
        E_high = reac.RE
        file.write('C forward temperature dependent reaction rate\n')
        file.write('      RF = %+.15E*T**%+.15E*EXP(%+.15E/T)\n' % (A_high, B_high, -E_high))
        file.write('C Third body information\n')
        N_third_body = reac.ITHB
        file.write('      Number of third-body enhanced species %g\n' % (N_third_body))
        for n in range(0, N_third_body):
            file.write('      Species index %g efficiency %g \n' % (reac.NKTB[n], reac.AIK[n]))

    # fall-off
    elif R_type == 4:
        file.write('C fall off\n')
        A_high = reac.RA * Convert_A
        B_high = reac.RB
        E_high = reac.RE

        A_low = reac.Fall[0]*1E3*Convert_A
        B_low = reac.Fall[1]
        E_low = reac.Fall[2]

        file.write('C forward temperature dependent reaction rate\n')
        file.write('      RF = %+.15E*T**%+.15E*EXP(%+.15E/T)\n' % (A_high, B_high, -E_high))
        file.write('C Third body information\n')
        N_third_body = reac.ITHB
        file.write('      Number of third-body enhanced species %g\n' % (N_third_body))
        for n in range(0, N_third_body):
            file.write('      Species index %g efficiency %g \n' % (reac.NKTB[n], reac.AIK[n]))

        if reac.isLindemann :
            file.write('C Lindemann 3-parameters\n      Low pressure limit ABE: %+.15E %+.15E %+.15E\n' % (A_low, B_low, E_low ) )
        if reac.isTroe6 :
            file.write('C Troe 6-parameters\n      Low pressure limit ABE: %+.15E %+.15E %+.15E\n' % (A_low, B_low, E_low))
            file.write('      alpha: %+.15E T3: %+.15E T1: %+.15E\n' % (reac.Fall[3], reac.Fall[4], reac.Fall[5] ) )

        if reac.isTroe7 :
            file.write('C Troe 7-parameters\n      Low pressure limit ABE: %+.15E %+.15E %+.15E\n' % (A_low, B_low, E_low))
            file.write('      alpha: %+.15E T3: %+.15E T1: %+.15E T2: %+.15E\n' % (reac.Fall[3], reac.Fall[4], reac.Fall[5], reac.Fall[6] ))
    elif R_type == 5:
        file.write('C PLOG\n')
        # convert list to array
        N_Plog = reac.NPLG
        P_plog = np.array(reac.P_plog)
        A_high = np.array(reac.A_plog)*Convert_A
        B_high = np.array(reac.B_plog)
        E_high = np.array(reac.E_plog)
        for n in range (0, N_Plog):
            file.write('      P = %g (Pa)\n' % (P_plog[n]))
            file.write('      RF = %+.15E*T**%+.15E*EXP(%+.15E/T)\n' % (A_high[n], B_high[n], -E_high[n]))
    else:
        print('Unknown reaction type, not supported\n')
        sys.exit()
    file.write('C\n')
file.close()

