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
import numpy as np

# class of reaction information
class ReactionInfo:
    def __init__(self):
        self.isReversible = False
        self.isThirdbody = False
        self.isFalloff = False
        self.isChemical = False
        self.isPLOG = False
        self.isSimple = False
        self.isLindemann = False
        self.isTroe = False
        self.isTroe6 = False
        self.isTroe7 = False

        # Forward reaction ABE factors
        self.RA = 0.0
        self.RB = 0.0
        self.RE = 0.0

        # Number of third body species
        self.ITHB = -1
        # Third body coefficients
        self.AIK = []
        # Third body species indices
        self.NKTB = []

        # Fall off reaction parameters
        # 1st:3rd : RKLOW ABE
        # 4th:7th : alpha, T3, T1, T2
        # Note that index starts from 0
        self.Fall = []


# function to get reaction information
def get_reaction_info(g, i):
# g is the initialized chemistry
# e.g.
# g = ct.Solution('h2.cti')
# i the reaction index - 1 (Python is 0 based)
# E is converted and divided by ruc in this function
# i.e.
# k = A*T**B*exp(-E/T)
# Note that A is kmol based!!!
    ruc = 8314.4621

    r1 = g.reaction(i)
    R_type = r1.reaction_type

    reac = ReactionInfo()

    if r1.reversible:
        reac.isReversible = True

    if R_type == 1:
        reac.isSimple = True

        rfhigh = r1.rate
        a_high = rfhigh.pre_exponential_factor
        b_high = rfhigh.temperature_exponent
        e_high = rfhigh.activation_energy / ruc

        reac.RA = a_high
        reac.RB = b_high
        reac.RE = e_high

    # simple third-body, not fall-off
    elif R_type == 2:
        reac.isThirdbody = True

        rfhigh = r1.rate
        a_high = rfhigh.pre_exponential_factor
        b_high = rfhigh.temperature_exponent
        e_high = rfhigh.activation_energy / ruc

        reac.RA = a_high
        reac.RB = b_high
        reac.RE = e_high

        efficiencies = r1.efficiencies
        third_body_species_keys = efficiencies.keys()
        third_body_efficiencies_values = efficiencies.values()
        third_body_species = list(third_body_species_keys)
        third_body_efficients = list(third_body_efficiencies_values)

        n_third = len(third_body_species)
        reac.ITHB = n_third

        for n in range(n_third):
            # index
            reac.NKTB.append( g.species_index(third_body_species[n]) )
            # coefficients
            reac.AIK.append ( third_body_efficients[n] )

    # fall-off
    elif R_type == 4:
        reac.isFalloff = True

        rfhigh = r1.high_rate
        a_high = rfhigh.pre_exponential_factor
        b_high = rfhigh.temperature_exponent
        e_high = rfhigh.activation_energy / ruc

        rflow = r1.low_rate
        a_low = rflow.pre_exponential_factor
        b_low = rflow.temperature_exponent
        e_low = rflow.activation_energy / ruc

        reac.RA = a_high
        reac.RB = b_high
        reac.RE = e_high

        reac.Fall.append(a_low)
        reac.Fall.append(b_low)
        reac.Fall.append(e_low)

        efficiencies = r1.efficiencies
        third_body_species_keys = efficiencies.keys()
        third_body_efficiencies_values = efficiencies.values()
        third_body_species = list(third_body_species_keys)
        third_body_efficients = list(third_body_efficiencies_values)

        n_third = len(third_body_species)
        reac.ITHB = n_third

        for n in range(n_third):
            # index
            reac.NKTB.append(g.species_index(third_body_species[n]))
            # coefficients
            reac.AIK.append(third_body_efficients[n])

        # Lindemann
        if r1.falloff.falloff_type == 100:
            reac.isLindemann = True
            # no further fall off parameters

        # Troe
        if r1.falloff.falloff_type == 110:
            reac.isTroe = True

            troe_parameters = r1.falloff.parameters
            alpha = troe_parameters[0]
            t3 = troe_parameters[1]
            t1 = troe_parameters[2]
            t2 = troe_parameters[3]

            if t2==0:
                reac.isTroe6 = True
            else:
                reac.isTroe7 = True

            reac.Fall.append(alpha)
            reac.Fall.append(t3)
            reac.Fall.append(t1)
            reac.Fall.append(t2)
    else:
        print ('Unknown reaction type, no supported!')
        return 0
    return reac
