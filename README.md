# Python-functions-based-on-Cantera

Author: Yang 


Email: yanggao.me@gmail.com

Extract reaction information from Cantera functions

including: reaction rate parameters, third-body species and efficiencies, fall-off parameters

reac.RA, reac.RB, reac.RE : forward temperature ABE factors

RE is converted and divided by RU, so

kf = A*T**B*exp(-E/T)

RA is kmole based, and converted to mole based in the test_get_reaction_info.py based on reaction order

reac.ITHB: number of enhanced third-body species, type: int

reac.NKTB: indices of enhanced third-body species, type: list

reac.AIK: efficiencies of enhanced third-body species, type: list

reac.Fall: first 3-> low pressure limit ABE, 4th-7th: alpha, T3, T1, T2 (T2 is only for 7 parameters Troe)

reaction type flag:

reac.isReversible = False

reac.isThirdbody = False

reac.isFalloff = False

reac.isChemical = False

reac.isPLOG = False

reac.isSimple = False

reac.isLindemann = False

reac.isTroe = False

reac.isTroe6 = False

reac.isTroe7 = False
