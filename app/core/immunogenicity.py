import sys


def predict(state):
    immunoscale = {"A": 0.127, "C": -0.175, "D": 0.072, "E": 0.325, "F": 0.380, "G": 0.110, "H": 0.105, "I": 0.432, "K": -0.700, "L": -
                   0.036, "M": -0.570, "N": -0.021, "P": -0.036, "Q": -0.376, "R": 0.168, "S": -0.537, "T": 0.126, "V": 0.134, "W": 0.719, "Y": -0.012}
    immunoweight = [0.00, 0.00, 0.10, 0.31, 0.30, 0.29, 0.26, 0.18, 0.00]
    state.imm = {}

    for pep in state.all_peps:
        peptide = pep.upper()
        peplen = len(peptide)
        cterm = peplen - 1
        score = 0
        count = 0
        mask_num = [0, 1, cterm]
        mask_out = [1, 2, "cterm"]

        if peplen > 9:
            pepweight = immunoweight[:5] + \
                ((peplen - 9) * [0.30]) + immunoweight[5:]
        else:
            pepweight = immunoweight

        try:
            for pos in peptide:
                if pos not in immunoscale.keys():
                    print("!!!", pos)
                    raise KeyError()
                elif count not in mask_num:
                    score += pepweight[count] * immunoscale[pos]
                    count += 1
                else:
                    count += 1
            state.imm[peptide] = round(score, 5)
        except IOError as e:
            print("I/O error({0}): {1}".format(e.errno, e.strerror))
            raise ("Error: Please make sure you are entering in correct amino acids.")
        except:
            print("Unexpected error:", sys.exc_info()[0])
            raise

    return [k for k, v in sorted(state.imm.items(), key=lambda item: -item[1])]
    # result_list.sort(key=lambda tup: tup[-1], reverse=True)
    # state.imm = result_list
