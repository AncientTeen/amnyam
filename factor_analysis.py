from __future__ import annotations

from correl_regress_Mult import *


def generality(buff: array[array[float]]) -> None:
    corrMatr = corrMatrix(buff)
    multiCorrelationCoff = multCorrCoff(corrMatr)

    h_sq = []
    for i in range(len(multiCorrelationCoff)):
        h_sq.append(multiCorrelationCoff[i] ** 2)




    pass
