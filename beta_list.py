__author__ = 'lomovskayami'

import random
import bisect
import collections


def cdf(weights):
    total = sum(weights)
    result = []
    cumsum = 0
    for w in weights:
        cumsum += w
        result.append(cumsum / total)
    return result


def res_choice():
    weights = [0.9, 0.74, 1.02, 0.97, 0.75, 0.8, 1.08, 0.77, 1.49, 1.45, 1.32, 1.25, 1.14, 1.21, 0.92, 0.95, 0.72, 0.76, 0.64, 0.99]
    population = 'ACLMEQHKVIFYWTGSDNPR'
    assert len(population) == len(weights)
    cdf_vals = cdf(weights)
    x = random.random()
    idx = bisect.bisect(cdf_vals, x)
    return population[idx]

weights = [0.9, 0.74, 1.02, 0.97, 0.75, 0.8, 1.08, 0.77, 1.49, 1.45, 1.32, 1.25, 1.14, 1.21, 0.92, 0.95, 0.72, 0.76, 0.64, 0.99]
population = 'ACLMEQHKVIFYWTGSDNPR'



def random_res(length):
    """
    Makes seq using probabilities for amino acids in beta-sheet
    :param length: integer
    :return: list
    """
    res_list = [res_choice(population, weights) for i in range(length)]
    return res_list


#print random_res(124)