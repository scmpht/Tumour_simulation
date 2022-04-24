# -----------------------------------------------------------
# Functions to implement genetic algorithm
#
# 2022, Morgan Thomas, University of Leeds
# scmpht@leeds.ac.uk
# Python version 3.8.12
# -----------------------------------------------------------

import numpy as np
import random


def crossover(cell1, cell2):
    """
    """
    cut = random.randint(0, len(cell1.genes))

    a1 = cell1.genes[0:cut]
    b1 = cell2.genes[0:cut]

    a2 = cell1.genes[cut:]
    b2 = cell1.genes[cut:]

    cell1.genes = np.append(a1, b2)
    cell2.genes = np.append(b1, a2)

    cell1.updateCell()
    cell2.updateCell()


def mutate(cell, rate=0.001):
    """
    """
    for i in range(len(cell.genes)):
        if random.random() < rate*cell.mutation_rate**2:
            cell.genes[i] = 1 - cell.genes[i]
        else:
            pass


    cell.updateCell()
