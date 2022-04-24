# -----------------------------------------------------------
# Functions to implement genetic algorithm.
#
# 2022, Morgan Thomas, University of Leeds
# scmpht@leeds.ac.uk
# Python version 3.8.12
# -----------------------------------------------------------

import numpy as np
import random


def crossover(cell1, cell2):
    """
    Takes in two Cell instances (see cancer_cell.py), splits their genes
    at the same random point, swaps the gene splits and
    rejoins them. Updates the two cells with their new crossed
    gene arrays.
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
    Takes a Cell instance (see cancer_cell.py),
    iterates through each gene, and with a probability
    of the 'rate' argument flips the value to either
    a zero or a one. Updates the cell with it's new genes.

    rate = float less than 1 that controls the rate of mutation.

    mutation rate is altered by the cells mutation_rate
    attribute determined by it's genes (see cancer_cell.py).
    This value is sqaured so mutation events are
    exponential not summative.
    """
    for i in range(len(cell.genes)):
        if random.random() < rate*(cell.mutation_rate**2):
            cell.genes[i] = 1 - cell.genes[i]
        else:
            pass

    cell.updateCell()
