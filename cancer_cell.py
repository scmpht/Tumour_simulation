# -----------------------------------------------------------
# Class to control cell features in the simulation.
#
# 2022, Morgan Thomas, University of Leeds
# scmpht@leeds.ac.uk
# Python version 3.8.12
# -----------------------------------------------------------

import numpy as np

class Cell:
    """Class to control cell features in the simulation."""

    def __init__(self):
        """
        Initialised the class with default cell traits
        and genome.
        """
        self.fitness = -20 # Fitness function when all 0.
        self.cancerous = False
        self.resistant = False
        self.mutation_rate = 1
        self.mutation_burden = 0

        # Genome and gene sets.
        self.genes = np.zeros(50)
        self.oncogene = self.genes[0:10]
        self.tumour_sup = self.genes[10:30]
        self.resistance = self.genes[28:32]
        self.mutagenes = self.genes[48:50]

    def calcFitness(self):
        """
        Returns fitness of the cell.
        This is effectively the fitness function 
        of the genetic algorithm.
        Sums the mutated (ones) oncogenes and subtracts 
        the sum of the non-mutated (zeroes) of tumour
        suppressors.
        """
        tumour_suppressors = np.count_nonzero(self.tumour_sup==0)
        oncogenes = np.count_nonzero(self.oncogene==1)
        fitness = oncogenes - tumour_suppressors
        return fitness

    def updateCell(self):
        """
        Called to update the cell attributes for 
        the next generation. Called after mutate() 
        and crossover() functions (see genetic_algorithm.py).
        """
        self.fitness = self.calcFitness()
        self.cancerous = self.fitness > 0
        self.resistant = (np.count_nonzero(self.resistance==0) == 4)
        self.mutation_rate = np.count_nonzero(self.mutagenes) + 1
        self.mutation_burden = np.count_nonzero(self.genes)

        self.oncogene = self.genes[0:10]
        self.tumour_sup = self.genes[10:30]
        self.oncogene = self.genes[0:10]
        self.resistance = self.genes[28:32]
        self.mutagenes = self.genes[48:50]


   
