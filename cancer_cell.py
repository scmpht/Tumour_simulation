# -----------------------------------------------------------
# Functions to implement genetic algorithm
#
# 2022, Morgan Thomas, University of Leeds
# scmpht@leeds.ac.uk
# Python version 3.8.12
# -----------------------------------------------------------

import numpy as np

class Cell:

    def __init__(self):
        """
        """
        self.fitness = -20
        self.cancerous = False
        self.genes = np.zeros(50)
        self.resistant = False

        self.mutation_rate = 1
        self.mutation_burden = 0

        self.oncogene = self.genes[0:10]
        self.tumour_sup = self.genes[10:30]
        self.resistance = self.genes[28:32]
        self.mutagenes = self.genes[32:34]


    def updateCell(self):
        """
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
        self.mutagenes = self.genes[32:34]


    def calcFitness(self):
        """
        """
        tumour_suppressors = np.count_nonzero(self.tumour_sup==0)
        oncogenes = np.count_nonzero(self.oncogene==1)
        fitness = oncogenes - tumour_suppressors
        return fitness
