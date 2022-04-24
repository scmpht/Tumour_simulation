# -----------------------------------------------------------
# Main script to be run for tumour growth simulation.
# Run this script by entering the following in the terminal 
# whilst in the directory where this file is stored:
# python cansim.py

# Additional flags:
#  --dimension (x, y) where x, y are dimensions of the output.
#  --mutation x where x < 1 and controls the rate of mutation.
#
# 2022, Morgan Thomas, University of Leeds
# scmpht@leeds.ac.uk
# Python version 3.8.12
# -----------------------------------------------------------

# Removing pygame welcome message.
import os
os.environ['PYGAME_HIDE_SUPPORT_PROMPT'] = "hide"

import pygame, sys, random, time
from pygame.locals import *
from cancer_cell import Cell
from genetic_algorithm import mutate, crossover
import numpy as np
from copy import copy
import matplotlib.pyplot as plt
import argparse

# Setting constants for simulation.
session = True
pause = False
drug = False
spacing = 0 
cell_size = 4
frame_rate = 30        
start_time = time.time()
dimensions = (400, 400) 
mutation_rate = 0.001


# Following section adds argument parsing in the terminal.
parser = argparse.ArgumentParser(description='Runs Cancer simulation.')

parser.add_argument('--dimensions', dest='dims', required=False)
parser.add_argument('--mutation', dest='mut', required=False)
args = parser.parse_args()

# Changing dimensions based on argument.
if args.dims and int(args.dims) > 8:
    dimensions = (int(args.dims), int(args.dims))

# Changing mutation rate based on argument.
if args.mut and float(args.mut) < 1:
    mutation_rate = float(args.mut)


# Instantiating grid of Cell objects where cells[x][y]
# is the cell at coordinates x, y.
cells_x = int(dimensions[0]/cell_size)
cells_y = int(dimensions[1]/cell_size)
totalCells = cells_x*cells_y
cells = [[Cell() for x in range(cells_y)] for y in range(cells_x)]

# The following are required to track stats for generated plots.
statsRec = [[] for i in range(7)]
totalCells = cells_x*cells_y

# Initialising pygame display.
pygame.init()
display = pygame.display.set_mode((cells_x*(cell_size+spacing), 
                                   cells_y*(cell_size+spacing)))
display.fill(pygame.Color('black'))
clock = pygame.time.Clock()


def getNeighbours(x, y):
    """
    Takes x, y coordinates of a cell and
    returns a list of neighbouring cells and 
    list of their positions, including the 
    cell itself.
    """ 
    neighbours = []
    positions = []
    for i in [-1,0,1]:
        for j in [-1, 0, 1]:
            # This if statement and try/except are required to 
            # prevent indexing cells not in the grid.
            if x + i >= 0 and y + j >= 0: 
                try:
                    neighbours.append(cells[x + i][y + j])
                    positions.append([x + i, y + j])
                except IndexError:
                    pass
    return neighbours, positions


def spread(x, y):
    """
    Takes x, y coordinates of a cell and
    converts the genes of the least fit 
    neighboring cell into a copy its own genes.
    """
    neighbours, positions = getNeighbours(x, y)
    neighbFit = np.array([cell.calcFitness() for cell in neighbours])
    minIdx = np.where(neighbFit == neighbFit.min())
    minIdx = np.random.choice(minIdx[0], 1)
    weakest = positions[minIdx[0]]
    cells[weakest[0]][weakest[1]].genes = copy(cells[x][y].genes) # Copy required to overcome mutability.


def cross(x, y):
    """
    Takes x, y coordinates of a cell and
    at performs the crossover function (see
    genetic_algorithm.py) with a random 
    neighbouring cancer cell.
    """
    neighbours, positions = getNeighbours(x, y)
    neighbCancer = np.array([cell.cancerous for cell in neighbours])
    idx = np.where(neighbCancer)
    idx = np.random.choice(idx[0], 1)
    partner = positions[idx[0]]
    if (x == partner[0] and y == partner[1]):
        return
    else:
        crossover(cells[x][y], cells[partner[0]][partner[1]])
    

def outsideCell(x, y):
    """
    Takes x, y coordinates of a cell and
    evaluates to true or false if it is on the
    tumour periphery.
    """
    neighbours, _ = getNeighbours(x, y)
    neighbCell = np.array([not cell.cancerous for cell in neighbours])
    # 50% chance that outside is defined by more than two adjacent healthy
    # cells, or more than four adjacent healthy cells.
    # This gave the tumour more life-like growth/shrink dynamics.
    if random.random() > 0.5:
        return sum(neighbCell)>2
    else:
        return sum(neighbCell)>4


def die(x, y):
    """Replaces a cell with a new instance of a cell."""
    cells[x][y] = Cell()


def cycle():
    """
    Main function of genetic algorithm,
    called once every generation.
    Iterates through each cell and determines behaviour
    based on rules (see README flow chart).
    """

    # Creating counters to track statistics.
    cellCount = 0
    cancerCount = 0
    resistantCount = 0
    sumFitness = 0
    cellFitness = 0
    cancerFitness = 0
    resistantFitness = 0 

    for x in range(cells_x):
        for y in range(cells_y):


            # Following section is require to track statistics.
            sumFitness = sumFitness + cells[x][y].fitness

            if cells[x][y].cancerous and cells[x][y].resistant:
                resistantFitness = resistantFitness + cells[x][y].fitness
                resistantCount = resistantCount + 1

            elif cells[x][y].cancerous:
                cancerFitness = cancerFitness + cells[x][y].fitness
                cancerCount = cancerCount + 1

            elif not cells[x][y].cancerous:
                cellFitness = cellFitness + cells[x][y].fitness
                cellCount = cellCount + 1
            
            # Following section is the rule-based behaviour,
            # see README for flow chart.
            if drug:
                if cells[x][y].cancerous and cells[x][y].resistant:
                    spread(x, y)
                    cross(x, y)
                elif cells[x][y].cancerous and outsideCell(x, y):
                    die(x, y)
            elif cells[x][y].cancerous:
                    spread(x, y)
                    cross(x, y)
            mutate(cells[x][y], rate=mutation_rate)

    # Appending statistics for this generation to the stats tracker.
    statsRec[0].append(sumFitness/totalCells)
    statsRec[1].append(cellFitness/totalCells)
    statsRec[2].append(cancerFitness/totalCells)
    statsRec[3].append(resistantFitness/totalCells)
    statsRec[4].append(cellCount)
    statsRec[5].append(cancerCount)
    statsRec[6].append(resistantCount)


def draw():
    """
    Called after each cycle to draw the cells 
    in the simulation by iterating through each cell
    drawing a rectangle for each cell.
    """

    for x in range(cells_x):
        for y in range(cells_y):
            # Higher mutational burden makes the cell brighter.
            intensity = cells[x][y].mutation_burden*4

            # Colour resitant cancerous cells green.
            if cells[x][y].cancerous and cells[x][y].resistant:
                pygame.draw.rect(display, 
                                (0, intensity, 0),
                                (x * (cell_size + spacing), 
                                 y * (cell_size + spacing), 
                                 cell_size, cell_size), 0)
            
            # Colour cancerous cells red.
            elif cells[x][y].cancerous:
                pygame.draw.rect(display, 
                                (intensity, 0, 0),
                                (x * (cell_size + spacing), 
                                 y * (cell_size + spacing), 
                                 cell_size, cell_size), 0)

            # Colour healthy cells blue.
            else:
                pygame.draw.rect(display, 
                                (0, 0, intensity),
                                (x * (cell_size + spacing), 
                                 y * (cell_size + spacing), 
                                 cell_size, cell_size), 0)
    
    pygame.display.update()
    clock.tick(frame_rate)

# Trackers to plot lines on graph when drug is toggled.
drugOn = []
drugOff = []

def makePlots(data):
    """
    Returns two plots from the stats data collected
    in the cycle function.

    cell_fitness.png tracks cell counts of each type 
    of cell over the generations.


    cell_counts.png tracks cell counts of each type 
    of cell over the generations.

    dashed lines indicate drug turning on and off.
    """

    fig1, ax = plt.subplots()

    # Necessary to remove lines when no cells of this type are present.
    data[2] = [np.nan if x <= 0 else x for x in data[2]]
    data[3] = [np.nan if x <= 0 else x for x in data[3]]


    # Cell fitness plot
    plt.plot(data[0], 'yellow', label = 'All cells')
    plt.plot(data[1], 'dodgerblue', label = 'Healthy cells')
    plt.plot(data[2], 'crimson', label = 'Cancer cells')
    plt.plot(data[3], 'limegreen', label = 'Resistant cells')
    # Plotting drug on and off lines.
    for i in drugOn:
        plt.axvline(x=i, ls='--', color='black', alpha=0.7)
    for i in drugOff:
        plt.axvline(x=i, ls='--', color='black', alpha=0.3)
    plt.xlabel('Generation')
    plt.ylabel('Average fitness')
    plt.title('Cell fitness per generation')
    plt.legend(loc='center left')
    fig1.savefig("cell_fitness.png")

    # Cell counts plot
    fig2, ax = plt.subplots()
    plt.plot(data[4], 'dodgerblue', label = 'Healthy cells')
    plt.plot(data[5], 'crimson', label = 'Cancer cells')
    plt.plot(data[6], 'limegreen', label = 'Resistant cells')
    # Plotting drug on and off lines.
    for i in drugOn:
        plt.axvline(x=i, ls='--', color='black', alpha=0.7)
    for i in drugOff:
        plt.axvline(x=i, ls='--', color='black', alpha=0.3)
    plt.xlabel('Generation')
    plt.ylabel('Cell counts')
    plt.title('Cell counts per generation')
    plt.legend(loc='center right')
    fig2.savefig('cell_counts.png')


def run_automaton():
    """Main function that runs the simulation."""
    generation = 0
    global drug
    global pause
    while session:
        pygame.display.set_caption('Tumour simulation: Generation {}  Time: {}'.format(str(generation), 
                                    round(time.time() - start_time,2)))
        display.fill('black')

        # Event tracking.
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                makePlots(statsRec)
                pygame.quit(); sys.exit();
            if event.type == pygame.KEYDOWN and event.key == pygame.K_SPACE:
                    drugOff.append(generation) if drug else drugOn.append(generation)
                    drug = not drug
                    print('Toggle therapy.')
            if event.type == pygame.KEYDOWN and event.key == pygame.K_RETURN:
                    pause = not pause
                    print('Toggle pause.')

        
        if not pause:
            # Update the generation.
            cycle()
            draw()
            generation += 1

if __name__ == '__main__':
    run_automaton()