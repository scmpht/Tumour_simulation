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

import time
import random
import sys
import argparse
import matplotlib.pyplot as plt
from copy import copy
import numpy as np
from genetic_algorithm import mutate, crossover
from cancer_cell import Cell
from pygame.locals import *
import pygame
import os
os.environ['PYGAME_HIDE_SUPPORT_PROMPT'] = "hide"


parser = argparse.ArgumentParser(description='Runs Cancer simulation.')

# add arguments
parser.add_argument('--dimensions', dest='dims', required=False)
parser.add_argument('--mutation', dest='mut', required=False)
args = parser.parse_args()

# set grid size
dimensions = (400, 400)
if args.dims and int(args.dims) > 8:
    dimensions = (int(args.dims), int(args.dims))

mutation_rate = 0.001
if args.mut and float(args.mut) < 1:
    mutation_rate = float(args.mut)


session = True
drug = False
spacing = 0                         # Sets space between each cell
cell_size = 4                       # Sets the drawn size of each cell
frame_rate = 30
start_time = time.time()

cells_x = int(dimensions[0]/cell_size)
cells_y = int(dimensions[1]/cell_size)
totalCells = cells_x*cells_y
cells = [[Cell() for x in range(cells_y)] for y in range(cells_x)]

statsRec = [[] for i in range(7)]
totalCells = cells_x*cells_y


pygame.init()
display = pygame.display.set_mode((cells_x*(cell_size+spacing),
                                   cells_y*(cell_size+spacing)))
display.fill(pygame.Color('black'))
clock = pygame.time.Clock()


def getNeighbours(x, y):
    """
    Returns a list of neighbouring cells including the cell itself
    """
    neighbours = []
    positions = []
    for i in [-1, 0, 1]:
        for j in [-1, 0, 1]:
            if x + i >= 0 and y + j >= 0:
                try:
                    neighbours.append(cells[x + i][y + j])
                    positions.append([x + i, y + j])
                except:
                    pass
    return neighbours, positions


def spread(x, y):
    """
    """
    neighbours, positions = getNeighbours(x, y)
    neighbFit = np.array([cell.calcFitness() for cell in neighbours])
    minIdx = np.where(neighbFit == neighbFit.min())
    minIdx = np.random.choice(minIdx[0], 1)
    weakest = positions[minIdx[0]]
    cells[weakest[0]][weakest[1]].genes = copy(cells[x][y].genes)


def cross(x, y):
    """
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
    """
    neighbours, _ = getNeighbours(x, y)
    neighbCell = np.array([not cell.cancerous for cell in neighbours])
    if random.random() > 0.5:
        return sum(neighbCell) > 2
    else:
        return sum(neighbCell) > 4


def die(x, y):
    """
    """
    cells[x][y] = Cell()


def cycle():
    """
    """
    cellCount = 0
    cancerCount = 0
    resistantCount = 0
    sumFitness = 0
    cellFitness = 0
    cancerFitness = 0
    resistantFitness = 0

    for x in range(cells_x):
        for y in range(cells_y):

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

            # counter for each type of cell

            if drug:  # Have to update attributes
                if cells[x][y].cancerous and cells[x][y].resistant:
                    spread(x, y)
                    cross(x, y)
                elif cells[x][y].cancerous and outsideCell(x, y):
                    die(x, y)
            elif cells[x][y].cancerous:
                spread(x, y)
                cross(x, y)
            mutate(cells[x][y], rate=mutation_rate)

    statsRec[0].append(sumFitness/totalCells)
    statsRec[1].append(cellFitness/totalCells)
    statsRec[2].append(cancerFitness/totalCells)
    statsRec[3].append(resistantFitness/totalCells)
    statsRec[4].append(cellCount)
    statsRec[5].append(cancerCount)
    statsRec[6].append(resistantCount)


def draw():
    """
    """
    for x in range(cells_x):
        for y in range(cells_y):
            intensity = cells[x][y].mutation_burden*4

            if cells[x][y].cancerous and cells[x][y].resistant:
                pygame.draw.rect(display,
                                 (0, intensity, 0),
                                 (x * (cell_size + spacing),
                                     y * (cell_size + spacing),
                                     cell_size, cell_size), 0)

            elif cells[x][y].cancerous:
                pygame.draw.rect(display,
                                 (intensity, 0, 0),
                                 (x * (cell_size + spacing),
                                     y * (cell_size + spacing),
                                     cell_size, cell_size), 0)

            else:
                pygame.draw.rect(display,
                                 (0, 0, intensity),
                                 (x * (cell_size + spacing),
                                     y * (cell_size + spacing),
                                     cell_size, cell_size), 0)

    pygame.display.update()
    clock.tick(frame_rate)


drugOn = []
drugOff = []


def makePlots(data):
    """
    """

    fig1, ax = plt.subplots()

    data[2] = [np.nan if x <= 0 else x for x in data[2]]
    data[3] = [np.nan if x <= 0 else x for x in data[3]]

    plt.plot(data[0], 'yellow', label='All cells')
    plt.plot(data[1], 'dodgerblue', label='Healthy cells')
    plt.plot(data[2], 'crimson', label='Cancer cells')
    plt.plot(data[3], 'limegreen', label='Resistant cells')
    for i in drugOn:
        plt.axvline(x=i, ls='--', color='black', alpha=0.7)
    for i in drugOff:
        plt.axvline(x=i, ls='--', color='black', alpha=0.3)
    plt.xlabel('Generation')
    plt.ylabel('Average fitness')
    plt.title('Cell fitness per generation')
    plt.legend(loc='center left')
    fig1.savefig("cell_fitness.png")

    fig2, ax = plt.subplots()
    plt.plot(data[4], 'dodgerblue', label='Healthy cells')
    plt.plot(data[5], 'crimson', label='Cancer cells')
    plt.plot(data[6], 'limegreen', label='Resistant cells')
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
    """
    """
    generation = 0
    global drug
    while session:
        pygame.display.set_caption('Tumour simulation: Generation {}  Time: {}'.format(str(generation),
                                                                                       round(time.time() - start_time, 2)))
        display.fill('black')
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                makePlots(statsRec)
                pygame.quit()
                sys.exit()
            if event.type == pygame.KEYDOWN and event.key == pygame.K_SPACE:
                drugOff.append(
                    generation) if drug else drugOn.append(generation)
                drug = not drug
                print('Toggle drug')

        cycle()
        draw()
        generation += 1


if __name__ == '__main__':
    run_automaton()
