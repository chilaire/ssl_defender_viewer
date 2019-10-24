#!/usr/bin/python3
import json
import sys
import pygame

from createSol import *

if (len(sys.argv) < 2) :
    sys.exit("Usage: " + sys.argv[0] + " <problem.json>")

problem_path = sys.argv[1]

with open(problem_path) as problem_file:
    problem = Problem(json.load(problem_file))

s = createSol(problem)
s.create_graph()
print(s.graph.adj_pos)
print("------------------------------------------")
print(s.graph.shot)
s.dom_ind_set(3)
print(s.solutionsList)

"""Convert to json
x = {
    "defenders" = [[x,y] for (x,y) in s.get_solution()]
}

solution = Solution(json.load(x))

b = Board(problem, solution)

b.run()
"""

sys.exit()
