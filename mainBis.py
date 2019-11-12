#!/usr/bin/python3
import json
import sys
import pygame

from createSol import *
from solution import *
from board import *

if (len(sys.argv) < 2) :
    sys.exit("Usage: " + sys.argv[0] + " <problem.json>")

problem_path = sys.argv[1]

with open(problem_path) as problem_file:
    problem = Problem(json.load(problem_file))

s = createSol(problem)
s.create_graph()
for k in range(1,len(s.problem.opponents[0])*2+1):
    print(k)
    if s.dom_ind_set_glouton(k):
        break

solutionPb = s.get_solution()
print("solution :",solutionPb)

if solutionPb != None :

    x = {
        "defenders" : solutionPb
    }
    solution = Solution(x)


    b = Board(problem, solution)
    b.run()

sys.exit()
