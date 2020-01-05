#!/usr/bin/python3
import json
import sys
import pygame
import time

from createSol import *
from solution import *
from board import *

if (len(sys.argv) < 2) :
    sys.exit("Usage: " + sys.argv[0] + " <problem.json>" + " [-greedy] (to use greedy aglo, exact algo by default)")

problem_path = sys.argv[1]

with open(problem_path) as problem_file:
    problem = Problem(json.load(problem_file))

tm = time.time()
s = createSol(problem)
possible = s.create_graph()
tm = time.time() - tm
print("creation graph : ", tm, "s")

tm = time.time()
if possible :
    if len(sys.argv) == 3 and sys.argv[2].lower() == "greedy" :
        print("Executing greedy method...")
        for k in range(len(s.problem.opponents[0])*100+1):
            print(k, " defenders")
            if s.dom_ind_set_glouton(k):
                break
    else :
        print("Executing exact method...")
        k = 0
        while not s.dom_ind_set(k):
            print(k , " defenders")
            k+=1

tm = time.time() - tm
print("solving : ", tm, "s")


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
