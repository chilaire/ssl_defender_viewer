#!/usr/bin/python3
import json
import sys
import pygame
import time

from createSol import *
from solution import *
from board import *

if (len(sys.argv) < 2) :
    sys.exit("Usage: " + sys.argv[0] + " <problem.json>" + " (greedy \ exact)" + " default case:exact")

problem_path = sys.argv[1]

with open(problem_path) as problem_file:
    problem = Problem(json.load(problem_file))



tm = time.time()
s = createSol(problem)
s.get_striking_shots()
print(s.striking_shots)
s.create_position_grid()

s.create_graph()
print(s.position_grid)
for p in range(len(s.graph.adj_pos)):
    print(s.graph.adj_pos[p])
for p in range(len(s.graph.shot)):
    print(s.graph.shot[p])

tm = time.time() - tm
print("creation graph : ", tm, "s")
for k in range(1,len(s.problem.opponents[0])*2+1):
    print(k)
    if len(sys.argv) == 3 :
        if sys.argv[2].lower() == "greedy" :
            print("Executing greedy method...")
            if s.dom_ind_set_glouton(k):
                break
        else :
            print("Executing exact method...")
            if s.dom_ind_set(k):
                break
    else :
        print("Executing exact method...")
        if s.dom_ind_set(k):
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
