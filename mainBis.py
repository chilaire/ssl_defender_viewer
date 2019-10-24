#!/usr/bin/python3
import json
import sys

from createSol import *

if (len(sys.argv) < 2) :
    sys.exit("Usage: " + sys.argv[0] + " <problem.json>")

problem_path = sys.argv[1]

with open(problem_path) as problem_file:
    problem = Problem(json.load(problem_file))

s = createSol(problem)
s.create_graph()
print(s.graph.adj_pos)
print(s.graph.shot)

sys.exit()
