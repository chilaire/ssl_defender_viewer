from graph import *
from problem import *
from math import sqrt, floor

class createSol :
    def __init__(self, problem):
        self.problem = problem
        self.graph = graph()
        self.step = floor(problem.robot_radius / problem.pos_step)

    def intercept(start_shot, end_shot): #hyp : start_shot and end_shot are integers
        #TODO : verify if x, y are in the field
        ret = []
        x_min = min(start_shot[0], end_shot[0]) - self.step
        x_max = max(start_shot[0], end_shot[0]) + self.step
        y_min = min(start_shot[1], end_shot[1]) - self.step
        y_max = max(start_shot[1], end_shot[1]) + self.step
        for x in range(x_min, x_max+1):
            for y in range(y_min, y_max+1):
                inter = segmentCircleIntersection(start_shot, end_shot, [x,y], self.problem.robot_radius)
                if inter != None:
                    add = False
                    for i in len(self.problem.opponents[0]):
                        d = sqrt((x-self.problem.opponents[0][i])**2 + (y-self.problem.opponents[1][i])**2)
                        add = add and (d > self.problem.robot_radius) #not optimal
                    if add:
                        ret.append([x,y])
        return ret
