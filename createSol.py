from graph import *
from problem import *
from math import sqrt, floor

class createSol :
    def __init__(self, problem):
        self.problem = problem
        self.graph = graph()
        self.step = floor(problem.robot_radius / problem.pos_step)

    #hyp1 : start_shot and end_shot are integers
    #hyp2 : the shot starts and ends inside the field
    def intercept(start_shot, end_shot):
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
                        add = add and not superposition([x,y], [problem.opponents[*][i]]) #not optimal
                    if add:
                        ret.append([x,y])
        return ret

    def superposition(pos1, pos2):
        d = sqrt((pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2)
        return d <= self.problem.robot_radius

    
