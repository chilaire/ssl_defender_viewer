from graph import *
from problem import *
from math import sqrt, floor , pi
from geometry import *

class createSol :
    def __init__(self, problem):
        self.problem = problem
        self.graph = graph()
        self.step = problem.robot_radius / problem.pos_step #pas floor ici

    #hyp1 : start_shot and end_shot are real_position (!= grid_position)
    #hyp2 : the shot starts and ends inside the field
    def intercept(self, start_shot, end_shot):
        ret = []
        i_min = int(floor( (min(start_shot[0], end_shot[0]) - self.problem.robot_radius) / self.problem.pos_step))
        i_max = int(floor( (max(start_shot[0], end_shot[0]) + self.problem.robot_radius) / self.problem.pos_step))
        j_min = int(floor( (min(start_shot[1], end_shot[1]) - self.problem.robot_radius) / self.problem.pos_step))
        j_max = int(floor( (max(start_shot[1], end_shot[1]) + self.problem.robot_radius) / self.problem.pos_step))
        for i in range(i_min, i_max+1):
            for j in range(j_min, j_max+1):
                (x,y) = (i*self.problem.pos_step,j*self.problem.pos_step)
                inter = segmentCircleIntersection(start_shot, end_shot, [x,y], self.problem.robot_radius)
                if inter is not None:
                    add = True
                    for k in range(len(self.problem.opponents[0])):
                        add = add and not self.superposition([x,y], [self.problem.opponents[0][k], self.problem.opponents[1][k]]) #not optimal
                    if add:
                        ret.append([i,j])
        return ret

    #real_positions
    def superposition(self, pos1, pos2):
        d = sqrt((pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2)
        return d <= 2*self.problem.robot_radius

    def add(self, pos, id_shot):
        (need_to_check_neighbours, index) = self.graph.add_pos(pos[0], pos[1], id_shot)
        if need_to_check_neighbours :
            neighbours = []
            for k in range (floor(2*self.step)):
                for l in range (floor(2*self.step)):
                    if k*k+l*l<= 4*self.step*self.step:
                        neighbours.append((k+pos[0],l+pos[1]))
            self.graph.add_pos_adja(index, pos, neighbours)


    def create_graph(self):
        for opp_id in range(self.problem.getNbOpponents()):
            opp_pos = self.problem.getOpponent(opp_id)
            theta = 0
            while theta < 2*pi :
                for goal in self.problem.goals :
                    kick_result = goal.kickResult(opp_pos, theta)
                    if not kick_result is None:
                        id_shot = self.graph.add_shot()
                        pos_to_add = self.intercept(opp_pos,kick_result)
                        for pos in pos_to_add :
                            self.add(pos,id_shot)
                theta += self.problem.theta_step
