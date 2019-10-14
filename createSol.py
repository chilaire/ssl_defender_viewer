from graph import *
from problem import *
from math import sqrt, floor

class createSol :
    def __init__(self, problem):
        self.problem = problem
        self.graph = graph()
        self.step = problem.robot_radius / problem.pos_step #pas floor ici

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
        return d <= 2*self.problem.robot_radius

    def add(pos, id_shot):
        (nfound, index) = self.graph.add_pos(pos[0], pos[1], id_shot)
        neighbours = [] #TODO est-ce qu'on a déjà les entiers à ce point de l'algo ?
        for k in range (floor(2*self.step)):
            for l in range (floor(2*self.step)):
                if k*k+l*l<= 4*self.step*self.step:
                    neighbours.append((k+pos[0],l+pos[1]))
        self.graph.add_pos_adja(index, pos, neighbours)

    def create_graph():
        for opp_id in range(self.problem.getNbOpponents()):
            opp_pos = self.problem.getOpponent(opp_id)
            for theta in range(0,2*pi,self.problem.theta_step) :
                for goal in self.problem.goals :
                    kick_result = self.problem.goal.kickResult(opp_pos, theta)
                    if not kick_result is None:
                        id_shot = self.graph.add_shot()
                        pos_to_add = intersept(opp_pos,kick_result)
                        for pos in pos_to_add :
                            add(pos,id_shot)
