from graph import *
from problem import *
from math import sqrt, floor , pi
from geometry import *

class createSol :
    def __init__(self, problem):
        self.problem = problem
        self.graph = graph()
        self.step = problem.robot_radius / problem.pos_step
        self.solution = []
        self.sol_inter = []
        self.current_state = 0

    """
    Check if two positions overlap (real_position)
    """
    def overlaps(self, pos1, pos2):
        d = sqrt((pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2)
        return d <= 2*self.problem.robot_radius


    """
    Gets all the positions that intercept a shot (= segment from start_shot to end_shot)
    Returns a list of [i,j] that are grid_position
    """
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
                        add = add and not self.overlaps([x,y], [self.problem.opponents[0][k], self.problem.opponents[1][k]]) #not optimal
                        for goal in self.problem.goals:
                            x1 = goal.posts[0][0]
                            y1 = goal.posts[1][0]
                            x2 = goal.posts[0][1]
                            y2 = goal.posts[1][1]
                            add = add and not (sqrt((x - x1)**2 + (y - y1)**2) < self.problem.robot_radius)
                            add = add and not (sqrt((x - x2)**2 + (y - y2)**2) < self.problem.robot_radius)
                    if add:
                        ret.append([i,j])
        return ret

    """
    Adds a position  pos (grid_position) as neighbour of a shot in the graph.
    Finds all the positions that could overlap pos, to treat the case of overlapping
    (if two positions overlap, they are neighbours).
    """
    def add(self, pos, id_shot):
        (need_to_treat_overlapping, index) = self.graph.add_pos(pos[0], pos[1], id_shot)
        if need_to_treat_overlapping :
            neighbours = []
            for k in range (floor(-2*self.step), floor(2*self.step)):
                for l in range (floor(-2*self.step), floor(2*self.step)):
                    if (k!=0 or l!=0) and k*k+l*l<= 4*self.step*self.step:
                        neighbours.append((k+pos[0],l+pos[1]))
            self.graph.add_pos_adja(index, pos, neighbours)

    """
    Creates the graph.
    """
    def create_graph(self):
        #we find all the shot drawn by opponents that might strike
        for opp_id in range(self.problem.getNbOpponents()):
            opp_pos = self.problem.getOpponent(opp_id)
            theta = 0
            while theta < 2*pi :
                for goal in self.problem.goals :
                    kick_result = goal.kickResult(opp_pos, theta)
                    if not kick_result is None:
                        # we add in G a shot-vertex, and a pos-vertex for each position that intercept the shot
                        id_shot = self.graph.add_shot()
                        pos_to_add = self.intercept(opp_pos,kick_result)
                        for pos in pos_to_add :
                            self.add(pos,id_shot)
                theta += self.problem.theta_step
        #now that G is finished and sorted, we convert the list of adjacencies
        self.graph.convert()
        self.sol_inter = [None for _ in range(2**sel.graph.get_nb_shot()]

    """
    Gets the solution found
    Converts the solution in real_position
    """
    def get_solution(self):
        if self.solution == []:
            return None
        return [[x*self.problem.pos_step,y*self.problem.pos_step] for (x,y) in self.solution]


    """
    Exact Algorithm :
    Finds an independent set of pos-vertices that dominates the shot-vertices
    Returns true and sets the solution with the index of the vertices in the solution if it finds it
         or false and sets the solution empty otherwise
    """
#here the adjacencies list are index
    def dom_ind_set(self, k):
        #if there is no more shot-vertex -> we have found the solution S -> True
        #else, we start from the neighbourhood of the first shot-vertex found
        shot_neighbours = self.graph.get_first_shot_neighbourhood()
        if shot_neighbours is None:
            return True
        #if k=0, no more defender to protect the shot found -> FALSE
        if k==0:
            return False
        #for each neighbour ui of v :
            #Gi = G\{ vi and its neighbours (pos and shot)}
        #we stop at the first vi where dom_ind_set(k-1) on Gi is true
        #   -> we can construct S with vi in S -> TRUE
        for n in shot_neighbours:
            (i,j) = self.graph.adj_pos[n][0]
            self.solution.append((i,j))
            self.graph.remove_vertex_and_neighbours(n,k)
            if self.dom_ind_set(k-1) :
                return True
            self.graph.revive_vertex_and_neighbours(n,k)
            self.solution.pop()
        #no neighbour or none can be in S
        #   -> impossible to find S -> FALSE
        return False


    """
    appoximate Algorithm
    """
    def dom_ind_set_glouton(self, k):
        #we get the best position with the heuristic
        #       and a bool that says if there is a shot to dominates
        (best_pos,found_shot) = self.graph.get_heuristic_pos()
        if not found_shot : #no shot to cover anymore => found S
            return True
        if k==0:
            return False
        if best_pos == None :#there is a shot but no position can protect it => fail
            return False
        (i,j) = self.graph.adj_pos[best_pos][0]
        self.solution.append((i,j))
        self.graph.remove_vertex_and_neighbours(best_pos,k)
        if self.dom_ind_set_glouton(k-1) :
            return True
        self.graph.revive_vertex_and_neighbours(best_pos,k)
        self.solution.pop()
        return False


    def dom_ind_set_dyn(self, k):
        #if there is no more shot-vertex -> we have found the solution S -> True
        #else, we start from the neighbourhood of the first shot-vertex found
        shot_neighbours = self.graph.get_first_shot_neighbourhood()
        if shot_neighbours is None:
            return True
        #if k=0, no more defender to protect the shot found -> FALSE
        if k==0:
            return False
        #for each neighbour ui of v :
            #Gi = G\{ vi and its neighbours (pos and shot)}
        #we stop at the first vi where dom_ind_set(k-1) on Gi is true
        #   -> we can construct S with vi in S -> TRUE
        for n in shot_neighbours:
            (i,j) = self.graph.adj_pos[n][0]
            self.solution.append((i,j))
            self.graph.remove_vertex_and_neighbours(n,k)
            self.current_state += 2**n
            tmp = self.sol_inter[self.current_state]
            if (tmp!=None and len(tmp)==k-1):
                self.sol_inter[self.current_state -2**n] = [(i,j)] + tmp
                return True
            elif self.dom_ind_set(k-1) :
                #TODO : update sol_inter
                return True
            sel.current_state -= 2**n
            self.graph.revive_vertex_and_neighbours(n,k)
            self.solution.pop()
        #no neighbour or none can be in S
        #   -> impossible to find S -> FALSE
        return False
