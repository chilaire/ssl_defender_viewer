from graph import *
from problem import *
from math import sqrt, floor , pi
from geometry import *
import numpy as np

class createSol :
    def __init__(self, problem):
        self.problem = problem
        self.striking_shots = None
        self.graph = graph()
        self.step = problem.robot_radius / problem.pos_step
        self.offset = int(floor( 2 * problem.robot_radius / problem.pos_step)) #nb d'indices max entre 2 positions de robots qui se touchent
        self.solution = []

    """
    Check if two positions overlap (real_position)
    """
    def overlaps(self, pos1, pos2):
        d = sqrt((pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2)
        return d <= 2*self.problem.robot_radius


    """
    Get the list of striking shots
    """
    def get_striking_shots(self):
        self.striking_shots = []
        for opp_id in range(self.problem.getNbOpponents()):
            opp_pos = self.problem.getOpponent(opp_id)
            for goal in self.problem.goals :
                theta = 0    #TODO:calc theta1 theta2 = arctan(yg-ya/xg_xa) blabla
                while theta < 2*pi :
                    kick_result = goal.kickResult(opp_pos, theta)
                    if not kick_result is None:
                        self.striking_shots.append(kick_result)
                    theta += self.problem.theta_step

    """
    takes a real position and gives the grid indexes
    """
    def real_to_grid (self, x, y):
        j = x - self.problem.field_limits[0][0]
        i = y - self.problem.field_limits[1][0]
        j += self.problem.robot_radius #a robot can be out as long as it is in contact with the line
        i += self.problem.robot_radius
        return ( int(floor(i/self.problem.pos_step)) , int(floor(j/self.problem.pos_step)) )

    """
    takes a grid indexes and gives the real position
    """
    def grid_to_real (self, i, j):
        x = j*self.problem.pos_step + self.problem.field_limits[0][0] - self.problem.robot_radius
        y = i*self.problem.pos_step + self.problem.field_limits[1][0] - self.problem.robot_radius
        return ( x , y )

    """
    create position_grid : -2=forbiden -1=not in graph id=index in graph
    """
    def create_position_grid(self): #TODO:peut mieux faire !!!!!!
        xmax = self.problem.field_limits[0][1] + self.problem.robot_radius
        ymax = self.problem.field_limits[1][1] + self.problem.robot_radius
        (N,M) = self.real_to_grid( xmax , ymax )
        self.position_grid = np.full((N+1,M+1),-1)
        for opp_id in range(self.problem.getNbOpponents()):
            [xo,yo] = self.problem.getOpponent(opp_id)
            (imin,jmin) = self.real_to_grid(xo-2*self.problem.robot_radius, yo-2*self.problem.robot_radius)
            (imax,jmax) = (imin+2*self.offset+1,jmin+2*self.offset+1)
            for i in range(imin,min(imax+1,N)):
                for j in range(jmin,min(jmax+1,M)):
                    if ( self.position_grid[i,j] == -1 ):
                        (xd,yd) = self.grid_to_real(i,j)
                        if self.overlaps([xd,yd], [xo,yo]):
                            self.position_grid[i,j] = -2
        for goal in self.problem.goals:
            [x1,x2] = goal.posts[0]
            [y1,y2] = goal.posts[1]
            (imin1,jmin1) = self.real_to_grid(x1-self.problem.robot_radius, y1-self.problem.robot_radius)
            (imax1,jmax1) = (imin1+self.offset+1,jmin1+self.offset+1)
            (imin2,jmin2) = self.real_to_grid(x2-self.problem.robot_radius, y2-self.problem.robot_radius)
            (imax2,jmax2) = (imin2+self.offset+1,jmin2+self.offset+1)
            for i in range(imin1,min(imax1+1,N)):
                for j in range(jmin1,min(jmax1+1,M)):
                    if ( self.position_grid[i,j] == -1 ):
                        (xd,yd) = self.grid_to_real(i,j)
                        if (sqrt((xd - x1)**2 + (yd - y1)**2) < self.problem.robot_radius):
                            self.position_grid[i,j] = -2
            for i in range(imin2,min(imax2+1,N)):
                for j in range(jmin2,min(jmax2+1,M)):
                    if ( self.position_grid[i,j] == -1 ):
                        (xd,yd) = self.grid_to_real(i,j)
                        if (sqrt((xd - x2)**2 + (yd - y2)**2) < self.problem.robot_radius):
                            self.position_grid[i,j] = -2


    """
    Gets all the positions that intercept a shot (= segment from start_shot to end_shot)
    Returns a list of [i,j] that are grid_position

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

    """
    Adds a position  pos (grid_position) as neighbour of a shot in the graph.
    Finds all the positions that could overlap pos, to treat the case of overlapping
    (if two positions overlap, they are neighbours).

    def add(self, pos, id_shot):
        (need_to_treat_overlapping, index) = self.graph.add_pos(pos[0], pos[1], id_shot)
        if need_to_treat_overlapping :
            neighbours = []
            for k in range (floor(-2*self.step), floor(2*self.step)):
                for l in range (floor(-2*self.step), floor(2*self.step)):
                    if (k!=0 or l!=0) and k*k+l*l<= 4*self.step*self.step:
                        neighbours.append((k+pos[0],l+pos[1]))
            self.graph.add_pos_adja(index, pos, neighbours)"""

    """
    Adds a position  pos (grid_position) as neighbour of a shot in the graph.
    Finds all the positions that could overlap pos, to treat the case of overlapping
    (if two positions overlap, they are neighbours).
    """
    """
        def add_intercept(self, pos, shot): #exclure avant les -2
            id_shot = self.graph.add_shot()

            if self.grid_position[i,j]== -1: #need to treat overlappings
                #ajouter le sommet
                #traiter vois
            else: #just add an edge between the pos and the shot
                #ajouter le sommet
    """







    """
    Creates the graph.
    """
    def create_graph(self):
        self.get_striking_shots()
        self.create_position_grid()
        for shot in striking_shots:
            id_shot = self.graph.add_shot()
            pos_to_add = self.intercept(opp_pos,shot)
            for pos in pos_to_add :
                self.add(pos,id_shot)


    """
    Gets the solution found
    Converts the solution in real_position
    """
    def get_solution(self):
        if self.solution == []:
            return None
        return [self.grid_to_real(i,j) for (i,j) in self.solution]


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
