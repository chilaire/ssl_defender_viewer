from graph import *
from problem import *
from math import sqrt, floor , ceil, pi, cos, tan, sqrt
from geometry import *
import numpy as np

class createSol :
    def __init__(self, problem):
        self.problem = problem
        self.striking_shots = None
        self.limits_grid = np.full((2,2),0)
        self.position_grid = None
        self.graph = graph()
        self.step = problem.robot_radius / problem.pos_step
        self.offset = int(floor( 2 * problem.robot_radius / problem.pos_step)) #nb d'indices max entre 2 positions de robots qui se touchent
        self.solution = []


    """
    ###################################
    axiliar functions
    ###################################
    """


    def range_in_grid_y(self, y0, y1, radius):#TODO: un vrai nom
        (i0,i1) = ( ceil((y0-radius)/self.problem.pos_step), floor((y1+radius)/self.problem.pos_step) )
        [imin,imax] = self.limits_grid[0]
        return (max(i0,imin),min(i1,imax))

    def range_in_grid_x(self, x0, x1, radius):#TODO: un vrai nom
        (j0,j1) = ( ceil((x0-radius)/self.problem.pos_step), floor((x1+radius)/self.problem.pos_step) )
        [jmin,jmax] = self.limits_grid[1]
        return (max(j0,jmin),min(j1,jmax))

    """
    Check if two positions overlap (real position)
    """
    def overlaps(self, pos1, pos2):
        d = sqrt((pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2)
        return d <= 2*self.problem.robot_radius


    """
    Get the list of striking shots
    Each shot is represented by (begining_coordinates, end_coordinate, angle)
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
                        self.striking_shots.append((opp_pos, kick_result, theta))
                    theta += self.problem.theta_step



    """
    Create position_grid : -2=forbiden -1=not in graph id=index in graph
    WARNING! we use negatives indexes for negatives positions, so the grid is a translation by (ilimMim,jlimMin) of the field
    [[  0,0     0,1   ...    0,imax       0,imin   ...   0,-2     0,-1 ]
     [  0,0     0,1   ...    0,jmax       0,jmin   ...   0,-2   imax,-1]
     [  1,0     1,1   ...    1,jmax       1,jmin   ...   1,-2     1,-1 ]
     [  ...     ...   ...     ...           ...    ...   ...      ...  ]
     [imax,0  imax,1  ...  imax,jmax     imax,jmin ... imax,-2  imax,-1]

     [imin,0  imin,1  ...  imin,jmax     imin,jmin ... imin,-2  imin,-1]
     [  ...     ...   ...     ...           ...    ...   ...      ...  ]
     [ -2,0    -2,1   ...   -2,jmax       -2,jmin  ...  -2,-2    -2,-1 ]
     [ -1,0    -1,1   ...   -1,jmax       -1,jmin  ...  -1,-2    -1,-1 ]]
    """
    def create_position_grid(self): #TODO:peut mieux faire !!!!!!
        M = self.problem.field_limits[0][1] - self.problem.field_limits[0][0]#the center of the robot has to be in the field
        N = self.problem.field_limits[1][1] - self.problem.field_limits[1][0]
        (N,M) = ( int(floor(N/self.problem.pos_step)) , int(floor(N/self.problem.pos_step)) )
        self.position_grid = np.full((N+1,M+1),-1)
        #get the limits of the grid.
        ilimMin = ceil(self.problem.field_limits[1][0]/self.problem.pos_step)
        ilimMax = floor(self.problem.field_limits[1][1]/self.problem.pos_step)
        jlimMin = ceil(self.problem.field_limits[0][0]/self.problem.pos_step)
        jlimMax = floor(self.problem.field_limits[0][1]/self.problem.pos_step)
        self.limits_grid = [[ilimMin,ilimMax],[jlimMin,jlimMax]]

        #get rid of all the positions tha overlap an opposent
        for opp_id in range(self.problem.getNbOpponents()):
            [xo,yo] = self.problem.getOpponent(opp_id)
            (imin,imax) = self.range_in_grid_y(yo, yo, 2*self.problem.robot_radius)
            (jmin,jmax) = self.range_in_grid_x(xo, xo, 2*self.problem.robot_radius)
            for i in range(imin,imax):
                for j in range(jmin,jmax):
                    if ( self.position_grid[i,j] == -1 ):
                        (xd,yd)= (j * self.problem.pos_step , i * self.problem.pos_step)
                        if self.overlaps([xd,yd], [xo,yo]):
                            self.position_grid[i,j] = -2
        #idem with the pole of the goal
        for goal in self.problem.goals:
            for g in range(2):
                [xp,yp] = goal.posts[g]
                (imin,imax) = self.range_in_grid_y(yp, yp, self.problem.robot_radius)
                (jmin,jmax) = self.range_in_grid_x(xp, xp, self.problem.robot_radius)
                for i in range(imin,imax):
                    for j in range(jmin,jmax):
                        if ( self.position_grid[i,j] == -1 ):
                            (xd,yd)= (j * self.problem.pos_step , i * self.problem.pos_step)
                            if (sqrt((xd - xp)**2 + (yd - yp)**2) < self.problem.robot_radius):
                                self.position_grid[i,j] = -2



    """
    Finds alle the grid_position that intercept the shot "shot" and add them to the graph
    """
    def add_intercept(self, shot): #exclure avant les -2
        ([x0,y0],[x1,y1],theta) = shot
        id_shot = self.graph.add_shot()
        if cos(theta) == 0: #TODO verifier
            if y0 > y1:
                (y0,y1) = (y1,y0)
            if x0 > x1:
                (x0,x1) = (x1,x0)
            (i0,i1) = self.range_in_grid_y(y0, y1,2*self.problem.robot_radius)
            (j0,j1) = self.range_in_grid_x(x0, x1,2*self.problem.robot_radius)
            for j in range(j0, j1+1):
                for i in range(i0, i1+1):
                    if self.position_grid[i,j] != -2:
                        self.add(i,j,id_shot)
        else:
            coef = tan(theta)
            if abs(coef) > 1:
                if x0 > x1:
                    ([x0,y0],[x1,y1]) = ([x1,y1],[x0,y0])
                (j0,j1) = self.range_in_grid_x(x0, x1,2*self.problem.robot_radius)
                for j in range(j0, j1+1):
                    y = coef * j + y0
                    (i0,i1) = self.range_in_grid_y(y, y,2*self.problem.robot_radius)
                    for i in range (i0, i1+1):
                        if self.position_grid[i,j] != -2:
                            self.add(i,j,id_shot)
            else:
                coef = 1/coef
                if y0 > y1:
                    ([x0,y0],[x1,y1]) = ([x1,y1],[x0,y0])
                (i0,i1) = self.range_in_grid_y(y0, y1,2*self.problem.robot_radius)
                for i in range(i0, i1+1):
                    x = coef * i + x0
                    (j0,j1) = self.range_in_grid_y(x, x,2*self.problem.robot_radius)
                    for j in range (j0, j1+1):
                        if self.position_grid[i,j] != -2:
                            self.add(i,j,id_shot)

    """
    Adds a position  pos (position_grid) as neighbour of a shot in the graph.
    """
    def add(self, i, j, id_shot):
        index = self.position_grid[i,j]
        if index == -2:
            self.fail('we cannot add a forbidden position (grid contains -2))')
        elif index== -1:
            index = self.graph.add_pos(i,j,id_shot)
            self.position_grid[i,j] = index
            self.graph.add_edge_shot(index,id_shot)
        else:
            self.graph.add_edge_shot(self.position_grid[i,j],id_shot)

    """
    Add an edge between all the grid_position in the graph that overlaps
    For every grid_position p in the graph, we look for the neighbor in the semi-circle under the position
    """
    def treat_overlapping(self):
        (imax,jmin,jmax) = (self.limits_grid[0][1],self.limits_grid[1][0],self.limits_grid[1][1])
        for p in range (self.graph.get_nb_pos()):
            (ip,jp) = self.graph.pos_in_grid(p)
            offset_i =  min( self.offset, (imax - ip))
            for k in range( offset_i ):
                i = ip + k
                tmp = k*k*self.problem.pos_step*self.problem.pos_step
                tmp = sqrt(4*self.problem.robot_radius*self.problem.robot_radius - tmp)
                j1 = min (floor((jp*self.problem.pos_step + tmp)/self.problem.pos_step), jmax)
                j0 = max (ceil((jp*self.problem.pos_step - tmp)/self.problem.pos_step), jmin)
                for j in range(j0,j1):
                    index = self.position_grid[i,j]
                    if index > -1:
                        self.graph.add_edge_pos(index,p)


    """
    Creates the graph.
    """
    def create_graph(self):
        self.get_striking_shots()
        self.create_position_grid()
        for shot in self.striking_shots:
            id_shot = self.graph.add_shot()
            self.add_intercept(shot)
        self.treat_overlapping()


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
