from graph import *
from problem import *
from math import sqrt, floor , ceil, pi, cos, sin, tan, sqrt
from geometry import *
import numpy as np

class createSol :
    def __init__(self, problem):
        self.problem = problem
        self.striking_shots = None
        self.index_lim = np.full((2,2),0)
        self.grid = None
        self.graph = graph()
        self.solution = []
        self.step = self.problem.pos_step
        self.radius = self.problem.robot_radius

    """
    Gets the solution found.
    Converts the solution in real_position.
    """
    def get_solution(self):
        if self.solution == []:
            return None
        return [[j*self.step,i*self.step] for (i,j) in self.solution]

    """
    Given two ordinates (y0 <= y1) and a radius, returns the range of indexes [i0,i1]
    that corresponds to the positions on the grid in [y0-radius, y1+radius]
    """
    def range_in_grid_y(self, y0, y1, radius):
        (i0,i1) = ( int(ceil((y0-radius)/self.step)), int(floor((y1+radius)/self.step)) )
        [imin,imax] = self.index_lim[0]
        return (max(i0,imin),min(i1,imax))
    """
    Given two abscissas (x0 <= x1) and a radius, returns the range of indexes [j0,j1]
    that corresponds to the positions on the grid in [x0-radius, x1+radius]
    """
    def range_in_grid_x(self, x0, x1, radius):
        (j0,j1) = ( int(ceil((x0-radius)/self.step)), int(floor((x1+radius)/self.step)) )
        [jmin,jmax] = self.index_lim[1]
        return (max(j0,jmin),min(j1,jmax))

    """
    Gets the list of striking shots
    Each shot is represented by (begining_coordinates, end_coordinate, angle)
    """
    def get_striking_shots(self):
        self.striking_shots = []
        for opp_id in range(self.problem.getNbOpponents()):
            opp_pos = self.problem.getOpponent(opp_id)
            theta_range = [0,2*pi]
            for goal in self.problem.goals :
                """for k in range(2):
                    x=opp_pos[0]-goal.posts[0,k]
                    y=opp_pos[1]-goal.posts[1,k]
                    if x==0 :
                        theta_range[k] = pi/2
                        if y < 0 :
                            theta_range[k] = -pi/2
                    else:
                        theta_range[k] = np.arctan2(y,x)
                if theta_range[0] > theta_range[1]:
                    theta_range = [theta_range[1], theta_range[0]]
                if theta_range[0] - theta_range[0] > pi :
                    theta_range = [theta_range[1], 2*pi + theta_range[0]]
                theta = ceil(theta_range[0]/self.problem.theta_step)*self.problem.theta_step
                #print(theta_range)"""
                theta = theta_range[0]
                while theta <= theta_range[1] :
                    #print(theta)
                    kick_result = goal.kickResult(opp_pos, theta)
                    if not kick_result is None:
                        self.striking_shots.append((opp_pos, kick_result, theta+pi))
                    theta += self.problem.theta_step



    """
    Creates the grid : -2=forbiden   -1=not in graph   id=index in graph
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
    def create_grid(self):
        #get the limits of the grid. (a position is valid as long as the robot touches the field)
        ilimMin = int(ceil((self.problem.field_limits[1][0]-self.radius)/self.step))
        ilimMax = int(floor((self.problem.field_limits[1][1]+self.radius)/self.step))
        jlimMin = int(ceil((self.problem.field_limits[0][0]-self.radius)/self.step))
        jlimMax = int(floor((self.problem.field_limits[0][1]+self.radius)/self.step))
        self.index_lim = [[ilimMin,ilimMax],[jlimMin,jlimMax]]
        (N,M) = ( ilimMax - ilimMin , jlimMax - jlimMin )
        self.grid = np.full((N+1,M+1),-1)
        #get rid of all the positions that overlap an opposent
        for opp_id in range(self.problem.getNbOpponents()):
            [xo,yo] = self.problem.getOpponent(opp_id)
            (imin,imax) = self.range_in_grid_y(yo, yo, 2*self.radius)
            (jmin,jmax) = self.range_in_grid_x(xo, xo, 2*self.radius)
            for i in range(imin,imax+1):
                for j in range(jmin,jmax+1):
                    if ( self.grid[i,j] == -1 ):
                        (xd,yd)= (j * self.step , i * self.step)
                        d = sqrt((xd - xo)**2 + (yd - yo)**2)
                        if d < 2*self.radius:
                            self.grid[i,j] = -2
        #idem with the pole of the goal
        for goal in self.problem.goals:
            for g in range(2):
                (xp,yp) = (goal.posts[0][g],goal.posts[1][g])
                (imin,imax) = self.range_in_grid_y(yp, yp, self.radius)
                (jmin,jmax) = self.range_in_grid_x(xp, xp, self.radius)
                for i in range(imin,imax+1):
                    for j in range(jmin,jmax+1):
                        if ( self.grid[i,j] == -1 ):
                            (xd,yd)= (j * self.step , i * self.step)
                            if (sqrt((xd - xp)*(xd - xp) + (yd - yp)*(yd - yp)) < self.radius):
                                self.grid[i,j] = -2



    """
    Finds all the grid_position that intercept the shot "shot" and add them to the graph
    Returns True if there is at least 1 position that can intercept the shot, False otherwise
    """
    def add_intercept(self, shot):
        at_least_one_found = False
        ([x0,y0],[x1,y1],theta) = shot
        id_shot = self.graph.add_shot()
        if x0 == x1: #y=f(x) is a vertical line
            if y0 > y1:
                (y0,y1) = (y1,y0)
            (i0,i1) = self.range_in_grid_y(y0, y1, self.radius)
            (j0,j1) = self.range_in_grid_x(x0, x0, self.radius)
            for j in range(j0, j1+1):
                for i in range(i0, i1+1):
                    if self.grid[i,j] != -2:
                        self.add(i,j,id_shot)
                        at_least_one_found = True
        else:
            coef = tan(theta)
            if abs(coef) <= 1: # y=f(x) between horizontal and oblic -> for each x (j), many i corresponding to y=f(x)
                if x0 > x1:
                    (x0,x1,y0,y1) = (x1,x0,y1,y0)
                (j0,j1) = self.range_in_grid_x(x0, x1,0)
                (a,b) = ( coef * self.step , y0 - coef * x0)
                for j in range(j0, j1+1):
                    y = a * j + b
                    (i0,i1) = self.range_in_grid_y(y, y,self.radius)
                    for i in range (i0, i1+1):
                        if self.grid[i,j] != -2:
                            self.add(i,j,id_shot)
                            at_least_one_found = True
            else: # y=f(x) between vertical and oblic -> for each y (i), many y corresponding to x=f-1(y)
                coef = 1/coef
                if y0 > y1:
                    (x0,x1,y0,y1) = (x1,x0,y1,y0)
                (i0,i1) = self.range_in_grid_y(y0, y1,0)
                (a,b) = ( coef * self.step , x0 - coef * y0)
                for i in range(i0, i1+1):
                    x = a * i + b
                    (j0,j1) = self.range_in_grid_x(x, x,self.radius)
                    for j in range (j0, j1+1):
                        if self.grid[i,j] != -2:
                            self.add(i,j,id_shot)
                            at_least_one_found = True
        return at_least_one_found


    """
    Adds a position (rpz by its oordinates (i,j)) as neighbour of a shot (rpz by its index id_shot) in the graph.
    """
    def add(self, i, j, id_shot):
        index = self.grid[i,j]
        if index == -2:
            self.fail('we cannot add a forbidden position (grid contains -2))')
        elif index == -1:
            index = self.graph.add_pos(i,j,id_shot)
            self.grid[i,j] = index
            self.graph.add_edge_shot(index,id_shot)
        else:
            self.graph.add_edge_shot(self.grid[i,j],id_shot)

    """
    Adds an edge between all the grid_position in the graph that overlaps
    For every grid_position p in the graph, we look for the neighbor in the semi-circle under the position
    """
    def treat_overlapping(self):
        (imax,jmin,jmax) = (self.index_lim[0][1],self.index_lim[1][0],self.index_lim[1][1])
        for p in range (self.graph.get_nb_pos()):
            (ip,jp) = self.graph.pos_in_grid(p)
            offset_i = int(floor( 2 * self.radius / self.step))
            offset_i =  min( offset_i, (imax - ip))
            for k in range( offset_i+1 ):
                i = ip + k
                tmp = k*k*self.step*self.step
                tmp = sqrt(4*self.radius*self.radius - tmp)
                j1 = min (int(floor((jp*self.step + tmp)/self.step)), jmax)
                j0 = max (int(ceil((jp*self.step - tmp)/self.step)), jmin)
                for j in range(j0,j1+1):
                    index = self.grid[i,j]
                    if index > -1:
                        self.graph.add_edge_pos(index,p)


    """
    Creates the graph.
    """
    def create_graph(self):
        self.get_striking_shots()
        self.create_grid()
        for shot in self.striking_shots:
            interceptable = self.add_intercept(shot)
            if not interceptable:
                return False
        self.treat_overlapping()
        return True



    """
    Exact Algorithm :
    Finds an independent set of pos-vertices that dominates the shot-vertices
    Returns true and sets the solution with the index of the vertices in the solution if it finds it
         or false and sets the solution empty otherwise
    """
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
