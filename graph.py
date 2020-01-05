import math

class graph :
    def __init__(self):
        self.adj_pos    = [] #list of [(i,j),[shot neighbourhood (index)],[pos neighbourhood (index)],"removed at round k" (bool)]
        self.shot       = [] #list of [[pos neighbourhood (index)],"removed at round k"(bool)]
        self.goalkeeper_position = [] #lost of pos_vertex (index) that are in the goalkeeper area
        self.goalkeeper_area = None #indexes min max of the goalkeeper_area

    """
    Returns the number of position_verticices in the graph
    """
    def get_nb_pos(self):
        return len(self.adj_pos)

    """
    Return the grid_position of a position_vertex
    """
    def pos_in_grid(self,index):
        return self.adj_pos[index][0]

    """
    Create a position_vertex in the graph's adjacency list if it is not already in the graph
    return the index of the pos in the graph (no eges yet)
    """
    def add_pos(self,i, j, idShot):
        self.adj_pos.append([(i,j), [], [], -1])
        return len(self.adj_pos)-1

    """
    Creates a shot_vertex and returns its index (no edge yet)
    """
    def add_shot(self,):
        self.shot.append([[], -1])
        return len(self.shot)-1


    """
    Add edge between two position_vertex
    """
    def add_edge_pos(self,index1, index2):
        self.adj_pos[index1][2].append( index2 )
        self.adj_pos[index2][2].append( index1 )


    """
    Add edge between a position_vetex and a shot_vertex
    """
    def add_edge_shot(self,indexPos, idShot):
        self.adj_pos[indexPos][1].append( idShot )
        self.shot[idShot][0].append( indexPos )



    """
    Gets the neighbourhood (position_vertex not removed) of the first not removed shot.
    Returns the neighbourhood or None if no shot is found
    """
    def get_first_shot_neighbourhood(self):
        for shot in self.shot:
            if (shot[1] < 0) :
                return [id_pos for id_pos in shot[0] if (self.adj_pos[id_pos][3] < 0)]
        return None

    """
    Gets the position_vertex with the most shot-neighbours, neighbour of the shot_vertex with the least neighbours
    (considering only non removed elements).
    Returns (position,True) if the position is found
            (None, True) if a shot with no neighbour is found
            (None, False) if no shot is found
    """
    def get_heuristic_pos(self):
        nb_nei_min = -1
        best_shot = None
        for shot in self.shot:           #gets "best_shot" = the shot with the least neighbor
            if shot[1] == -1:
                nb_nei = 0
                pos = None
                for n in shot[0]:
                    if self.adj_pos[n][3] == -1:
                        nb_nei += 1
                        pos = n
                if nb_nei == 0:         #if 0 neighbour -> game over
                    return (None, True)
                elif nb_nei == 1:       #if 1 neighbour -> take this one
                    return (pos, True)
                elif  nb_nei < nb_nei_min or nb_nei_min == -1 :
                    best_shot = shot
                    nb_nei_min = nb_nei
        if best_shot == None :          #gets the neighbour of best_shot covering the most shots
            return (None, False)
        best_pos = None
        nb_shots_max = 0
        for n in best_shot[0]:
            if self.adj_pos[n][3] == -1:
                nb_shots = 0
                for s in self.adj_pos[n][1]:
                    if self.shot[s][1] == -1:
                        nb_shots += 1
                if nb_shots > nb_shots_max:
                    best_pos = n
                    nb_shots_max = nb_shots
        return (best_pos, True)

    """
        Removes the pos-vertex at the index "index" and all its pos-neighbours at the turn "turn"
        Sets the component 3 (which says if the vertex is dead or alive) of the vertices at turn if they are alive
    """
    def remove_vertex_and_neighbours(self,index, turn):
        self.adj_pos[index][3] = turn           #Removes vertex
        for s in self.adj_pos[index][1] :       #Removes shot-vertex neighbour
            shot = self.shot[s]
            if shot[1] < 0 :
                shot[1] = turn
        for p in self.adj_pos[index][2]:       #Removes pos-vertex neighbour
            pos = self.adj_pos[p]
            if pos[3] < 0 :
                pos[3] = turn
        if not self.goalkeeper_area is None :  #Removes pos-vertex "neighbour" from goalkeeper_area
            (i,j) = self.adj_pos[index][0]
            ((i0,i1),(j0,j1)) = self.goalkeeper_area
            if i>=i0 and i<=i1 and j>=j0 and j<=j1 :
                for p in self.goalkeeper_position :
                    if self.adj_pos[p][3] < 0 :
                        self.adj_pos[p][3] = turn


    """
        Revives the pos-vertex at the index "index" and all its pos-neighbours removed at the turn "turn"
        Sets the component 3 (which says if the vertex is dead or alive) of the vertices at -1 if it has been kill at turn turn
    """
    def revive_vertex_and_neighbours(self,index, turn):
        self.adj_pos[index][3] = -1           #revives vertex
        for s in self.adj_pos[index][1] :       #revives shot-vertex neighbour
            shot = self.shot[s]
            if shot[1] == turn :
                shot[1] = -1
        for p in self.adj_pos[index][2] :       #revives pos-vertex neighbour
            pos = self.adj_pos[p]
            if pos[3] == turn :
                pos[3] = -1
        if not self.goalkeeper_area is None :  #revives pos-vertex "neighbour" from goalkeeper_area
            (i,j) = self.adj_pos[index][0]
            ((i0,i1),(j0,j1)) = self.goalkeeper_area
            if i>=i0 and i<=i1 and j>=j0 and j<=j1 :
                for p in self.goalkeeper_position :
                    if self.adj_pos[p][3] == turn :
                        self.adj_pos[p][3] = -1
