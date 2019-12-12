import math

class graph :
    def __init__(self):
        self.adj_pos    = [] #liste de [(i,j),[shot neighbourhood (index)],[pos neighbourhood ((i',j') then index)],"removed at round k"]
        self.shot       = []#liste de [[pos neighbourhood ((i',j') then index)],"removed at round k"]


        #########################################################
        ### [pos neighbourhood] is a list of position (i',j') ###
        #########################################################

    def get_nb_pos(self):
        return len(self.adj_pos)

    """
    Adds a position in the sorted grah's adjacency list if it is not already in the graph
    return true if we need to treat the overlapping, false otherwise
    and the index of the pos in the graph
    """
    #if the position was already in the graph, the overlapping has already been treated
    def add_pos(self,i, j, id_shot):
        self.shot[id_shot][0].append((i,j))
        (found,index) = self.find(i, j)
        if not found:
            self.adj_pos.insert(index, [(i,j), [], [], -1])
        self.adj_pos[index][1].append(id_shot)
        return (not found, index)

    """
    Creates a shot and returns its index (no edge yet)
    """
    def add_shot(self,): #modified
        self.shot.append([[], -1])
        return len(self.shot)-1


    """
    Add edge between two positions
    """
    def add_edgePos(self,index1, index2, i, j, k, l):
        self.adj_pos[index1][2].append( (k, l) )
        self.adj_pos[index2][2].append( (i, j) )

    """
    Find the index of the pos (i,j)
    Returns true and its index if the pos is already in the graph
            false and the position where it should be inserted
    """
    def find(self, i, j):
        if self.adj_pos == []:
            return (False, 0)
        start = 0
        end = len(self.adj_pos)-1
        if self.adj_pos[end][0] < (i,j):
            return (False, end+1)
        if self.adj_pos[start][0] == (i,j):
            return (True, start)
        if self.adj_pos[end][0] == (i,j):
            return (True, end)
        while end-start > 1:
            middle = start + (end-start)//2
            if self.adj_pos[middle][0] == (i,j):
                return (True, middle)
            else:
                if self.adj_pos[middle][0][0] < i:
                    start = middle
                elif self.adj_pos[middle][0][0] > i:
                    end = middle
                elif self.adj_pos[middle][0][1] < j:
                    start = middle
                else:
                    end = middle
        return (False, end)

    """
    Add an edge between a position and each potentiel neighbours
    """
    def add_pos_adja(self,ind_pos, pos, neighbours): #TODO: faire sans liste
        for (i,j) in neighbours:
            (found,index) = self.find(i,j)
            if found:
                self.add_edgePos(ind_pos, index, pos[0], pos[1], i, j)





    """
    Converts the adjacency lists of position (initially list of (i',j')) in a list of index
    """
    def convert(self):
        for index in range(len(self.adj_pos)):
            (i,j) = self.adj_pos[index][0]
            for p in range(len(self.adj_pos)):
                pos = self.adj_pos[p]
                for k in range(len(pos[2])):
                    if pos[2][k] == (i,j):
                        pos[2][k] = index
            for s in self.adj_pos[index][1]:
                shot = self.shot[s]
                for k in range(len(shot[0])):
                    if shot[0][k] == (i,j):
                        shot[0][k] = index






        #########################################################
        ### [pos neighbourhood] is a list of index of position###
        #########################################################

    """
    Gets the neighbourhood (position not removed) of the first not removed shot.
    Returns the neighbourhood or None if no shot is found
    """
    def get_first_shot_neighbourhood(self):
        for shot in self.shot:
            if (shot[1] < 0) :
                return [id_pos for id_pos in shot[0] if (self.adj_pos[id_pos][3] < 0)]
        return None

    """
    Gets the position with the more shot-neighbours, neighbour of the shot with fewer neighbours
    (considering only non removed elements).
    Returns (position,True) if the position is found
            (None, True) if a shot with no neighbour is found
            (None, False) if no shot is found
    """
    def get_heuristic_pos(self):
        nb_pos_min = math.inf
        best_shot = None
        for shot in self.shot:
            if shot[1] == -1:
                len_nei = 0
                pos = None
                for n in shot[0]:
                    if self.adj_pos[n][3] == -1:
                        len_nei += 1
                        pos = n
                if len_nei == 0:
                    return (None, True)
                elif len_nei == 1:
                    return (pos, True)
                elif len_nei < nb_pos_min:
                    best_shot = shot
        if best_shot == None :
            return (None, False)
        best_pos = None
        nb_shots_max = 0
        for n in best_shot[0]:
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
        self.adj_pos[index][3] = turn
        for p in self.adj_pos[index][2]:
            pos = self.adj_pos[p]
            if (pos[3] < 0) :
                pos[3] = turn
        for s in self.adj_pos[index][1] :
            shot = self.shot[s]
            if shot[1] < 0 :
                shot[1] = turn


    """
        Removes the pos-vertex at the index "index" and all its pos-neighbours at the turn "turn"
        Sets the component 3 (which says if the vertex is dead or alive) of the vertices at -1 if it has been kill at turn turn
    """
    def revive_vertex_and_neighbours(self,index, turn):
        self.adj_pos[index][3] = -1
        for p in self.adj_pos[index][2] :
            pos = self.adj_pos[p]
            if (pos[3] == turn) :
                pos[3] = -1
        for s in self.adj_pos[index][1] :
            shot = self.shot[s]
            if shot[1] == turn :
                shot[1] = -1
