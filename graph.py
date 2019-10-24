import math

class graph :
    def __init__(self):
        self.adj_pos    = []
        self.shot       = []

    """
    Adds a position
    return true if we need to check overlapping
    """
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
    Add edge between a shot and a pos

    def add_edgeShot(id_shot, i, j): #en fait pas besoin
        index = find(i, j)
        adj_pos[index][1].add(id_shot)
    """

    """
    Check if two positions are overlapping
    i,j first pos
    k,l second pos
    r radius of the robot

    def check_near(i, j, k, l, r):
        dist = math.sqrt((k - i)**2 + (l - j)**2)
        if(dist < r):
            return False
        return True
    """
    """
    Add edge between two shots
    """
    def add_edgePos(self,index1, index2, i, j, k, l):
        self.adj_pos[index1][2].append( (k, l) )
        self.adj_pos[index2][2].append( (i, j) )

    """
    Find index of (i, j) pos
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
        """for k in range(len(adj_pos)):
            if( adj_pos[k][0] == (i, j)): #Tuple comparaison
                return k
        return -1 #not in list"""

    """
    Add pos with empty value for other element

    def insert_pos(i, j):# modified Ben alors, la liste est triee, non mais
        adj_pos.append([(i,j), [], []])
        return len(adj_pos)-1
    """
#ici on a des index pour les pos
    def add_pos_adja(self,ind_pos, pos, neighbours): #TODO: faire sans liste
        for (i,j) in neighbours:
            (found,index) = self.find(i,j)
            #print("pos : ", ind_pos, ", ", pos, " ; vois : ", index, ", ", i, ", ", j, "\n\n")
            if found:
                self.add_edgePos(ind_pos, index, pos[0], pos[1], i, j)

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



    def get_first_shot(self):
        for shot in self.shot:
            if (shot[1] < 0) :
                return [id_pos for id_pos in shot[0] if (self.adj_pos[id_pos][3] < 0)]
        return None



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

    def revive_vertex_and_neighbours(self,index, turn):
        self.adj_pos[index][3] = -1
        for p in range(len(self.adj_pos[index][2])):
            pos = self.adj_pos[p]
            if (pos[3] == turn) :
                pos[3] = -1
        for s in range(len(self.adj_pos[index][1])) :
            shot = self.shot[s]
            if shot[1] == turn :
                shot[1] = -1

'''
# Est-ce que par hasard on ne devrait pas stocker ce qu'on enleve ? Histoire de pouvoir remonter...
    def remove_vertex_and_neighbours(index):
            for (k,l) in self.adj_pos[index][2] :
                (found2,index2)=find(k,l)
                for (kk,ll) in self.adj_pos[index2][2] :
                    #supprimer l'arrete (kk,ll)->(k,l) (pas opti)
                    (found3,index3)=find(kk,ll)
                    self.adj_pos[index2][2] = [(x,y) for (x,y) in self.adj_pos[index2][2] if (x!=k and y!=l)]
                for shot in self.adj_pos[index2][1] :
                    #supprimer l'arrete shot->(k,l)
                    self.shot[shot] = [(x,y) for (x,y) in self.shot[shot] if (x!=k and y!=l)]
                #supprimer ligne (k,l)
                self.adj_pos.pop(index2) # /!\ operation lineaire
            for shot in self.adj_pos[index][1] :
                for (kk,ll) in self.shot[shot] :
                    #supprimer l'arrete (kk,ll)->shot
                    (found4, index4) = find(kk, ll)
                    self.adj_pos[index4][2] = [(x,y) for (x,y) in self.adj_pos[index4][2] if (x!=kk and y!=ll)]
                #supprimer ligne shot
                self.shot.pop(shot)
            #supprimer ligne (i,j)
            self.adj_pos.pop(index)
'''
