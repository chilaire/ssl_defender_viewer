import math

class graph :
    def __init__(self):
        self.adj_pos    = []
        self.shot       = []

    """
    Adds a position
    return true if we need to check overlapping
    """
    def add_pos(i, j, id_shot):
        shot[id_shot].append((i,j))
        (found,index) = find(i, j)
        if not found:
            adj_pos.insert(i,j)
        adj_pos[index][1].append(id_shot)
        return (not found, index)

    """
    Creates a shot and returns its index (no edge yet)
    """
    def add_shot(): #modified
        shot.append([])
        return len(shot)-1

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
    def add_edgePos(index1, index2, i, j, k, l):
        adj_pos[index1][2].append( (k, l) )
        adj_pos[index2][2].append( (i, j) )

    """
    Find index of (i, j) pos
    """
    def find(i, j):
        start = 0
        end = len(adj_pos)-1
        if adj_pos[start] == (i,j):
            return (True, start)
        if adj_pos[end] == (i,j):
            return (True, end)
        while start-end > 1:
            middle = (end-start)//2
            if adj_pos[middle][0] == (i,j):
                return (True, middle)
            else:
                if adj_pos[middle][0][0] < i:
                    start = middle
                elif adj_pos[middle][0][0] > i:
                    end = middle
                elif adj_pos[middle][0][1] < j:
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

    def insert_pos(i, j):# modified Ben alors, la liste est triée, non mais
        adj_pos.append([(i,j), [], []])
        return len(adj_pos)-1
    """

    def add_pos_adja(ind_pos, pos, neighbours):
        for (k,l) in neighbours:
            (found,index) = find(i,j)
            if found:
                add_edgePos(ind_pos, index, pos[0], pos[1], k, l)



""" Est-ce que par hasard on ne devrait pas stocker ce qu'on enlève ? Histoire de pouvoir remonter...
    def remove_vertex_and_neighbours(i,j):
        (found,index) = find(i,j)
        if found:
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
                self.adj_pos.pop(index2) # /!\ opération linéaire
            for shot in self.adj_pos[index][1] :
                for (kk,ll) in self.shot[shot] :
                    #supprimer l'arrete (kk,ll)->shot
                    (found4, index4) = find(kk, ll)
                    self.adj_pos[index4][2] = [(x,y) for (x,y) in self.adj_pos[index4][2] if (x!=kk and y!=ll)]
                #supprimer ligne shot
                self.shot.pop(shot) /!\ ARG
            #supprimer ligne (i,j)
            self.adj_pos.pop(index)

"""
