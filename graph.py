import math

class graph :
    def __init__(self):
        self.adj_pos    = []
        self.shot       = []

    """
    Add position and create edge between existing pos if necessary
    """
    def add_pos(i, j, id_shot, d):
        insert_pos(i,j)

    """
    Add shot and create edge between pos if necessary
    """
    def add_shot():
        pos.add(len(pos))
        return len(pos)+1

    """
    Add edge between a shot and a pos
    """
    def add_edgeShot(id_shot, i, j):
        index = find(i, j)
        adj_pos[index][1].add(id_shot)

    """
    Check if two positions are overlapping
    i,j first pos
    k,l second pos
    r radius of the robot
    """
    def check_near(i, j, k, l, r):
        dist = math.sqrt((k - i)**2 + (l - j)**2)
        if(dist < r):
            return False
        return True

    """
    Add edge between two shots
    """
    def add_edgePos(i, j, k ,l):
        index1 = find(i, j)
        index2 = find(k, l)
        adj_pos[index1][2].add( (k, l) )
        adj_pos[index2][2].add( (i, j) )

    """
    Find index of (i, j) pos
    """
    def find(i, j):
        for k in range(len(adj_pos)):
            if( adj_pos[k][0] == (i, j)): #Tuple comparaison
                return k
        return -1 #not in list

    """
    Add pos with empty value for other element
    """
    def insert_pos(i, j):
        if find(i, j) == -1:
            adj_pos.add([(i,j), [], []])
