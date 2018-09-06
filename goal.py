import math
import numpy
import sys

class Goal:
    def __init__(self, data):
        mandatory_keys = ["posts", "direction"]
        for key in mandatory_keys:
            if key not in data:
                raise ValueError("Cannot find '" + key + "'")
        self.posts = numpy.array(data["posts"]).transpose()
        self.direction = numpy.array(data["direction"])
        if (self.posts.shape != (2,2)):
            raise ValueError("Invalid shape for 'posts': "
                             + str(self.field_limits.shape) + " expecting (2, 2)")
        if (self.direction.shape != (2,)):
            raise ValueError("Invalid shape for 'direction': "
                             + str(self.direction.shape) + " expecting (2,)")

    """ Return [shotOnTarget, finalPos] for a kick from 'pos' with direction 'theta' """
    def kickResult(self, pos, theta):
        kick_dir = numpy.array([math.cos(theta),math.sin(theta)])
        if (kick_dir.dot(self.direction) <= 0):
            return [False,pos]
        # Line-segment intersection from points
        # https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection#Given_two_points_on_each_line
        x1 = self.posts[0,0]
        x2 = self.posts[0,1]
        y1 = self.posts[1,0]
        y2 = self.posts[1,1]
        x3 = pos[0]
        x4 = x3 + kick_dir[0]
        y3 = pos[1]
        y4 = y3 + kick_dir[1]
        t = ((x1-x3)*(y3-y4)-(y1-y3)*(x3-x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4))
        # Intersection is between two points
        return [0 <= t <= 1, (x1+t*(x2-x1), y1+t*(y2-y1))]
        
