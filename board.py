import numpy
import pygame

from problem import *
from solution import *

class Board :
    def __init__(self, problem, solution):
        self.problem = problem
        self.solution = solution
        self.size = numpy.array([1280,960])
        self.goal_thickness = 5
        # colors
        self.background_color = (0,0,0)
        self.opponent_color = (255,0,255)
        self.defender_color = (255,255,0)
        self.goal_color = (255,255,255)
        self.success_color = (0,255,0)
        self.failure_color = (255,0,0)

    """ Return the position of the center of the image """
    def getImgCenter(self):
        return self.size / 2

    """ Return the ratio between image and field size [px/m]"""
    def getRatio(self):
        return 0.95 * min(self.size[0] / self.problem.getFieldWidth(),
                          self.size[1] / self.problem.getFieldHeight())

    """ From field referential to img position """
    def getPixelFromField(self, pos_in_field):
        ratio = self.getRatio();
        offset_field = pos_in_field - self.problem.getFieldCenter()
        offset_pixel = self.getRatio() * offset_field
        offset_pixel[1] *= -1 # Y axis is inverted to get the Z-axis pointing outside of the screen
        pixel = self.getImgCenter() + offset_pixel
        return [int(pixel[0]), int(pixel[1])]

    def drawSegmentInField(self, screen, color, pos1, pos2, thickness):
        start = self.getPixelFromField(pos1)
        end = self.getPixelFromField(pos2)
        pygame.draw.line(screen, color, start, end, thickness)
    
    def drawKickRay(self, screen, robot_pos, kick_dir):
        # Getting closest goal to score
        kick_end = None
        best_dist = None
        for goal in self.problem.goals:
            kick_result = goal.kickResult(robot_pos, kick_dir)
            if not kick_result is None:
                goal_dist = numpy.linalg.norm(robot_pos - kick_result)
                if best_dist == None or goal_dist < best_dist:
                    best_dist = goal_dist
                    kick_end = kick_result
        if not kick_end is None:
            # Checking if kick is intercepted by one of the opponent and which one is the first
            intercepted = False
            for def_id in range(self.solution.getNbDefenders()):
                defender = self.solution.getDefender(def_id)
                collide_point = segmentCircleIntersection(robot_pos, kick_end,
                                                          defender, self.problem.robot_radius)
                if not collide_point is None:
                    kick_end = collide_point
                    intercepted = True
            ##### TODO
            color = self.failure_color
            if intercepted:
                color = self.success_color
            self.drawSegmentInField(screen, color, robot_pos, kick_end, 1)

    def drawKickRays(self, screen):
        for opp_id in range(self.problem.getNbOpponents()):
            kick_dir = 0
            while kick_dir < 2 * math.pi:
                self.drawKickRay(screen, self.problem.getOpponent(opp_id), kick_dir)
                kick_dir += self.problem.theta_step
            

    def drawGoals(self, screen):
        for goal in self.problem.goals:
            self.drawSegmentInField(screen, self.goal_color,
                                    goal.posts[:,0], goal.posts[:,1],
                                    self.goal_thickness)
    
    def drawRobots(self, screen, robots, color):
        for robot_id in range(robots.shape[1]):
            pygame.draw.circle(screen, color,
                               self.getPixelFromField(robots[:, robot_id]),
                               int(self.problem.robot_radius * self.getRatio()))
        
    def drawOpponents(self, screen):
        self.drawRobots(screen, self.problem.opponents, self.opponent_color)

    def drawDefenders(self, screen):
        self.drawRobots(screen, self.solution.defenders, self.defender_color)
        
    def draw(self, screen):
        self.drawKickRays(screen)
        self.drawGoals(screen)
        self.drawOpponents(screen)
        self.drawDefenders(screen)

    def run(self):
        pygame.init()
        screen = pygame.display.set_mode(self.size)
        running = True
        while running:
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    running=False
            keys=pygame.key.get_pressed()
            if (keys[pygame.K_ESCAPE]): running = False
        
            screen.fill(self.background_color)
            self.draw(screen);
            pygame.display.flip()
        
