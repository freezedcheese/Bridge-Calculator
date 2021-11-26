import math
import turtle

def find_hyp(x, y):
    return math.sqrt(x**2 + y**2)

class Vector():
    def __init__(self, x_dir, y_dir, mag=1):
        self.x_dir = x_dir
        self.y_dir = y_dir
        self.mag = mag

    def get_opposite(self):
        return Vector(-self.x_dir, -self.y_dir, self.mag)
    
    def get_mag(self):
        return self.mag
    
    def get_x_mag(self):
        return self.mag * (self.x_dir / find_hyp(self.x_dir, self.y_dir))
    
    def get_y_mag(self):
        return self.mag * (self.y_dir / find_hyp(self.x_dir, self.y_dir))
    
    def set_mag(self, mag):
        self.mag = mag
    
    def set_mag_by_comp(self, x_mag_comp, y_mag_comp):
        if self.x_dir != 0:
            self.mag = x_mag_comp * (find_hyp(self.x_dir, self.y_dir) / self.x_dir)
        elif self.y_dir != 0:
            self.mag = y_mag_comp * (find_hyp(self.x_dir, self.y_dir) / self.y_dir)
        else:
            self.mag = 0

class TrussNode():
    def __init__(self, loc=(0,0)):
        self.loc = loc
        self.neighbours = []
        self.external_forces = []
        self.internal_forces = {}

    def new_neighbour_rel(self, rel):
        neighbour_loc = (self.loc[0] + rel[0], self.loc[1] + rel[1])
        self.neighbours.append(TrussNode(neighbour_loc))
        self.neighbours[-1].neighbours.append(self)

        return self.neighbours[-1]
    
    def add_neighbour(self, neighbour):
        self.neighbours.append(neighbour)
        neighbour.neighbours.append(self)
    
    def add_external_force(self, force):
        self.external_forces.append(force)
    
    def calc_rel_pos(self, neighbour):
        return Vector(neighbour.loc[0] - self.loc[0], neighbour.loc[1] - self.loc[1])

    def calc_internal_forces(self, debug_print=False):
        unknown_internal_forces = []

        f_x_tot, f_y_tot = 0, 0

        for neighbour in self.neighbours:
            if self in neighbour.internal_forces:
                self.internal_forces[neighbour] = neighbour.internal_forces[self].get_opposite()

                f_x_tot += self.internal_forces[neighbour].get_x_mag()
                f_y_tot += self.internal_forces[neighbour].get_y_mag()
            else:
                unknown_internal_forces.append(neighbour)
                self.internal_forces[neighbour] = Vector(neighbour.loc[0] - self.loc[0], neighbour.loc[1] - self.loc[1])
        
        for force in self.external_forces:
            f_x_tot += force.get_x_mag()
            f_y_tot += force.get_y_mag()

        if len(unknown_internal_forces) > 2:
            print("ERROR: node at", self.loc, "contains more than 2 unknown forces")
            return
        
        elif len(unknown_internal_forces) == 2:
            f_a = self.internal_forces[unknown_internal_forces[0]]
            f_b = self.internal_forces[unknown_internal_forces[1]]

            if f_a.get_x_mag() == 0:
                self.internal_forces[unknown_internal_forces[1]].set_mag(-(f_b.get_mag()/f_b.get_x_mag()) * f_x_tot)
                if f_b.get_mag() == 0:
                    self.internal_forces[unknown_internal_forces[0]].set_mag_by_comp(-f_x_tot, -f_y_tot)
                else:
                    self.internal_forces[unknown_internal_forces[0]].set_mag(-f_y_tot - self.internal_forces[unknown_internal_forces[1]].get_mag() * (f_b.get_y_mag()/f_b.get_mag()))
            elif f_a.get_y_mag() == 0:
                self.internal_forces[unknown_internal_forces[1]].set_mag(-(f_b.get_mag()/f_b.get_y_mag()) * f_y_tot)
                if f_b.get_mag() == 0:
                    self.internal_forces[unknown_internal_forces[0]].set_mag_by_comp(-f_x_tot, -f_y_tot)
                else:
                    self.internal_forces[unknown_internal_forces[0]].set_mag(-f_x_tot - self.internal_forces[unknown_internal_forces[1]].get_mag() * (f_b.get_x_mag()/f_b.get_mag()))
            else:
                self.internal_forces[unknown_internal_forces[1]].set_mag((f_b.get_mag() * (f_y_tot/f_a.get_y_mag() - f_x_tot/f_a.get_x_mag())) / (f_b.get_x_mag()/f_a.get_x_mag() - f_b.get_y_mag()/f_a.get_y_mag()))
                self.internal_forces[unknown_internal_forces[0]].set_mag((f_a.get_mag()/f_a.get_x_mag()) * (-f_x_tot - self.internal_forces[unknown_internal_forces[1]].get_mag() * (f_b.get_x_mag()/f_b.get_mag())))
        
        elif len(unknown_internal_forces) == 1:
            self.internal_forces[unknown_internal_forces[0]].set_mag_by_comp(-f_x_tot, -f_y_tot)

        if debug_print:
            for thing in self.neighbours:
                print(self.internal_forces[thing].get_mag())

class Truss():
    def __init__(self, flattop, csa, young_mod):
        '''
        flattop: Flattop (bool)
        csa: Cross-Sectional Area (mm^2)
        yound_mod: Young's Modulus (MPa)
        '''
        self.flattop = flattop
        self.csa = csa
        self.young_mod = young_mod

        self.truss_list = []
        self.lengths = {}
        self.internal_forces = {}
        self.deformations = {}

    def gen_warren_truss(self, triangle_count, diagonal):
        self.truss_list = []
        self.truss_list.append(TrussNode())

        for i in range(triangle_count):
            self.truss_list.append(self.truss_list[-1].new_neighbour_rel(diagonal))
            self.truss_list.append(self.truss_list[-1].new_neighbour_rel((diagonal[0], -diagonal[1])))

        for i in range(0,triangle_count*2,2):
            self.truss_list[i].add_neighbour(self.truss_list[i+2])

        for i in range(1,triangle_count*2-1,2):
            self.truss_list[i].add_neighbour(self.truss_list[i+2])

        if self.flattop:
            self.truss_list.append(self.truss_list[-1].new_neighbour_rel((0,diagonal[1])))
            self.truss_list.append(self.truss_list[0].new_neighbour_rel((0,diagonal[1])))

            self.truss_list[-2].add_neighbour(self.truss_list[-4])
            self.truss_list[-1].add_neighbour(self.truss_list[1])

        for node in self.truss_list:
            for neighbour in node.neighbours:
                key = tuple(sorted([self.truss_list.index(node), self.truss_list.index(neighbour)]))

                self.lengths[key] = find_hyp(node.calc_rel_pos(neighbour).x_dir, node.calc_rel_pos(neighbour).y_dir)

    def update_forces(self, force_dict):
        for key in force_dict:
            self.truss_list[key].add_external_force(force_dict[key])

        if self.flattop:
            self.truss_list[-1].calc_internal_forces()
            self.truss_list[-2].calc_internal_forces()

        for node in self.truss_list:
            node.calc_internal_forces()

        for node in self.truss_list:
            for neighbour in node.internal_forces:
                key = tuple(sorted([self.truss_list.index(node), self.truss_list.index(neighbour)]))

                if node.calc_rel_pos(neighbour).get_x_mag() != 0:
                    self.internal_forces[key] = abs(node.internal_forces[neighbour].get_mag()) * (-1 if node.internal_forces[neighbour].get_x_mag() / node.calc_rel_pos(neighbour).get_x_mag() < 0 else 1)
                else:
                    self.internal_forces[key] = abs(node.internal_forces[neighbour].get_mag()) * (-1 if node.internal_forces[neighbour].get_y_mag() / node.calc_rel_pos(neighbour).get_y_mag() < 0 else 1)
    
    def calc_deformations(self):
        for key in self.internal_forces:
            self.deformations[key] = ((self.internal_forces[key] / self.csa) / self.young_mod) * self.lengths[key]
    
    def dist_loads(self):
        for node in self.truss_list:
            for neighbour in node.neighbours:
                pass
        
    def draw_truss(self, colour=False):
        SCALE = 15
        COLOUR_MULT = 0.5
        turtle.speed(0)
        turtle.colormode(255)
        turtle.Screen().tracer(0)

        max_x = 0
        for node in self.truss_list:
            max_x = node.loc[0] if node.loc[0] > max_x else max_x
        
        x_offset = max_x * SCALE / 2

        for node in self.truss_list:
            for neighbour in node.neighbours:
                turtle.penup()
                if colour:
                    if node.calc_rel_pos(neighbour).get_x_mag() != 0:
                        RED_MULT = COLOUR_MULT if node.internal_forces[neighbour].get_x_mag() / node.calc_rel_pos(neighbour).get_x_mag() < 0 else 0
                        BLUE_MULT = COLOUR_MULT if node.internal_forces[neighbour].get_x_mag() / node.calc_rel_pos(neighbour).get_x_mag() > 0 else 0
                    else:
                        RED_MULT = COLOUR_MULT if node.internal_forces[neighbour].get_y_mag() / node.calc_rel_pos(neighbour).get_y_mag() < 0 else 0
                        BLUE_MULT = COLOUR_MULT if node.internal_forces[neighbour].get_y_mag() / node.calc_rel_pos(neighbour).get_y_mag() > 0 else 0
                    turtle.pencolor(int(abs(node.internal_forces[neighbour].get_mag()) * RED_MULT), 0, int(abs(node.internal_forces[neighbour].get_mag()) * BLUE_MULT))
                turtle.goto((node.loc[0] * SCALE - x_offset, node.loc[1] * SCALE))
                turtle.pendown()
                turtle.goto((neighbour.loc[0] * SCALE - x_offset, neighbour.loc[1] * SCALE))
        
        turtle.Screen().update()
        
        turtle.mainloop()

truss = Truss(True, 1,1)
truss.gen_warren_truss(16, (3,4))

truss_forces = {0: Vector(0,1,40), 32: Vector(0,1,40), 16: Vector(0,-1,80)}

truss.update_forces(truss_forces)
truss.draw_truss(True)