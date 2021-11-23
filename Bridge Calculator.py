import math
import turtle

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
    
    def set_mag_by_x_mag(self, x_mag_comp):
        self.mag = x_mag_comp * (find_hyp(self.x_dir, self.y_dir) / self.x_dir)
    
    def set_mag_by_y_mag(self, y_mag_comp):
        self.mag = y_mag_comp * (find_hyp(self.x_dir, self.y_dir) / self.y_dir)

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

        self.internal_forces[self.neighbours[-1]] = []

        return self.neighbours[-1]
    
    def add_neighbour(self, neighbour):
        self.neighbours.append(neighbour)
        neighbour.neighbours.append(self)

        self.internal_forces[self.neighbours[-1]] = []
    
    def add_external_force(self, force):
        self.external_forces.append(force)

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
                self.internal_forces[unknown_internal_forces[0]].set_mag(-f_y_tot - self.internal_forces[unknown_internal_forces[1]].get_mag() * (f_b.get_y_mag()/f_b.get_mag()))
            elif f_a.get_y_mag() == 0:
                self.internal_forces[unknown_internal_forces[1]].set_mag(-(f_b.get_mag()/f_b.get_y_mag()) * f_y_tot)
                self.internal_forces[unknown_internal_forces[0]].set_mag(-f_x_tot - self.internal_forces[unknown_internal_forces[1]].get_mag() * (f_b.get_x_mag()/f_b.get_mag()))
            else:
                self.internal_forces[unknown_internal_forces[1]].set_mag((f_b.get_mag() * (f_y_tot/f_a.get_y_mag() - f_x_tot/f_a.get_x_mag())) / (f_b.get_x_mag()/f_a.get_x_mag() - f_b.get_y_mag()/f_a.get_y_mag()))
                self.internal_forces[unknown_internal_forces[0]].set_mag((f_a.get_mag()/f_a.get_x_mag()) * (-f_x_tot - self.internal_forces[unknown_internal_forces[1]].get_mag() * (f_b.get_x_mag()/f_b.get_mag())))
        
        elif len(unknown_internal_forces) == 1:
            self.internal_forces[unknown_internal_forces[0]].set_mag(find_hyp(f_x_tot, self.f_y_tot))

        if debug_print:
            for thing in self.neighbours:
                print(self.internal_forces[thing].get_x_mag(), self.internal_forces[thing].get_y_mag(), "tension" if self.internal_forces[thing].get_x_mag() / (thing.loc[0] - self.loc[0]) > 0 else "compression")

def find_hyp(x, y):
    return math.sqrt(x**2 + y**2)

def draw_truss(truss):
    SCALE = 20
    turtle.speed(10)
    for node in truss:
        for neighbour in node.neighbours:
            turtle.penup()
            turtle.goto((node.loc[0] * SCALE, node.loc[1] * SCALE))
            turtle.pendown()
            turtle.goto((neighbour.loc[0] * SCALE, neighbour.loc[1] * SCALE))
        
    turtle.mainloop()

#simple truss
node_1 = TrussNode((0,0))
node_2 = node_1.new_neighbour_rel((6,0))
node_3 = node_2.new_neighbour_rel((6,0))
node_4 = node_3.new_neighbour_rel((6,0))
node_5 = node_4.new_neighbour_rel((6,0))

node_6 = TrussNode((3,4))
node_7 = node_6.new_neighbour_rel((6,0))
node_8 = node_7.new_neighbour_rel((6,0))
node_9 = node_8.new_neighbour_rel((6,0))

node_1.add_neighbour(node_6)
node_6.add_neighbour(node_2)
node_2.add_neighbour(node_7)
node_7.add_neighbour(node_3)
node_3.add_neighbour(node_8)
node_8.add_neighbour(node_4)
node_4.add_neighbour(node_9)
node_9.add_neighbour(node_5)

node_1.add_external_force(Vector(0, 1, 120))
node_2.add_external_force(Vector(0, 1, -80))
node_3.add_external_force(Vector(0, 1, -80))
node_4.add_external_force(Vector(0, 1, -80))
node_5.add_external_force(Vector(0, 1, 120))

node_1.calc_internal_forces()
node_6.calc_internal_forces()
node_2.calc_internal_forces()
node_7.calc_internal_forces()
node_3.calc_internal_forces(True)

truss = [node_1, node_2, node_3, node_4, node_5, node_6, node_7, node_8, node_9]

draw_truss(truss)