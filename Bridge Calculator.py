import math
import turtle
import matplotlib.pyplot as plt
import numpy as np

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
        if find_hyp(self.x_dir, self.y_dir) == 0:
            return 0
        else:
            return self.mag * (self.x_dir / find_hyp(self.x_dir, self.y_dir))
    
    def get_y_mag(self):
        if find_hyp(self.x_dir, self.y_dir) == 0:
            return 0
        else:
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
    
    def set_mag_by_current_dir_comps(self):
        self.mag = find_hyp(self.x_dir, self.y_dir)
    
    def set_dir_comps_by_current_mag(self):
        self.x_dir = self.get_x_mag()
        self.y_dir = self.get_y_mag()
    
    def add_vector(self, vector):
        self.set_dir_comps_by_current_mag()
        vector.set_dir_comps_by_current_mag()

        self.x_dir += vector.x_dir
        self.y_dir += vector.y_dir

        self.set_mag_by_current_dir_comps()


class TrussNode():
    def __init__(self, loc=(0,0)):
        self.loc = loc
        self.neighbours = set()
        self.external_forces = []
        self.internal_forces = {}

    def new_neighbour_rel(self, rel):
        neighbour_loc = (self.loc[0] + rel[0], self.loc[1] + rel[1])
        neighbour = TrussNode(neighbour_loc)
        self.neighbours.add(neighbour)
        neighbour.neighbours.add(self)

        return neighbour
    
    def add_neighbour(self, neighbour):
        self.neighbours.add(neighbour)
        neighbour.neighbours.add(self)
    
    def add_external_force(self, force):
        self.external_forces.append(force)
    
    def calc_rel_pos(self, neighbour):
        return Vector(neighbour.loc[0] - self.loc[0], neighbour.loc[1] - self.loc[1])

    def calc_internal_forces(self, debug_print=False):
        unknown_internal_forces = []

        f_x_net, f_y_net = 0, 0

        for neighbour in self.neighbours:
            if self in neighbour.internal_forces:
                self.internal_forces[neighbour] = neighbour.internal_forces[self].get_opposite()

                f_x_net += self.internal_forces[neighbour].get_x_mag()
                f_y_net += self.internal_forces[neighbour].get_y_mag()
            else:
                unknown_internal_forces.append(neighbour)
                self.internal_forces[neighbour] = Vector(neighbour.loc[0] - self.loc[0], neighbour.loc[1] - self.loc[1])
        
        for force in self.external_forces:
            f_x_net += force.get_x_mag()
            f_y_net += force.get_y_mag()

        if len(unknown_internal_forces) > 2:
            print("ERROR: node at", self.loc, "contains more than 2 unknown forces")
            return
        
        elif len(unknown_internal_forces) == 2:
            f_a = self.internal_forces[unknown_internal_forces[0]]
            f_b = self.internal_forces[unknown_internal_forces[1]]

            if f_a.get_x_mag() == 0:
                self.internal_forces[unknown_internal_forces[1]].set_mag(-(f_b.get_mag()/f_b.get_x_mag()) * f_x_net)
                if f_b.get_mag() == 0:
                    self.internal_forces[unknown_internal_forces[0]].set_mag_by_comp(-f_x_net, -f_y_net)
                else:
                    self.internal_forces[unknown_internal_forces[0]].set_mag(-f_y_net - self.internal_forces[unknown_internal_forces[1]].get_mag() * (f_b.get_y_mag()/f_b.get_mag()))
            elif f_a.get_y_mag() == 0:
                self.internal_forces[unknown_internal_forces[1]].set_mag(-(f_b.get_mag()/f_b.get_y_mag()) * f_y_net)
                if f_b.get_mag() == 0:
                    self.internal_forces[unknown_internal_forces[0]].set_mag_by_comp(-f_x_net, -f_y_net)
                else:
                    self.internal_forces[unknown_internal_forces[0]].set_mag(-f_x_net - self.internal_forces[unknown_internal_forces[1]].get_mag() * (f_b.get_x_mag()/f_b.get_mag()))
            else:
                self.internal_forces[unknown_internal_forces[1]].set_mag((f_b.get_mag() * (f_y_net/f_a.get_y_mag() - f_x_net/f_a.get_x_mag())) / (f_b.get_x_mag()/f_a.get_x_mag() - f_b.get_y_mag()/f_a.get_y_mag()))
                self.internal_forces[unknown_internal_forces[0]].set_mag((f_a.get_mag()/f_a.get_x_mag()) * (-f_x_net - self.internal_forces[unknown_internal_forces[1]].get_mag() * (f_b.get_x_mag()/f_b.get_mag())))
        
        elif len(unknown_internal_forces) == 1:
            self.internal_forces[unknown_internal_forces[0]].set_mag_by_comp(-f_x_net, -f_y_net)

        if debug_print:
            for thing in self.neighbours:
                print(self.internal_forces[thing].get_mag())

class Truss():
    def __init__(self, flattop, young_mod, poisson, tension_ult, comp_ult, cs_properties, a_rxn_loc=None, b_rxn_loc=None):
        '''
        flattop: Flattop (bool)
        cs_properties: Cross-Sectional Properties
        
        young_mod: Young's Modulus (MPa)
        poisson: Poisson's Ratio

        a_rxn_loc: Pivot reaction point A_x, A_y
        b_rxn_loc: Roller reaction point B_y
        '''

        self.flattop = flattop

        self.tension_ult = tension_ult
        self.comp_ult = comp_ult

        self.cs_properties = cs_properties        
        self.width = cs_properties["width"]
        self.thickness = cs_properties["thickness"]

        self.young_mod = young_mod
        self.poisson = poisson

        self.truss_list = []
        self.lengths = {}
        self.internal_forces = {}
        self.external_forces = {}
        self.deformed_lengths = {}

        self.failures = {}

        self.a_rxn_loc = a_rxn_loc
        self.b_rxn_loc = b_rxn_loc

    def gen_periodic_warren_truss(self, triangle_count, diagonal):
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
    
    def gen_nonperiodic_warren_truss(self, relative_base_xs, height):
        self.truss_list = []
        self.truss_list.append(TrussNode())

        for x in relative_base_xs:
            self.truss_list.append(self.truss_list[-1].new_neighbour_rel((x/2, height)))
            self.truss_list.append(self.truss_list[-1].new_neighbour_rel((x/2, -height)))

        for i in range(0,len(relative_base_xs)*2,2):
            self.truss_list[i].add_neighbour(self.truss_list[i+2])

        for i in range(1,len(relative_base_xs)*2-1,2):
            self.truss_list[i].add_neighbour(self.truss_list[i+2])

        if self.flattop:
            self.truss_list.append(self.truss_list[-1].new_neighbour_rel((0,height)))
            self.truss_list.append(self.truss_list[0].new_neighbour_rel((0,height)))

            self.truss_list[-2].add_neighbour(self.truss_list[-4])
            self.truss_list[-1].add_neighbour(self.truss_list[1])

    def update_lengths(self):
        for node in self.truss_list:
            for neighbour in node.neighbours:
                key = tuple(sorted([self.truss_list.index(node), self.truss_list.index(neighbour)]))

                self.lengths[key] = find_hyp(node.calc_rel_pos(neighbour).x_dir, node.calc_rel_pos(neighbour).y_dir)
                

    def update_internal_forces(self):
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
    
    def calc_deformed_lengths(self):
        self.update_lengths()
        for key in self.internal_forces:
            self.deformed_lengths[key] = ((self.internal_forces[key] / (self.width * self.thickness)) / self.young_mod) * self.lengths[key]
    
    def add_load(self, loc, load):
        key = None

        for node in self.truss_list:
            for neighbour in node.neighbours:
                #ifs only run when node is on the "correct" side of neighbour, but this is guaranteed to happen since for any nodes a and b, a -> b and b -> a are both checked
                if node.calc_rel_pos(neighbour).y_dir == 0 and node.loc[1] - loc[1] == 0 and loc[0] <= node.loc[0] and loc[0] >= neighbour.loc[0]:
                    key = tuple(sorted([self.truss_list.index(node), self.truss_list.index(neighbour)]))
                elif node.calc_rel_pos(neighbour).y_dir != 0 and node.loc[1] - loc[1] != 0 and loc[0]:
                    if node.calc_rel_pos(neighbour).x_dir / node.calc_rel_pos(neighbour).y_dir == (loc[0] - node.loc[0]) / (loc[1] - node.loc[1]) and loc[0] <= node.loc[0] and loc[0] >= neighbour.loc[0]:
                        key = tuple(sorted([self.truss_list.index(node), self.truss_list.index(neighbour)]))
        
        if key == None:
            return None
        else:
            node_1 = self.truss_list[key[0]]
            node_2 = self.truss_list[key[1]]
            dist_to_node_1 = find_hyp(node_1.loc[0] - loc[0], node_1.loc[1] - loc[1])
            dist_to_node_2 = find_hyp(node_2.loc[0] - loc[0], node_2.loc[1] - loc[1])

            node_1_force = Vector(load.x_dir, load.y_dir, load.get_mag() * dist_to_node_2 / (dist_to_node_1 + dist_to_node_2))
            node_2_force = Vector(load.x_dir, load.y_dir, load.get_mag() * dist_to_node_1 / (dist_to_node_1 + dist_to_node_2))

            node_1.add_external_force(node_1_force)
            node_2.add_external_force(node_2_force)

            if node_1 in self.external_forces:
                self.external_forces[node_1].add_vector(node_1_force)
            else:
                self.external_forces[node_1] = node_1_force
            
            if node_2 in self.external_forces:
                self.external_forces[node_2].add_vector(node_2_force)
            else:
                self.external_forces[node_2] = node_2_force
    
    def set_rxn_locs(self, a_rxn_loc, b_rxn_loc):
        '''
        a_rxn_loc: Pivot reaction point A_x, A_y
        b_rxn_loc: Roller reaction point B_y
        '''

        self.a_rxn_loc = a_rxn_loc
        self.b_rxn_loc = b_rxn_loc

    def update_rxn_forces(self):
        f_x_net = 0
        f_y_net = 0
        m_a_net = 0

        for node in self.external_forces:
            self.external_forces[node].set_dir_comps_by_current_mag()
            f_x_net += self.external_forces[node].x_dir
            f_y_net += self.external_forces[node].y_dir

            m_a_net += (node.loc[0] - self.a_rxn_loc[0]) * self.external_forces[node].y_dir
            m_a_net += (node.loc[1] - self.a_rxn_loc[1]) * -self.external_forces[node].x_dir
        
        a_rxn = Vector(-f_x_net, -f_y_net + (m_a_net / (self.b_rxn_loc[0] - self.a_rxn_loc[0])))
        b_rxn = Vector(0, -(m_a_net / (self.b_rxn_loc[0] - self.a_rxn_loc[0])))
        
        a_rxn.set_mag_by_current_dir_comps()
        b_rxn.set_mag_by_current_dir_comps()

        self.add_load(self.a_rxn_loc, a_rxn)
        self.add_load(self.b_rxn_loc, b_rxn)

    def calc_displacement_at_node(self, point_loc, force_direction):
        self.calc_deformed_lengths()

        virtual_truss = Truss(self.flattop, self.young_mod, self.poisson, self.tension_ult, self.comp_ult, truss_cs_properties)

        for node in self.truss_list:
            virtual_truss.truss_list.append(TrussNode((node.loc[0], node.loc[1])))
        
        for i in range(len(self.truss_list)):
            for neighbour in self.truss_list[i].neighbours:
                virtual_truss.truss_list[i].add_neighbour(virtual_truss.truss_list[self.truss_list.index(neighbour)])
        
        force_direction.set_mag(1000)

        virtual_truss.add_load(point_loc, force_direction)
        virtual_truss.set_rxn_locs(self.a_rxn_loc, self.b_rxn_loc)
        virtual_truss.update_rxn_forces()
        virtual_truss.update_internal_forces()

        work = 0

        for key in self.deformed_lengths:
            work += self.deformed_lengths[key] * virtual_truss.internal_forces[key]
        
        return work / 1000

    def force_tension_failure_member(self):
        return self.tension_ult * self.width * self.thickness
    
    def force_compression_failure_member(self):
        return -self.comp_ult * self.width * self.thickness
    
    def force_buckling_member(self):
        stress_crit = (4*math.pi**2*self.young_mod) / (12*(1-self.poisson**2))*(self.thickness/self.width)**2
        failure_force = -stress_crit * self.width * self.thickness

        return failure_force
    
    def update_failures(self):
        self.failures["force_tension_failure_member"] = {}
        self.failures["force_compression_failure_member"] = {}
        self.failures["force_buckling_member"] = {}

        for key in self.internal_forces:
            self.failures["force_tension_failure_member"][key] = self.force_tension_failure_member()
            self.failures["force_compression_failure_member"][key] = self.force_compression_failure_member()
            self.failures["force_buckling_member"][key] = self.force_buckling_member()
    
    def draw_internal_forces(self):
        self.update_failures()

        x_values = []

        for i in range(len(self.internal_forces)):
            x_values.append(i)

        force_tension_failure_list = []
        force_compression_failure_list = []
        force_buckling_member_list = []
        internal_forces_list = []

        labels = []

        for key in self.internal_forces:
            force_tension_failure_list.append(self.failures["force_tension_failure_member"][key])
            force_compression_failure_list.append(self.failures["force_compression_failure_member"][key])
            force_buckling_member_list.append(self.failures["force_buckling_member"][key])
            internal_forces_list.append(self.internal_forces[key])

            labels.append(key)
        
        plt.figure()
        plt.bar(x_values, force_tension_failure_list, color="tab:blue", zorder=1)
        plt.bar(x_values, force_compression_failure_list, color="tab:red", zorder=2)
        plt.bar(x_values, force_buckling_member_list, color="tab:orange", zorder=3)
        plt.bar(x_values, self.internal_forces.values(), color="tab:green", zorder=4)

        plt.xticks(rotation=90)
        plt.gca().set_xticks(x_values)
        plt.gca().set_xticklabels(labels)
        
    def draw_truss(self, colour=False):
        max_mag = 0
        for key in self.internal_forces:
            max_mag = abs(self.internal_forces[key]) if abs(self.internal_forces[key]) > max_mag else max_mag
        
        SCALE = 0.5
        COLOUR_MULT = 255 / max_mag
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

class BeamSlice():
    def __init__(self, dx, last_slice, young_mod, poisson, tension_ult, comp_ult, shear_ult, glue_shear_ult, cs_properties):
        self.dx = dx

        self.young_mod = young_mod
        self.poisson = poisson

        self.top_thickness = cs_properties["top_thickness"]
        self.top_width = cs_properties["top_width"]
        self.web_thickness = cs_properties["web_thickness"]
        self.web_height = cs_properties["web_height"]
        self.web_spacing = cs_properties["web_spacing"]
        self.tab_width = cs_properties["tab_width"]
        self.bottom_thickness = cs_properties["bottom_thickness"]
        self.bottom_width = cs_properties["bottom_width"]

        self.tension_ult = tension_ult
        self.comp_ult = comp_ult
        self.shear_ult = shear_ult
        self.glue_shear_ult = glue_shear_ult

        self.centroidal_axis = None
        self.second_moment_area = None

        self.external_force = Vector(0, 0, 0)

        self.shear_force = 0
        self.b_moment = 0
        self.curvature = 0

        self.q_from_centroid = None
        self.q_from_glue_top = None
        self.q_from_glue_bot = None

        self.flexural_stress_top = None
        self.flexural_stress_bot = None
        self.shear_stress_max = None
        self.shear_stress_glue_top = None
        self.shear_stress_glue_bot = None

        self.last_slice = last_slice

        if self.last_slice != None:
            self.propagate_neighbour_properties()
    
    def propagate_neighbour_properties(self):
        if self.last_slice != None:
            self.shear_force = self.last_slice.shear_force
            self.b_moment = self.last_slice.b_moment
            self.curvature = self.last_slice.curvature
        else:
            self.shear_force = 0
            self.b_moment = 0
            self.curvature = 0
    
    def set_external_force(self, force):
        self.external_force = force

    def add_external_force(self, force):
        self.external_force.add_vector(force)
    
    def calc_centroidal_axis_and_second_moment_of_area(self):
        A1 = self.web_thickness*self.web_height
        y1 = self.web_height/2 + self.bottom_thickness
        A2 = self.web_thickness*self.tab_width
        y2 = self.web_height - self.web_thickness/2 + self.bottom_thickness
        A3 = self.top_width*self.top_thickness
        y3 = self.web_height + self.top_thickness/2 + self.bottom_thickness
        A4 = A2
        y4 = y2
        A5 = A1
        y5 = y1
        A6 = self.bottom_thickness*self.bottom_width
        y6 = self.bottom_thickness/2
        A7 = 0
        y7 = 0
        A8 = 0
        y8 = 0

        #COMMENT OUT TO SOLVE DESIGN 0
        if A6 != 0:
            A7 = A2
            y7 = self.web_thickness/2 + self.bottom_thickness
            A8 = A2
            y8 = y7

        self.centroidal_axis = (A1*y1 + A2*y2 + A3*y3 + A4*y4 + A5*y5 + A6*y6 + A7*y7 + A8*y8)/(A1 + A2 + A3 + A4 + A5 + A6 + A7 + A8)

        I1 = self.web_thickness*self.web_height**3/12
        I2 = self.tab_width*self.web_thickness**3/12
        I3 = self.top_width*self.top_thickness**3/12
        I4 = I2
        I5 = I1
        I6 = self.bottom_width*self.bottom_thickness**3/12
        I7 = 0
        I8 = 0

        #COMMENT OUT TO SOLVE DESIGN 0
        if A6 != 0:
            I7 = self.tab_width*self.web_thickness**3/12
            I8 = self.tab_width*self.web_thickness**3/12

        self.second_moment_area = I1 + I2 + I3 + I4 + I5 + I6 + I7 + I8 + A1*(self.centroidal_axis-y1)**2 + A2*(self.centroidal_axis-y2)**2 + A3*(self.centroidal_axis-y3)**2 + A4*(self.centroidal_axis-y4)**2 + A5*(self.centroidal_axis-y5)**2 + A6*(self.centroidal_axis-y6)**2 + + A7*(self.centroidal_axis-y7)**2 + + A8*(self.centroidal_axis-y8)**2

    def calc_Qs(self):
        self.calc_centroidal_axis_and_second_moment_of_area()
        A1 = self.top_width*self.top_thickness
        d1 = self.bottom_thickness + self.web_height + self.top_thickness/2 - self.centroidal_axis
        A2 = (self.bottom_thickness + self.web_height - self.centroidal_axis)*self.web_thickness
        d2 = (self.bottom_thickness + self.web_height - self.centroidal_axis) / 2
        A3 = self.web_thickness*self.tab_width
        d3 = (self.bottom_thickness + self.web_height - self.web_thickness/2 - self.centroidal_axis)
        A4 = A3
        d4 = d3
        A5 = A2
        d5 = d2

        self.q_from_centroid = A1*d1 + A2*d2 + A3*d3 + A4*d4 + A5*d5

        self.q_from_glue_top = A1*d1

        A4 = self.bottom_width*self.bottom_thickness
        d4 = self.bottom_thickness/2 - self.centroidal_axis

        self.q_from_glue_bot = A4*d4
    
    def update_internal_properties(self):
        self.propagate_neighbour_properties()
        self.calc_centroidal_axis_and_second_moment_of_area()

        self.external_force.set_dir_comps_by_current_mag()
        self.shear_force += self.external_force.get_y_mag()
        
        self.b_moment += self.shear_force * self.dx
        self.curvature = self.b_moment / (self.young_mod * self.second_moment_area)

        self.flexural_stress_top = -self.b_moment * (self.bottom_thickness + self.web_height + self.top_thickness - self.centroidal_axis) / self.second_moment_area
        self.flexural_stress_bot = -self.b_moment * (-self.centroidal_axis) / self.second_moment_area

        self.calc_Qs()

        self.shear_stress_max = (self.shear_force * self.q_from_centroid) / (self.second_moment_area * self.web_thickness * 2)
        self.shear_stress_glue_top = (self.shear_force * self.q_from_glue_top) / (self.second_moment_area * self.web_thickness * 2)
        self.shear_stress_glue_bot = (self.shear_force * self.q_from_glue_bot) / (self.second_moment_area * self.web_thickness * 2)
    
    def moment_tension_failure_wall(self):
        if self.b_moment < 0:
            failure_moment = (self.tension_ult * self.second_moment_area) / self.centroidal_axis
        else:
            failure_moment = (self.tension_ult * self.second_moment_area) / (self.bottom_thickness + self.web_height + self.top_thickness - self.centroidal_axis)
        
        return failure_moment

    def moment_compression_failure_wall(self):
        if self.b_moment < 0:
            failure_moment = (self.comp_ult * self.second_moment_area) / (self.bottom_thickness + self.web_height + self.top_thickness - self.centroidal_axis)
        else:
            failure_moment = (self.comp_ult * self.second_moment_area) / self.centroidal_axis
        
        return failure_moment
    
    def force_shear_failure_wall(self):
        self.calc_Qs()
        failure_force = self.shear_ult * 2 * self.web_thickness * self.second_moment_area / self.q_from_centroid

        return failure_force

    def force_shear_failure_glue_top(self):
        self.calc_Qs()
        failure_force = self.glue_shear_ult * 2 * self.web_thickness * self.second_moment_area / self.q_from_glue_top
        
        return failure_force
    
    def force_shear_failure_glue_bot(self):
        self.calc_Qs()
        failure_force = self.glue_shear_ult * 2 * self.web_thickness * self.second_moment_area / self.q_from_glue_bot
        
        return failure_force
    
    def moment_buckling_compressive_top(self):
        stress_crit = (4*math.pi**2*self.young_mod) / (12*(1-self.poisson**2)) * (self.top_thickness/self.web_spacing)**2
        failure_moment = stress_crit * self.second_moment_area / ((self.bottom_thickness + self.web_height + self.top_thickness/2) - self.centroidal_axis)
        
        return failure_moment

    def moment_buckling_compressive_bot(self):
        stress_crit = (4*math.pi**2*self.young_mod) / (12*(1-self.poisson**2)) * (self.bottom_thickness/self.bottom_width)**2
        failure_moment = stress_crit * self.second_moment_area / (self.bottom_thickness/2 - self.centroidal_axis)
        
        return failure_moment

    def moment_buckling_flanges(self): 
        stress_crit = (0.425*math.pi**2*self.young_mod) / (12*(1-self.poisson**2)) * (self.top_thickness/((self.top_width - self.bottom_width)/2))**2
        failure_moment = stress_crit * self.second_moment_area / ((self.bottom_thickness + self.web_height + self.top_thickness/2) - self.centroidal_axis)
        
        return failure_moment

    def moment_buckling_webs_flexural(self): #web_thickness -> top_thickness?
        stress_crit = (6*math.pi**2*self.young_mod) / (12*(1-self.poisson**2)) * (self.web_thickness/(self.bottom_thickness + self.web_height - self.centroidal_axis))**2
        failure_moment = stress_crit * self.second_moment_area / (self.centroidal_axis - self.bottom_thickness)
        
        return failure_moment

class Beam():
    def __init__(self, young_mod, poisson, tension_ult, comp_ult, shear_ult, glue_shear_ult, segments, diaphragms, resolution, a_rxn_loc=None, b_rxn_loc=None):
        '''
        young_mod: Young's Modulus (MPa)
        poisson: Poisson's Ratio
        tension_ult: Tensile Strength (MPa)
        comp_ult: Compressive Strength (MPa)
        shear_ult: Shear Strength (MPa)
        glue_shear_ult: Glue Shear Strength (MPa)

        segments: Cross-Sectional Segment Properties
        diaphragms: Diaphragm Locations

        resolution: Number of Beam Slices (higher -> more accuracy)

        a_rxn_loc: Pivot reaction point A_x, A_y
        b_rxn_loc: Roller reaction point B_y
        '''

        self.beam_list = []

        self.young_mod = young_mod
        self.poisson = poisson

        self.tension_ult = tension_ult
        self.comp_ult = comp_ult
        self.shear_ult = shear_ult
        self.glue_shear_ult = glue_shear_ult

        self.segments = sorted(segments)
        self.diaphragms = diaphragms
        self.length = segments[-1][0][1]

        self.external_forces = {}

        self.a_rxn_loc = a_rxn_loc
        self.b_rxn_loc = b_rxn_loc

        self.sfd = []
        self.bmd = []
        self.curv = []

        self.failures = {}

        self.resolution = resolution
        self.dx = self.length / self.resolution
    
    def gen_beam(self):
        for i in range(self.resolution):
            cs_properties = None
            for segment in self.segments:
                if i*self.dx >= segment[0][0] and i*self.dx < segment[0][1]:
                    cs_properties = segment[1]
            
            if cs_properties == None:
                print("ERROR: location", i*self.dx, "has no defined segment")
                return

            if len(self.beam_list) == 0:
                self.beam_list.append(BeamSlice(self.dx, None, self.young_mod, self.poisson, self.tension_ult, self.comp_ult, self.shear_ult, self.glue_shear_ult, cs_properties))
            else:
                self.beam_list.append(BeamSlice(self.dx, self.beam_list[-1], self.young_mod, self.poisson, self.tension_ult, self.comp_ult, self.shear_ult, self.glue_shear_ult, cs_properties))
            
            self.beam_list[-1].update_internal_properties()
    
    def update_external_forces(self):
        for i in range(self.resolution):
            self.beam_list[i].set_external_force(Vector(0,0,0))
            for loc in self.external_forces:
                if loc[0] >= i*self.dx and loc[0] < (i+1)*self.dx:
                    self.beam_list[i].add_external_force(self.external_forces[loc])
            
            self.beam_list[i].update_internal_properties()

    def add_load(self, loc, load):
        if loc in self.external_forces:
            self.external_forces[loc].add_vector(load)
        else:
            self.external_forces[loc] = load

        self.update_external_forces()
    
    def set_rxn_locs(self, a_rxn_loc, b_rxn_loc):
        '''
        a_rxn_loc: Pivot reaction point A_x, A_y
        b_rxn_loc: Roller reaction point B_y
        '''

        self.a_rxn_loc = a_rxn_loc
        self.b_rxn_loc = b_rxn_loc
    
    def update_rxn_forces(self):
        f_x_net = 0
        f_y_net = 0
        m_a_net = 0

        for loc in self.external_forces:
            self.external_forces[loc].set_dir_comps_by_current_mag()
            f_x_net += self.external_forces[loc].x_dir
            f_y_net += self.external_forces[loc].y_dir

            m_a_net += (loc[0] - self.a_rxn_loc[0]) * self.external_forces[loc].y_dir
            m_a_net += (loc[1] - self.a_rxn_loc[1]) * -self.external_forces[loc].x_dir
        
        a_rxn = Vector(-f_x_net, -f_y_net + (m_a_net / (self.b_rxn_loc[0] - self.a_rxn_loc[0])))
        b_rxn = Vector(0, -(m_a_net / (self.b_rxn_loc[0] - self.a_rxn_loc[0])))
        
        a_rxn.set_mag_by_current_dir_comps()
        b_rxn.set_mag_by_current_dir_comps()

        self.add_load(self.a_rxn_loc, a_rxn)
        self.add_load(self.b_rxn_loc, b_rxn)
    
    def forces_shear_buckling_webs(self):
        web_shear_failure_forces = {}

        for i in range(len(self.diaphragms)-1):
            span = self.diaphragms[i], self.diaphragms[i+1]
            width = span[1] - span[0]

            initial_slice = None
            for i in range(self.resolution):
                if span[0] >= i*self.dx and span[0] < (i+1)*self.dx:
                    initial_slice = self.beam_list[i]
            
            if initial_slice == None:
                print("ERROR: no defined cross-sectional properties for diaphragm creation at", span[0])

            initial_slice.calc_Qs()
            stress_crit = (5*math.pi**2*self.young_mod)/(12*(1-self.poisson**2))*((initial_slice.web_thickness/(width*2/3))**2 +((initial_slice.web_thickness)/(initial_slice.web_height))**2)
            
            failure_force = stress_crit * 2 * initial_slice.web_thickness * initial_slice.second_moment_area / initial_slice.q_from_centroid

            web_shear_failure_forces[span] = failure_force
        
        return web_shear_failure_forces
    
    def update_internal_properties(self):
        self.sfd = []
        self.bmd = []
        self.curv = []

        for i in range(self.resolution):
            self.sfd.append(self.beam_list[i].shear_force)
            self.bmd.append(self.beam_list[i].b_moment)
            self.curv.append(self.beam_list[i].curvature)

    def update_failures(self):
        self.failures["moment_tension_failure_wall"] = []
        self.failures["moment_compression_failure_wall"] = []
        self.failures["force_shear_failure_wall"] = []
        self.failures["force_shear_failure_glue_top"] = []
        self.failures["force_shear_failure_glue_bot"] = []
        self.failures["moment_buckling_compressive_top"] = []
        self.failures["moment_buckling_compressive_bot"] = []
        self.failures["moment_buckling_flanges"] = []
        self.failures["moment_buckling_webs_flexural"] = []
        self.failures["forces_buckling_webs_shear"] = []

        for i in range(self.resolution):
            beam_slice = self.beam_list[i]

            self.failures["moment_tension_failure_wall"].append(beam_slice.moment_tension_failure_wall())
            self.failures["moment_compression_failure_wall"].append(beam_slice.moment_compression_failure_wall())
            self.failures["force_shear_failure_wall"].append(beam_slice.force_shear_failure_wall())
            self.failures["force_shear_failure_glue_top"].append(beam_slice.force_shear_failure_glue_top())
            self.failures["moment_buckling_compressive_top"].append(beam_slice.moment_buckling_compressive_top())

            if beam_slice.bottom_thickness != 0 and beam_slice.bottom_width != 0:
                self.failures["force_shear_failure_glue_bot"].append(beam_slice.force_shear_failure_glue_bot())
                self.failures["moment_buckling_compressive_bot"].append(beam_slice.moment_buckling_compressive_bot())
            else:
                self.failures["force_shear_failure_glue_bot"].append(np.nan)
                self.failures["moment_buckling_compressive_bot"].append(np.nan)

            self.failures["moment_buckling_flanges"].append(beam_slice.moment_buckling_flanges())
            self.failures["moment_buckling_webs_flexural"].append(beam_slice.moment_buckling_webs_flexural())

            web_shear_failure_forces = self.forces_shear_buckling_webs()
            web_shear_failure_force = np.nan
            for span in web_shear_failure_forces:
                if i*self.dx >= span[0] and i*self.dx < span[1]:
                    web_shear_failure_force = web_shear_failure_forces[span]
                
            self.failures["forces_buckling_webs_shear"].append(web_shear_failure_force)
    
    def draw_internal_properties(self):
        self.update_internal_properties()
        self.update_failures()

        x_values = []

        for i in range(self.resolution):
            x_values.append(i*self.dx)
        
        modes = ["force_shear_failure_wall", "force_shear_failure_glue_top", "force_shear_failure_glue_bot", "forces_buckling_webs_shear"]
        plt.figure()
        for i in range(4):
            plt.subplot(2,2,i+1)
            plt.plot(x_values, self.sfd)
            plt.plot(x_values, self.failures[modes[i]])
            plt.title("SFD vs. " +modes[i])

        modes = ["moment_tension_failure_wall", "moment_compression_failure_wall", "moment_buckling_compressive_top", "moment_buckling_compressive_bot", "moment_buckling_flanges", "moment_buckling_webs_flexural"]
        plt.figure()
        for i in range(6):
            plt.subplot(3,2,i+1)
            plt.gca().invert_yaxis()
            plt.plot(x_values, self.bmd)
            plt.plot(x_values, self.failures[modes[i]])
            plt.title("BMD vs. " +modes[i])

        plt.figure()
        plt.gca().invert_yaxis()
        plt.plot(x_values, self.curv)
        plt.title("Curvature")

class Bridge():
    def __init__(self, truss_flattop, young_mod, poisson, tension_ult, comp_ult, shear_ult, glue_shear_ult, truss_cs_properties, beam_segments, beam_diaphragms, resolution=1000):
        self.truss = Truss(truss_flattop, young_mod, poisson, tension_ult, comp_ult, truss_cs_properties)
        self.beam = Beam(young_mod, poisson, tension_ult, comp_ult, shear_ult, glue_shear_ult, beam_segments, beam_diaphragms, resolution)
    
    def gen_bridge_periodic_truss(self, triangle_count, diagonal):
        self.beam.gen_beam()
        self.truss.gen_periodic_warren_truss(triangle_count, diagonal)

    def gen_bridge_nonperiodic_truss(self, relative_base_xs, height):
        self.beam.gen_beam()
        self.truss.gen_nonperiodic_warren_truss(relative_base_xs, height)

    def set_rxn_locs(self, a_rxn_loc, b_rxn_loc):
        self.truss.set_rxn_locs(a_rxn_loc, b_rxn_loc)
        self.beam.set_rxn_locs(a_rxn_loc, b_rxn_loc)
    
    def update_bridge(self):
        self.truss.update_rxn_forces()
        self.beam.update_rxn_forces()

        self.truss.update_internal_forces()
        self.truss.update_failures()
        self.beam.update_internal_properties()
        self.beam.update_failures()
    
    def add_load(self, loc, load):
        self.truss.add_load(loc, Vector(load.x_dir, load.y_dir, load.get_mag() * 0.5))
        self.beam.add_load(loc, Vector(load.x_dir, load.y_dir, load.get_mag() * 0.5))

'''--edit--'''
bridge_height = 140

truss_cs_properties =  {"width": 80, "thickness": 1.27}

beam_segments =    [((0,1280),      {"top_thickness": 1.27,
                                    "top_width": 100,
                                    "web_thickness": 1.27,
                                    "web_height": 140-1.27*2,
                                    "web_spacing": 75-1.27*2,
                                    "tab_width": 10,
                                    "bottom_thickness": 1.27,
                                    "bottom_width": 75})
                    ]

truss_relative_base_xs = [30, 245, 245, 90, 210, 210, 90, 130, 30]

'''--do not edit--'''
beam_diaphragms = []

x=0
for i in truss_relative_base_xs:
    beam_diaphragms.append(x + i/4)
    beam_diaphragms.append(x + 3*i/4)
    x += i

bridge = Bridge(True, 4000, 0.2, 30, 6, 4, 2, truss_cs_properties, beam_segments, beam_diaphragms, 1000)
bridge.gen_bridge_nonperiodic_truss(truss_relative_base_xs, bridge_height)

bridge.set_rxn_locs((15,0), (1075,0))
bridge.add_load((565,0), Vector(0,-1,1000))
bridge.add_load((1265,0), Vector(0,-1,1000))

bridge.update_bridge()

bridge.beam.draw_internal_properties()
bridge.truss.draw_internal_forces()

print(bridge.truss.calc_displacement_at_node((640), Vector(0, -1, 1)))

plt.show()

bridge.truss.draw_truss(True)