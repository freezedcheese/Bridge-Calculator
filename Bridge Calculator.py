import math
import turtle
import matplotlib.pyplot as plt
import numpy as np

def find_hyp(x, y):
    '''
    Helper function to find hypoteneuse between points
    '''

    return math.sqrt(x**2 + y**2)

class Vector():
    def __init__(self, x_dir, y_dir, mag=1):
        '''
        Vector helper class to simplify repeated vector calculations
        Directions are set using x and y components, which are NOT related to the vector's magnitude
        The vector's magnitude is stored separately
        This system simplifies some calculations

        x_dir: x component of direction vector
        y_dir: y component of direction vector
        mag: Magnitude of vector
        '''

        self.x_dir = x_dir
        self.y_dir = y_dir
        self.mag = mag

    def get_opposite(self):
        '''
        Return vector with opposite direction and equal magnitude
        '''
        
        return Vector(-self.x_dir, -self.y_dir, self.mag)
    
    def get_mag(self):
        '''
        Return vector magnitude
        '''

        return self.mag
    
    def get_x_mag(self):
        '''
        Return vector x magnitude
        '''

        if find_hyp(self.x_dir, self.y_dir) == 0:
            return 0
        else:
            return self.mag * (self.x_dir / find_hyp(self.x_dir, self.y_dir))
    
    def get_y_mag(self):
        '''
        Return y vector magnitude
        '''

        if find_hyp(self.x_dir, self.y_dir) == 0:
            return 0
        else:
            return self.mag * (self.y_dir / find_hyp(self.x_dir, self.y_dir))
    
    def set_mag(self, mag):
        '''
        Set vector magnitude
        '''

        self.mag = mag
    
    def set_mag_by_comp(self, x_mag_comp, y_mag_comp):
        '''
        Set vector magnitude by setting magnitude components, while maintaining direction vector
        '''

        if self.x_dir != 0:
            self.mag = x_mag_comp * (find_hyp(self.x_dir, self.y_dir) / self.x_dir)
        elif self.y_dir != 0:
            self.mag = y_mag_comp * (find_hyp(self.x_dir, self.y_dir) / self.y_dir)
        else:
            self.mag = 0
    
    def set_mag_by_current_dir_comps(self):
        '''
        Set vector magnitude to magnitude of current direction vector
        '''

        self.mag = find_hyp(self.x_dir, self.y_dir)
    
    def set_dir_comps_by_current_mag(self):
        '''
        Set direction vector component magnitudes to current component magnitudes
        '''

        self.x_dir = self.get_x_mag()
        self.y_dir = self.get_y_mag()
    
    def add_vector(self, vector):
        '''
        Add another vector to this one by adding magnitude components
        This also modifies the current direction vector - the resulting direction vector has a magnitude equal to that of the resulting vector
        '''

        self.set_dir_comps_by_current_mag()
        vector.set_dir_comps_by_current_mag()

        self.x_dir += vector.x_dir
        self.y_dir += vector.y_dir

        self.set_mag_by_current_dir_comps()


class TrussNode():
    def __init__(self, loc=(0,0)):
        '''
        TrussNode class (each node represents a connection point between multiple members)
        Used as members of a truss shaped data structure in the Truss class
        Used by the Truss class to calculate internal forces

        loc: (x,y) coordinates of node [tuple]
        '''

        self.loc = loc
        self.neighbours = set()
        self.external_forces = []
        self.internal_forces = {}

    def new_neighbour_rel(self, rel):
        '''
        Create new neighbouring node positioned relative to this node, linked as a neighbour to this node
        '''

        neighbour_loc = (self.loc[0] + rel[0], self.loc[1] + rel[1])
        neighbour = TrussNode(neighbour_loc)
        self.neighbours.add(neighbour)
        neighbour.neighbours.add(self)

        return neighbour
    
    def add_neighbour(self, neighbour):
        '''
        Add preexisting node as a neighbour to this node
        
        neighbour: Neighbour node reference [TrussNode]
        '''

        self.neighbours.add(neighbour)
        neighbour.neighbours.add(self)
    
    def add_external_force(self, force):
        '''
        Add truss external force vector
        '''

        self.external_forces.append(force)
    
    def calc_rel_pos(self, neighbour):
        '''
        Calculate relative position between this node and another

        neighbour: Neighbour node reference [TrussNode]
        '''

        return Vector(neighbour.loc[0] - self.loc[0], neighbour.loc[1] - self.loc[1])

    def calc_internal_forces(self, debug_print=False):
        '''
        Calculate truss internal forces required between this node and all neighbouring nodes to satisfy equilibrium
        '''

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
                if self.internal_forces[unknown_internal_forces[1]].get_mag() == 0:
                    self.internal_forces[unknown_internal_forces[0]].set_mag_by_comp(-f_x_net, -f_y_net)
                else:
                    self.internal_forces[unknown_internal_forces[0]].set_mag((f_a.get_mag()/f_a.get_x_mag()) * (-f_x_net - self.internal_forces[unknown_internal_forces[1]].get_mag() * (f_b.get_x_mag()/f_b.get_mag())))
        
        elif len(unknown_internal_forces) == 1:
            self.internal_forces[unknown_internal_forces[0]].set_mag_by_comp(-f_x_net, -f_y_net)

        if debug_print:
            for thing in self.neighbours:
                print(self.internal_forces[thing].get_mag())

class Truss():
    def __init__(self, flattop, young_mod, poisson, tension_ult, comp_ult, glue_shear_ult, cs_properties, modified_cs_members=None, modified_tab_length_nodes=None, a_rxn_loc=None, b_rxn_loc=None):
        '''
        Truss data structure consisting of TrussNodes
        Contains TrussNode elements representing truss joints
        Can solve for all internal forces within any shape of truss using method of joints - including the bat shape from Friday afternoon's quiz 8! :)
        This is extremely useful when solving asymmetrical, nonperiodic trusses (like the final design we settled on)

        flattop: True if bridge has a flat top (and bottom) [bool]

        young_mod: Young's modulus (MPa) [float]
        poisson: Poisson's ratio [float]

        tension_ult: Tensile strength (MPa) [float]
        comp_ult: Compressive strength (MPa) [float]
        glue_shear_ult: Glue shear strength (MPa) [float]

        cs_properties: Cross-sectional properties [dict]

        modified_cs_members: Modified cross sectional properties [dict]
        modifies_tab_length_nodes: Modified glue tab lengths [dict]

        a_rxn_loc: (x,y) pivot reaction point A_x, A_y [tuple]
        b_rxn_loc: (x,y) roller reaction point B_y [tuple]
        '''

        self.flattop = flattop

        self.tension_ult = tension_ult
        self.comp_ult = comp_ult
        self.glue_shear_ult = glue_shear_ult

        self.cs_properties = cs_properties
        self.width = cs_properties["width"]
        self.thickness = cs_properties["thickness"]
        self.tab_width = cs_properties["tab_width"]
        self.tab_length = cs_properties["tab_length"]

        self.modified_cs_members = modified_cs_members
        self.modified_tab_length_nodes = modified_tab_length_nodes

        self.young_mod = young_mod
        self.poisson = poisson

        self.truss_list = []
        self.internal_forces = {}
        self.external_forces = {}
        self.node_forces = []

        self.lengths = {}
        self.deformed_lengths = {}

        self.failures = {}
        self.min_FOS = {}

        self.a_rxn_loc = a_rxn_loc
        self.b_rxn_loc = b_rxn_loc

    def gen_periodic_warren_truss(self, triangle_count, diagonal, flipped=False):
        '''
        Generate Truss data structure consisting of TrussNodes that represent a periodic and symmetrical Warren truss

        triangle_count: Number of triangles in Warren truss [int]
        diagonal: (x, y) distance between first 2 nodes [tuple]
        flipped: True if bridge is upside-down [bool]
        '''

        self.truss_list = []
        self.truss_list.append(TrussNode((0, 0) if not flipped else (0, diagonal[1])))

        diagonal = diagonal if not flipped else (diagonal[0], -diagonal[1])

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
    
    def gen_nonperiodic_warren_truss(self, relative_base_xs, height, flipped=False):
        '''
        Generate Truss data structure consisting of TrussNodes that represent a nonperiodic and asymmetrical Warren truss

        relative_base_xs: Relative x distances between the nodes at the base of the truss [list]
        height: height of bridge [float]
        flipped: True if bridge is upside-down [bool]
        '''

        self.truss_list = []
        self.truss_list.append(TrussNode((0, 0) if not flipped else (0, height)))

        height = height if not flipped else -height

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
        '''
        Update lengths of Truss members
        '''

        for node in self.truss_list:
            for neighbour in node.neighbours:
                key = tuple(sorted([self.truss_list.index(node), self.truss_list.index(neighbour)]))

                self.lengths[key] = find_hyp(node.calc_rel_pos(neighbour).x_dir, node.calc_rel_pos(neighbour).y_dir)
                

    def update_internal_forces(self):
        '''
        Update internal member forces in Truss
        Makes calls to each TrussNode to calculate the internal forces at each joint, then consolidates them

        Forces:
        - Internal Member Forces
        - Joint Shear Forces
        '''

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

        self.node_forces = []
        for node in self.truss_list:
            f_x_net = 0
             
            for neighbour in node.internal_forces:
                if neighbour.loc[0] - node.loc[0] < 0:
                    f_x_net += node.internal_forces[neighbour].get_x_mag()

            self.node_forces.append(abs(f_x_net))

    def calc_deformed_lengths(self):
        '''
        Calculate deformed lengths of Truss members
        '''

        self.update_lengths()
        for key in self.internal_forces:
            self.deformed_lengths[key] = ((self.internal_forces[key] / (self.width * self.thickness)) / self.young_mod) * self.lengths[key]
    
    def add_load(self, loc, load):
        '''
        Add external load to Truss

        loc: (x,y) location of force [tuple]
        load: Load force [Vector]
        '''

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
        Set locations of reaction forces on Truss

        a_rxn_loc: (x,y) pivot reaction point A_x, A_y [tuple]
        b_rxn_loc: (x,y) roller reaction point B_y [tuple]
        '''

        self.a_rxn_loc = a_rxn_loc
        self.b_rxn_loc = b_rxn_loc

    def update_rxn_forces(self):
        '''
        Update reaction forces on Truss based on current acting external loads
        '''

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

    def calc_displacement(self, point_loc, force_direction):
        '''
        Calculate displacement at point of interest using method of virtual work

        point_loc: (x,y) point of interest [tuple]
        force_direction: Displacement direction vector [Vector]
        '''

        self.calc_deformed_lengths()

        virtual_truss = Truss(self.flattop, self.young_mod, self.poisson, self.tension_ult, self.comp_ult, self.glue_shear_ult, self.cs_properties)

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

    def force_tension_failure_member(self, key):
        '''
        Calculate tensile force required to cause tensile failure of member with id "key"

        key: Member key of form (node_a_index, node_b_index) [tuple]
        '''

        width = self.width if key not in self.modified_cs_members else self.modified_cs_members[key]["width"]
        thickness = self.thickness if key not in self.modified_cs_members else self.modified_cs_members[key]["thickness"]

        return self.tension_ult * width * thickness
    
    def force_compression_failure_member(self, key):
        '''
        Calculate compressive force required to cause compressive yielding failure of member with id "key"

        key: Member key of form (node_a_index, node_b_index) [tuple]
        '''

        width = self.width if key not in self.modified_cs_members else self.modified_cs_members[key]["width"]
        thickness = self.thickness if key not in self.modified_cs_members else self.modified_cs_members[key]["thickness"]

        return -self.comp_ult * width * thickness
    
    def force_buckling_member(self, key):
        '''
        Calculate compressive force required to cause global buckling of member with id "key"

        key: Member key of form (node_a_index, node_b_index) [tuple]
        '''

        width = self.width if key not in self.modified_cs_members else self.modified_cs_members[key]["width"]
        thickness = self.thickness if key not in self.modified_cs_members else self.modified_cs_members[key]["thickness"]

        stress_crit = (4*math.pi**2*self.young_mod) / (12*(1-self.poisson**2))*(thickness/width)**2
        failure_force = -stress_crit * width * thickness

        return failure_force
    
    def force_shear_failure_glue_tab(self, i):
        '''
        Calculate shear force required to cause shear failure at joint (node) with id "i" 

        i: Node id [int]
        '''

        tab_length = self.tab_length if i not in self.modified_tab_length_nodes else self.modified_tab_length_nodes[i]

        failure_force = self.glue_shear_ult * (self.tab_width * tab_length**3 / 12) * self.tab_width / (tab_length * self.tab_width * tab_length / 4)

        return failure_force
    
    def update_failures(self):
        '''
        Update forces required to cause failure by each failure mode for every single truss joint and member

        Forces:
        - Internal Member Forces
        - Joint Shear Forces

        Failure modes:
        - Truss Member Tensile Yield Failure
        - Truss Member Compressive Yield Failure
        - Truss Member Local Buckling Failure
        - Truss Joint Glue Tab Shear Failure
        '''

        self.failures["force_tension_failure_member"] = {}
        self.failures["force_compression_failure_member"] = {}
        self.failures["force_buckling_member"] = {}
        self.failures["force_shear_failure_glue_tab"] = {}

        for key in self.internal_forces:
            self.failures["force_tension_failure_member"][key] = self.force_tension_failure_member(key)
            self.failures["force_compression_failure_member"][key] = self.force_compression_failure_member(key)
            self.failures["force_buckling_member"][key] = self.force_buckling_member(key)

        for i in range(len(self.truss_list)):            
            self.failures["force_shear_failure_glue_tab"][i] = self.force_shear_failure_glue_tab(i)
    
    def update_min_FOS(self):
        '''
        Update FOS for each failure mode by iterating through the forces in every joint and member and finding the joints and members with the lowest FOS in each mode

        Forces:
        - Internal Member Forces
        - Joint Shear Forces

        Failure modes:
        - Truss Member Tensile Yield Failure
        - Truss Member Compressive Yield Failure
        - Truss Member Local Buckling Failure
        - Truss Joint Glue Tab Shear Failure
        '''

        member_modes = ["force_tension_failure_member", "force_compression_failure_member", "force_buckling_member"]

        for mode in member_modes:
            self.min_FOS[mode] = None
            for key in self.internal_forces:
                if self.internal_forces[key] != 0:
                    if self.failures[mode][key] / self.internal_forces[key] > 0:
                        if self.min_FOS[mode] == None:
                            self.min_FOS[mode] = self.failures[mode][key] / self.internal_forces[key]
                        else:
                            self.min_FOS[mode] = min(self.failures[mode][key] / self.internal_forces[key], self.min_FOS[mode])
            
        node_mode = "force_shear_failure_glue_tab"

        self.min_FOS[node_mode] = None
        for i in range(len(self.node_forces)):
            if self.node_forces[i] != 0:
                    if self.failures[node_mode][i] / self.node_forces[i] > 0:
                        if self.min_FOS[node_mode] == None:
                            self.min_FOS[node_mode] = self.failures[node_mode][i] / self.node_forces[i]
                        else:
                            self.min_FOS[node_mode] = min(self.failures[node_mode][i] / self.node_forces[i], self.min_FOS[node_mode])
    
    def update_min_envelope_FOS(self, force_envelopes):
        '''
        Update FOS for each failure mode given truss force envelopes by iterating through and finding the joints and members with the lowest FOS at any point in time

        force_envelopes: Dictionary containining force envelope data for each force [dict]

        Forces:
        - Internal Member Forces
        - Joint Shear Forces

        Failure modes:
        - Truss Member Tensile Yield Failure
        - Truss Member Compressive Yield Failure
        - Truss Member Local Buckling Failure
        - Truss Joint Glue Tab Shear Failure
        '''

        member_modes = ["force_tension_failure_member", "force_compression_failure_member", "force_buckling_member"]

        for mode in member_modes:
            self.min_FOS[mode] = None
            for key in self.internal_forces:
                for internal_forces_envelope in (force_envelopes["truss_internal_forces_positive_envelope"], force_envelopes["truss_internal_forces_negative_envelope"]):
                    if internal_forces_envelope[key] != 0:
                        if self.failures[mode][key] / internal_forces_envelope[key] > 0:
                            if self.min_FOS[mode] == None:
                                self.min_FOS[mode] = self.failures[mode][key] / internal_forces_envelope[key]
                            else:
                                self.min_FOS[mode] = min(self.failures[mode][key] / internal_forces_envelope[key], self.min_FOS[mode])
            
        node_mode = "force_shear_failure_glue_tab"

        self.min_FOS[node_mode] = None
        for i in range(len(self.node_forces)):
            for node_forces_envelope in (force_envelopes["truss_node_forces_positive_envelope"], force_envelopes["truss_node_forces_negative_envelope"]):
                if node_forces_envelope[i] != 0:
                    if self.failures[node_mode][i] / node_forces_envelope[i] > 0:
                        if self.min_FOS[node_mode] == None:
                            self.min_FOS[node_mode] = self.failures[node_mode][i] / node_forces_envelope[i]
                        else:
                            self.min_FOS[node_mode] = min(self.failures[node_mode][i] / node_forces_envelope[i], self.min_FOS[node_mode])

    def draw_internal_forces(self, override_internal_forces_positive_dict=None, override_internal_forces_negative_dict=None, override_node_forces_positive_list=None, override_node_forces_negative_list=None):
        '''
        Plot all joint and member internal forces along with corresponding FOS for each failure mode

        Forces:
        - Internal Member Forces
        - Joint Shear Forces

        Failure modes:
        - Truss Member Tensile Yield Failure
        - Truss Member Compressive Yield Failure
        - Truss Member Local Buckling Failure
        - Truss Joint Glue Tab Shear Failure
        '''

        self.update_failures()

        member_x_values = []

        for i in range(len(self.internal_forces)):
            member_x_values.append(i)

        force_tension_failure_list = []
        force_compression_failure_list = []
        force_buckling_member_list = []
        
        internal_forces_list = []

        if override_internal_forces_positive_dict != None and override_internal_forces_negative_dict != None:
            override_internal_forces_positive_list = []
            override_internal_forces_negative_list = []

        member_labels = []

        for key in self.internal_forces:
            force_tension_failure_list.append(self.failures["force_tension_failure_member"][key])
            force_compression_failure_list.append(self.failures["force_compression_failure_member"][key])
            force_buckling_member_list.append(self.failures["force_buckling_member"][key])
            internal_forces_list.append(self.internal_forces[key])
            if override_internal_forces_positive_dict != None and override_internal_forces_negative_dict != None:
                override_internal_forces_positive_list.append(override_internal_forces_positive_dict[key])
                override_internal_forces_negative_list.append(override_internal_forces_negative_dict[key])

            member_labels.append(key)
        
        plt.figure()
        plt.bar(member_x_values, force_tension_failure_list, color="tab:blue", zorder=1)
        plt.bar(member_x_values, force_compression_failure_list, color="tab:red", zorder=2)
        plt.bar(member_x_values, force_buckling_member_list, color="tab:orange", zorder=3)
        if override_internal_forces_positive_dict != None and override_internal_forces_negative_dict != None:
            plt.bar(member_x_values, override_internal_forces_positive_list, color="tab:green", zorder=2)
            plt.bar(member_x_values, override_internal_forces_negative_list, color="tab:green", zorder=2)
        else:
            plt.bar(member_x_values, internal_forces_list, color="tab:green", zorder=2)
        plt.title("Member Loads")
        plt.ylabel("Force (N)")

        plt.xticks(rotation=90)
        plt.gca().set_xticks(member_x_values)
        plt.gca().set_xticklabels(member_labels)

        node_x_values = []

        force_shear_failure_glue_tab_list = []

        node_labels = []

        for i in range(len(self.truss_list)):
            node_x_values.append(i)

            force_shear_failure_glue_tab_list.append(self.failures["force_shear_failure_glue_tab"][i])
            node_labels.append(i)
        
        node_forces_list = self.node_forces
        
        plt.figure()
        plt.bar(node_x_values, force_shear_failure_glue_tab_list, color="tab:red", zorder=1)
        if override_node_forces_positive_list != None and override_node_forces_negative_list != None:
            plt.bar(node_x_values, override_node_forces_positive_list, color="tab:green", zorder=2)
            plt.bar(node_x_values, override_node_forces_negative_list, color="tab:green", zorder=2)
        else:
            plt.bar(node_x_values, node_forces_list, color="tab:green", zorder=2)
        plt.title("Node Glue Tab Shear Force")
        plt.ylabel("Shear Force (N)")

        plt.xticks(rotation=90)
        plt.gca().set_xticks(node_x_values)
        plt.gca().set_xticklabels(node_labels)
        
    def draw_truss(self):
        '''
        Draw image of truss using locations and neighbour connections of each TrussNode
        '''

        max_mag = 0
        for key in self.internal_forces:
            max_mag = max(abs(self.internal_forces[key]), max_mag)
        
        SCALE = 0.5
        COLOUR_MULT = 255 / max_mag
        turtle.speed(0)
        turtle.colormode(255)
        turtle.Screen().tracer(0)
        
        max_x = 0
        for node in self.truss_list:
            max_x = max(node.loc[0], max_x)
        
        x_offset = max_x * SCALE / 2

        for node in self.truss_list:
            for neighbour in node.neighbours:
                turtle.penup()
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
        '''
        BeamSlice class (each slice represents a small, discrete segment of the beam)
        Used as members of the Beam class
        Used to calculate changing beam properties along the length of the beam

        dx: Relative x spacing between this slice and the previous one [float]
        last_slice: Previous slice reference (used to propagate data through the beam) [BeamSlice]

        young_mod: Young's modulus (MPa) [float]
        poisson: Poisson's ratio [float]

        tension_ult: Tensile strength (MPa) [float]
        comp_ult: Compressive strength (MPa) [float]
        shear_ult: Shear strength (MPa) [float]
        glue_shear_ult: Glue shear strength (MPa) [float]

        cs_properties: Cross-sectional properties [dict]
        '''
        self.x = 0

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

        self.propagate_last_slice_properties()
    
    def propagate_last_slice_properties(self):
        '''
        Propagate properties from last slice to this slice
        This allows data to propagate through the beam, from slice to slice
        '''

        if self.last_slice != None:
            self.x = self.last_slice.x + self.dx
            self.shear_force = self.last_slice.shear_force
            self.b_moment = self.last_slice.b_moment
            self.curvature = self.last_slice.curvature
        else:
            self.x = 0
            self.shear_force = 0
            self.b_moment = 0
            self.curvature = 0
    
    def set_external_force(self, force):
        '''
        Set external load on this slice

        force: Load force [Vector]
        '''

        self.external_force = force

    def add_external_force(self, force):
        '''
        Add external load to the existing load on this slice

        force: Load force [Vector]
        '''

        self.external_force.add_vector(force)
    
    def update_centroidal_axis_and_second_moment_of_area(self):
        '''
        Update the height of the centroidal axis and the second moment of area for this slice
        This is calculated for each slice to account for varying cross sections across the beam

        Assumes the following general cross sectional shape - the size of each component can be set by changing the slice's cs_properties
        
        --------------  <-- top plate
          |-      -|<-- glue tab
          |        |    <-- web
          |-      -|
          ----------    <-- bottom plate
        '''

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
        

    def update_Qs(self):
        '''
        Update first moment of area from locations of interest to top/bottom of beam slice
        
        Locations of interest:
        - Height of centroidal axis
        - Height of top glue tab connection
        - Height of bottom glue tab connection
        '''

        self.update_centroidal_axis_and_second_moment_of_area()
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
        '''
        Update internal properties of this beam slice

        Properties:
        - Shear Force
        - Bending Moment
        - Curvature
        - Flexural Stress at Top
        - Flexural Stress at Bottom
        - Shear Stress at Centroidal Axis
        - Shear Stress at Top Glue Joint
        - Shear Stress at Bottom Glue Joint
        '''

        self.propagate_last_slice_properties()
        self.update_centroidal_axis_and_second_moment_of_area()
        
        self.external_force.set_dir_comps_by_current_mag()
        self.shear_force += self.external_force.get_y_mag()
        
        self.b_moment += self.shear_force * self.dx
        self.curvature = self.b_moment / (self.young_mod * self.second_moment_area)

        self.flexural_stress_top = -self.b_moment * (self.bottom_thickness + self.web_height + self.top_thickness - self.centroidal_axis) / self.second_moment_area
        self.flexural_stress_bot = -self.b_moment * (-self.centroidal_axis) / self.second_moment_area

        self.update_Qs()

        self.shear_stress_max = (self.shear_force * self.q_from_centroid) / (self.second_moment_area * self.web_thickness * 2)
        self.shear_stress_glue_top = (self.shear_force * self.q_from_glue_top) / (self.second_moment_area * self.web_thickness * 2)
        self.shear_stress_glue_bot = (self.shear_force * self.q_from_glue_bot) / (self.second_moment_area * self.web_thickness * 2)
    
    def moment_tension_failure_wall(self):
        '''
        Calculate bending moment required to cause a flexural stress resulting in a tensile failure at this beam slice
        '''

        if self.b_moment < 0:
            failure_moment = (self.tension_ult * self.second_moment_area) / self.centroidal_axis
        else:
            failure_moment = (self.tension_ult * self.second_moment_area) / (self.bottom_thickness + self.web_height + self.top_thickness - self.centroidal_axis)
        
        return failure_moment

    def moment_compression_failure_wall(self):
        '''
        Calculate bending moment required to cause a flexural stress resulting in a compressive failure at this beam slice
        '''

        if self.b_moment < 0:
            failure_moment = (self.comp_ult * self.second_moment_area) / (self.bottom_thickness + self.web_height + self.top_thickness - self.centroidal_axis)
        else:
            failure_moment = (self.comp_ult * self.second_moment_area) / self.centroidal_axis
        
        return failure_moment
    
    def force_shear_failure_wall(self):
        '''
        Calculate shear force required to cause a shear stress failure at the centroidal axis of this beam slice
        '''

        self.update_Qs()
        failure_force = self.shear_ult * 2 * self.web_thickness * self.second_moment_area / self.q_from_centroid

        return failure_force

    def force_shear_failure_glue_top(self):
        '''
        Calculate shear force required to cause a shear stress failure at the top glue joint of this beam slice
        '''

        self.update_Qs()
        failure_force = self.glue_shear_ult * 2 * self.web_thickness * self.second_moment_area / self.q_from_glue_top
        
        return failure_force
    
    def force_shear_failure_glue_bot(self):
        '''
        Calculate shear force required to cause a shear stress failure at the bottom glue joint of this beam slice
        '''

        self.update_Qs()
        failure_force = self.glue_shear_ult * 2 * self.web_thickness * self.second_moment_area / self.q_from_glue_bot
        
        return failure_force
    
    def moment_buckling_compressive_top(self):
        '''
        Calculate bending moment required to cause local buckling in the top plate of this beam slice
        '''

        stress_crit = (4*math.pi**2*self.young_mod) / (12*(1-self.poisson**2)) * (self.top_thickness/self.web_spacing)**2
        failure_moment = stress_crit * self.second_moment_area / ((self.bottom_thickness + self.web_height + self.top_thickness/2) - self.centroidal_axis)
        
        return failure_moment

    def moment_buckling_compressive_bot(self):
        '''
        Calculate bending moment required to cause local buckling in the bottom plate of this beam slice
        '''

        stress_crit = (4*math.pi**2*self.young_mod) / (12*(1-self.poisson**2)) * (self.bottom_thickness/self.bottom_width)**2
        failure_moment = stress_crit * self.second_moment_area / (self.bottom_thickness/2 - self.centroidal_axis)
        
        return failure_moment

    def moment_buckling_flanges(self):
        '''
        Calculate bending moment required to cause local buckling in the top flanges of this beam slice
        '''

        stress_crit = (0.425*math.pi**2*self.young_mod) / (12*(1-self.poisson**2)) * (self.top_thickness/((self.top_width - self.bottom_width)/2))**2
        failure_moment = stress_crit * self.second_moment_area / ((self.bottom_thickness + self.web_height + self.top_thickness/2) - self.centroidal_axis)
        
        return failure_moment

    def moment_buckling_webs_flexural(self):
        '''
        Calculate bending moment required to cause a flexural stress resulting in local buckling in the webs of this beam slice
        '''

        stress_crit = (6*math.pi**2*self.young_mod) / (12*(1-self.poisson**2)) * (self.web_thickness/(self.bottom_thickness + self.web_height - self.centroidal_axis))**2
        failure_moment = stress_crit * self.second_moment_area / (self.centroidal_axis - self.bottom_thickness)
        
        return failure_moment

class Beam():
    def __init__(self, young_mod, poisson, tension_ult, comp_ult, shear_ult, glue_shear_ult, segments, diaphragms, resolution, a_rxn_loc=None, b_rxn_loc=None):
        '''
        Beam data structure consisting of BeamSlices
        Contains BeamSlice elements representing discrete segments of the beam
        Various discrete beam slices allow for simple/consistent calculations given a varying cross section

        young_mod: Young's modulus (MPa) [float]
        poisson: Poisson's ratio [float]
        tension_ult: Tensile strength (MPa) [float]
        comp_ult: Compressive strength (MPa) [float]
        shear_ult: Shear strength (MPa) [float]
        glue_shear_ult: Glue shear strength (MPa) [float]

        segments: Cross-sectional segment properties -> also implies beam length [dict]
        diaphragms: Diaphragm locations [list]

        resolution: Number of beam slices (higher -> more accuracy) [int]

        a_rxn_loc: (x,y) pivot reaction point A_x, A_y [tuple]
        b_rxn_loc: (x,y) roller reaction point B_y [tuple]
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
        self.min_FOS = {}

        self.resolution = resolution
        self.dx = self.length / self.resolution
    
    def gen_beam(self):
        '''
        Generate Beam data structure consisting of linked BeamSlices
        Number of beam slices is determined by resolution of beam
        '''

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
        '''
        Update external forces on Beam by applying an external force to the beam slice at the location of each force
        '''

        for i in range(self.resolution):
            self.beam_list[i].set_external_force(Vector(0,0,0))
            for loc in self.external_forces:
                if loc[0] >= i*self.dx and loc[0] < (i+1)*self.dx:
                    self.beam_list[i].add_external_force(self.external_forces[loc])
            
            self.beam_list[i].update_internal_properties()

    def add_load(self, loc, load):
        '''
        Add external load to Beam

        loc: (x,y) location of force [tuple]
        load: Load force [Vector]
        '''

        if loc in self.external_forces:
            self.external_forces[loc].add_vector(load)
        else:
            self.external_forces[loc] = load

        self.update_external_forces()
    
    def set_rxn_locs(self, a_rxn_loc, b_rxn_loc):
        '''
        Set locations of reaction forces

        a_rxn_loc: (x,y) pivot reaction point A_x, A_y [tuple]
        b_rxn_loc: (x,y) roller reaction point B_y [tuple]
        '''

        self.a_rxn_loc = a_rxn_loc
        self.b_rxn_loc = b_rxn_loc
    
    def update_rxn_forces(self):
        '''
        Update reaction forces on Beam based on current acting external loads
        '''

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
    
    def calc_deflection(self, loc):
        '''
        Calculate deflection at point of interest using MAT #2
        Assumes no known horizontal tangents
        Finds tangential deflection between supports, and tangential deflection between support A and point of interest, then uses geometry to find deflection at point of interest

        point_loc: (x,y) point of interest [tuple]

        *Assumes that a_rxn_loc(x) < point_loc(x) < b_rxn_loc(x)
        '''

        tangential_deflection_at_b = 0
        tangential_deflection_at_loc = 0

        for i in range(self.resolution):
            if self.a_rxn_loc[0] < i*self.dx:
                if self.b_rxn_loc[0] >= i*self.dx:
                    tangential_deflection_at_b += i*self.dx * self.beam_list[i].curvature * self.dx
                if loc[0] >= i*self.dx:
                    tangential_deflection_at_loc += i*self.dx * self.beam_list[i].curvature * self.dx

        displacement_at_loc = (loc[0] - self.a_rxn_loc[0]) * (tangential_deflection_at_b / (self.b_rxn_loc[0] - self.a_rxn_loc[0])) - tangential_deflection_at_loc

        return displacement_at_loc       
    
    def forces_shear_buckling_webs(self):
        '''
        Calculate shear force required to cause local buckling in each span of the web between diaphragms
        '''

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
                return

            initial_slice.update_Qs()
            stress_crit = (5*math.pi**2*self.young_mod)/(12*(1-self.poisson**2))*((initial_slice.web_thickness/(width*2/3))**2 +((initial_slice.web_thickness)/(initial_slice.web_height))**2)
            
            failure_force = stress_crit * 2 * initial_slice.web_thickness * initial_slice.second_moment_area / initial_slice.q_from_centroid

            web_shear_failure_forces[span] = failure_force
        
        return web_shear_failure_forces
    
    def update_internal_properties(self):
        '''
        Update internal properties of Beam across its length
        Properties vary across the length of the beam

        Properties:
        - Shear Force
        - Bending Moment
        - Curvature
        '''

        self.sfd = []
        self.bmd = []
        self.curv = []

        for i in range(self.resolution):
            self.sfd.append(self.beam_list[i].shear_force)
            self.bmd.append(self.beam_list[i].b_moment)
            self.curv.append(self.beam_list[i].curvature)

    def update_failures(self):
        '''
        Update forces required to cause failure by each failure mode for every single beam slice

        Forces:
        - Shear Force
        - Internal Bending Moment

        Failure modes:
        - Beam Web Shear Failure
        - Beam Top Joint Shear Failure
        - Beam Bottom Joint Shear Failure
        - Beam Web Local Buckling Shear Failure
        - Beam Web Flexural Tensile Yield Failure
        - Beam Web Flexural Compressive Yield Failure
        - Beam Top Plate Local Buckling Failure
        - Beam Bottom Plate Local Buckling Failure
        - Beam Top Flange Local Buckling Failure
        - Beam Web Flexural Local Buckling Failure
        '''

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
    
    def update_min_FOS(self):
        '''
        Update FOS for each failure mode by iterating through the forces in every beam slice and finding the beam slices with the lowest FOS in each mode

        Forces:
        - Shear Force
        - Internal Bending Moment

        Failure modes:
        - Beam Web Shear Failure
        - Beam Top Joint Shear Failure
        - Beam Bottom Joint Shear Failure
        - Beam Web Local Buckling Shear Failure
        - Beam Web Flexural Tensile Yield Failure
        - Beam Web Flexural Compressive Yield Failure
        - Beam Top Plate Local Buckling Failure
        - Beam Bottom Plate Local Buckling Failure
        - Beam Top Flange Local Buckling Failure
        - Beam Web Flexural Local Buckling Failure
        '''

        sfd_modes = ["force_shear_failure_wall", "force_shear_failure_glue_top", "force_shear_failure_glue_bot", "forces_buckling_webs_shear"]

        for mode in sfd_modes:
            self.min_FOS[mode] = None
            for i in range(len(self.sfd)):
                if self.sfd[i] != 0:
                    if self.failures[mode][i] / self.sfd[i] > 0:
                        if self.min_FOS[mode] == None:
                            self.min_FOS[mode] = self.failures[mode][i] / self.sfd[i]
                        else:
                            self.min_FOS[mode] = min(self.failures[mode][i] / self.sfd[i], self.min_FOS[mode])
        
        bmd_modes = ["moment_tension_failure_wall", "moment_compression_failure_wall", "moment_buckling_compressive_top", "moment_buckling_compressive_bot", "moment_buckling_flanges", "moment_buckling_webs_flexural"]

        for mode in bmd_modes:
            self.min_FOS[mode] = None
            for i in range(len(self.bmd)):
                if self.bmd[i] != 0:
                    if self.failures[mode][i] / self.bmd[i] > 0:
                        if self.min_FOS[mode] == None:
                            self.min_FOS[mode] = self.failures[mode][i] / self.bmd[i]
                        else:
                            self.min_FOS[mode] = min(self.failures[mode][i] / self.bmd[i], self.min_FOS[mode])
    
    def update_min_envelope_FOS(self, property_envelope):
        '''
        Update FOS for each failure mode given beam property envelopes by iterating through and finding the beam slice with the lowest FOS at any point in time

        Forces:
        - Shear Force
        - Internal Bending Moment

        Failure modes:
        - Beam Web Shear Failure
        - Beam Top Joint Shear Failure
        - Beam Bottom Joint Shear Failure
        - Beam Web Local Buckling Shear Failure
        - Beam Web Flexural Tensile Yield Failure
        - Beam Web Flexural Compressive Yield Failure
        - Beam Top Plate Local Buckling Failure
        - Beam Bottom Plate Local Buckling Failure
        - Beam Top Flange Local Buckling Failure
        - Beam Web Flexural Local Buckling Failure
        '''

        sfd_modes = ["force_shear_failure_wall", "force_shear_failure_glue_top", "force_shear_failure_glue_bot", "forces_buckling_webs_shear"]

        for mode in sfd_modes:
            self.min_FOS[mode] = None
            for i in range(len(self.sfd)):
                for sfd_envelope in (property_envelope["beam_sfd_positive_envelope"], property_envelope["beam_sfd_negative_envelope"]):
                    if sfd_envelope[i] != 0:
                        if self.failures[mode][i] / sfd_envelope[i] > 0:
                            if self.min_FOS[mode] == None:
                                self.min_FOS[mode] = self.failures[mode][i] / sfd_envelope[i]
                            else:
                                self.min_FOS[mode] = min(self.failures[mode][i] / sfd_envelope[i], self.min_FOS[mode])
        
        bmd_modes = ["moment_tension_failure_wall", "moment_compression_failure_wall", "moment_buckling_compressive_top", "moment_buckling_compressive_bot", "moment_buckling_flanges", "moment_buckling_webs_flexural"]

        for mode in bmd_modes:
            self.min_FOS[mode] = None
            for i in range(len(self.bmd)):
                for bmd_envelope in (property_envelope["beam_bmd_positive_envelope"], property_envelope["beam_bmd_negative_envelope"]):
                    if bmd_envelope[i] != 0:
                        if self.failures[mode][i] / bmd_envelope[i] > 0:
                            if self.min_FOS[mode] == None:
                                self.min_FOS[mode] = self.failures[mode][i] / bmd_envelope[i]
                            else:
                                self.min_FOS[mode] = min(self.failures[mode][i] / bmd_envelope[i], self.min_FOS[mode])
    
    def draw_internal_properties(self, override_positive_sfd_list=None, override_negative_sfd_list=None, override_positive_bmd_list=None, override_negative_bmd_list=None, override_positive_curv_list=None, override_negative_curv_list=None):
        '''
        Plot all internal forces along beam along with corresponding FOS for each failure mode

        Forces:
        - Shear Force
        - Internal Bending Moment

        Failure modes:
        - Beam Web Shear Failure
        - Beam Top Joint Shear Failure
        - Beam Bottom Joint Shear Failure
        - Beam Web Local Buckling Shear Failure
        - Beam Web Flexural Tensile Yield Failure
        - Beam Web Flexural Compressive Yield Failure
        - Beam Top Plate Local Buckling Failure
        - Beam Bottom Plate Local Buckling Failure
        - Beam Top Flange Local Buckling Failure
        - Beam Web Flexural Local Buckling Failure
        '''

        self.update_internal_properties()
        self.update_failures()

        x_values = []

        for i in range(self.resolution):
            x_values.append(i*self.dx)
        
        sfd_modes = ["force_shear_failure_wall", "force_shear_failure_glue_top", "force_shear_failure_glue_bot", "forces_buckling_webs_shear"]
        plt.figure()
        for i in range(4):
            plt.subplot(2,2,i+1)
            if override_positive_sfd_list != None and override_positive_sfd_list != None:
                plt.plot(x_values, override_positive_sfd_list, color="tab:blue")
                plt.plot(x_values, override_negative_sfd_list, color="tab:blue")
            else:
                plt.plot(x_values, self.sfd, color="tab:blue")
            plt.plot(x_values, self.failures[sfd_modes[i]], color="tab:orange")
            plt.title("SFD vs. " +sfd_modes[i])
            plt.ylabel("Shear Force (N)")

        bmd_modes = ["moment_tension_failure_wall", "moment_compression_failure_wall", "moment_buckling_compressive_top", "moment_buckling_compressive_bot", "moment_buckling_flanges", "moment_buckling_webs_flexural"]
        plt.figure()
        for i in range(6):
            plt.subplot(3,2,i+1)
            plt.gca().invert_yaxis()
            if override_positive_bmd_list != None and override_positive_bmd_list != None:
                plt.plot(x_values, override_positive_bmd_list, color="tab:blue")
                plt.plot(x_values, override_negative_bmd_list, color="tab:blue")
            else:
                plt.plot(x_values, self.bmd, color="tab:blue")
            plt.plot(x_values, self.failures[bmd_modes[i]], color="tab:orange")
            plt.title("BMD vs. " +bmd_modes[i])
            plt.ylabel("Internal Moment (Nmm)")

        plt.figure()
        plt.gca().invert_yaxis()
        if override_positive_curv_list != None and override_positive_curv_list != None:
            plt.plot(x_values, override_positive_curv_list, color="tab:blue")
            plt.plot(x_values, override_negative_curv_list, color="tab:blue")
        else:
            plt.plot(x_values, self.curv, color="tab:blue")
        plt.title("Curvature")
        plt.ylabel("Curvature (rad/mm)")

class Bridge():
    def __init__(self, truss_flattop, young_mod, poisson, tension_ult, comp_ult, shear_ult, glue_shear_ult, truss_cs_properties, truss_modified_cs_members, truss_modified_tab_length_nodes, beam_segments, beam_diaphragms, resolution=1000):
        '''
        Bridge class to provide a common interface with a bridge consisting of a Beam and Truss
        Beam and Truss are treated independently, but Bridge ensures that common physical properties are kept consistent between the two
        '''
        
        self.truss = Truss(truss_flattop, young_mod, poisson, tension_ult, comp_ult, glue_shear_ult, truss_cs_properties, truss_modified_cs_members, truss_modified_tab_length_nodes)
        self.beam = Beam(young_mod, poisson, tension_ult, comp_ult, shear_ult, glue_shear_ult, beam_segments, beam_diaphragms, resolution)

        self.periodic = None
        self.truss_triangle_count = None
        self.truss_diagonal = None
        self.truss_relative_base_xs = None
        self.truss_height = None
        self.truss_flipped = False

        self.min_FOS = {}
    
    def reset(self):
        '''
        Reset Bridge to initial state by reinitializing and regenerating Beam and Truss
        '''

        self.beam.__init__(self.beam.young_mod, self.beam.poisson, self.beam.tension_ult, self.beam.comp_ult, self.beam.shear_ult, self.beam.glue_shear_ult, self.beam.segments, self.beam.diaphragms, self.beam.resolution, self.beam.a_rxn_loc, self.beam.b_rxn_loc)
        self.truss.__init__(self.truss.flattop, self.truss.young_mod, self.truss.poisson, self.truss.tension_ult, self.truss.comp_ult, self.truss.glue_shear_ult, self.truss.cs_properties, self.truss.modified_cs_members, self.truss.modified_tab_length_nodes, self.truss.a_rxn_loc, self.truss.b_rxn_loc)
        
        if self.periodic:
            self.gen_bridge_periodic_truss(self.truss_triangle_count, self.truss_diagonal, self.truss_flipped)
        else:
            self.gen_bridge_nonperiodic_truss(self.truss_relative_base_xs, self.truss_height, self.truss_flipped)

    def gen_bridge_periodic_truss(self, triangle_count, diagonal, flipped=False):
        '''
        Generate Beam and Truss where Truss is a periodic and symmetrical Warren truss

        triangle_count: Number of triangles in Warren truss [int]
        diagonal: (x, y) distance between first 2 nodes [tuple]
        flipped: True if bridge is upside-down [bool]
        '''

        self.beam.gen_beam()
        self.truss.gen_periodic_warren_truss(triangle_count, diagonal, flipped)

        self.periodic = True
        self.truss_triangle_count = triangle_count
        self.truss_diagonal = diagonal
        self.truss_flipped = flipped

    def gen_bridge_nonperiodic_truss(self, relative_base_xs, height, flipped=False):
        '''
        Generate Beam and Truss where Truss is a nonperiodic and asymmetrical Warren truss

        relative_base_xs: Relative x distances between the nodes at the base of the truss [list]
        height: height of bridge [float]
        flipped: True if bridge is upside-down [bool]
        '''

        self.beam.gen_beam()
        self.truss.gen_nonperiodic_warren_truss(relative_base_xs, height, flipped)

        self.periodic = False
        self.truss_relative_base_xs = relative_base_xs
        self.truss_height = height
        self.truss_flipped = flipped

    def set_rxn_locs(self, a_rxn_loc, b_rxn_loc):
        '''
        Set locations of reaction forces on Bridge (Truss + Beam)

        a_rxn_loc: (x,y) pivot reaction point A_x, A_y [tuple]
        b_rxn_loc: (x,y) roller reaction point B_y [tuple]
        '''

        self.truss.set_rxn_locs(a_rxn_loc, b_rxn_loc)
        self.beam.set_rxn_locs(a_rxn_loc, b_rxn_loc)
    
    def update_bridge(self):
        '''
        Update reaction forces and internal forces/properties for Bridge (Truss + Beam)
        '''

        self.truss.update_rxn_forces()
        self.beam.update_rxn_forces()
        
        self.truss.update_internal_forces()
        self.truss.update_failures()
        self.beam.update_internal_properties()
        self.beam.update_failures()

    def add_load(self, loc, load):
        '''
        Add external load to Bridge (Truss + Beam)

        loc: (x,y) location of force [tuple]
        load: Load force [Vector]
        '''

        self.truss.add_load(loc, Vector(load.x_dir, load.y_dir, load.get_mag() * 0.5))
        self.beam.add_load(loc, Vector(load.x_dir, load.y_dir, load.get_mag()))
    
    def update_bridge_FOS(self):
        '''
        Update FOS for Truss and Beam
        '''
        self.truss.update_min_FOS()
        self.beam.update_min_FOS()

    def display_bridge_FOS(self):
        '''
        Display FOS for each failure mode by consolidating the lowest FOS from Truss and Beam
        '''

        self.update_bridge()

        for mode in self.truss.min_FOS:
            self.min_FOS[mode] = self.truss.min_FOS[mode]
        
        for mode in self.beam.min_FOS:
            self.min_FOS[mode] = self.beam.min_FOS[mode]
        
        print("Min. FOS by Failure Mode:")
        for mode in self.min_FOS:
            print(str(mode) +":", self.min_FOS[mode])

    def draw_properties(self, override_properties_and_forces=None):
        if override_properties_and_forces != None:
            self.truss.draw_internal_forces(override_properties_and_forces["truss_internal_forces_positive_envelope"], override_properties_and_forces["truss_internal_forces_negative_envelope"], override_properties_and_forces["truss_node_forces_positive_envelope"], override_properties_and_forces["truss_node_forces_negative_envelope"])
            self.beam.draw_internal_properties(override_properties_and_forces["beam_sfd_positive_envelope"], override_properties_and_forces["beam_sfd_negative_envelope"], override_properties_and_forces["beam_bmd_positive_envelope"], override_properties_and_forces["beam_bmd_negative_envelope"], override_properties_and_forces["beam_curv_positive_envelope"], override_properties_and_forces["beam_curv_negative_envelope"])
        else:
            self.truss.draw_internal_forces()
            self.beam.draw_internal_properties()

def simulate_train_load(bridge, x, train_weight, train_wheel_spacing_rel):
    '''
    Simulate train loading at variable positions by resetting bridge and applying external loads at locations of train wheels
    
    bridge: Reference to Bridge [Bridge]
    x: x location of front of train [float]
    train_weight: Total weight of train (N) [float]
    train_wheel_spaing_rel: Relative x spacing of wheels from front to back of train [list]
    '''

    bridge.reset()

    bridge_height = abs(bridge.truss.truss_list[1].loc[1] - bridge.truss.truss_list[0].loc[1])
    bridge_length = bridge.beam.length

    train_load_locs = []
    for i in range(len(train_wheel_spacing_rel)):
        load_loc = x
        for j in range(i+1):
            load_loc -= train_wheel_spacing_rel[j]
        
        if load_loc >= 0 and load_loc <= bridge_length:
            train_load_locs.append(load_loc)

    for load_loc in train_load_locs:
        bridge.add_load((load_loc,bridge_height), Vector(0,-1,train_weight / len(train_wheel_spacing_rel)))
    
    bridge.update_bridge()

    truss_internal_forces = bridge.truss.internal_forces
    truss_node_forces = bridge.truss.node_forces

    beam_sfd = bridge.beam.sfd
    beam_bmd = bridge.beam.bmd
    beam_curv = bridge.beam.curv

    return truss_internal_forces, truss_node_forces, beam_sfd, beam_bmd, beam_curv

def run_train(bridge, train_weight, train_wheel_spacing_rel, resolution):
    '''
    Simulate running train across bridge by simulating train loading at various discrete positions
    All of the "worst case" bridge properties are saved, so that the final result is an envelope of the bridge properties closest to failure (which do not necessarily all occur at once)

    bridge: Reference to Bridge [Bridge]
    train_weight: Total weight of train (N) [float]
    train_wheel_spaing_rel: Relative x spacing of wheels from front to back of train [list]
    resolution: Number of discrete steps to simulate train loading over [int]
    '''

    train_length = sum(train_wheel_spacing_rel)

    dx = (bridge.beam.length + train_length) / resolution

    truss_internal_forces_positive_envelope = {}
    truss_internal_forces_negative_envelope = {}
    truss_node_forces_positive_envelope = []
    truss_node_forces_negative_envelope = []

    beam_sfd_positive_envelope = []
    beam_sfd_negative_envelope = []
    beam_bmd_positive_envelope = []
    beam_bmd_negative_envelope = []
    beam_curv_positive_envelope = []
    beam_curv_negative_envelope = []

    for i in range(resolution):
        truss_internal_forces, truss_node_forces, beam_sfd, beam_bmd, beam_curv = simulate_train_load(bridge, i*dx, train_weight, train_wheel_spacing_rel)
        for key in truss_internal_forces:
            if key not in truss_internal_forces_positive_envelope:
                truss_internal_forces_positive_envelope[key] = truss_internal_forces[key]
            elif truss_internal_forces[key] > truss_internal_forces_positive_envelope[key]:
                truss_internal_forces_positive_envelope[key] = truss_internal_forces[key]
            
            if key not in truss_internal_forces_negative_envelope:
                truss_internal_forces_negative_envelope[key] = truss_internal_forces[key]
            elif truss_internal_forces[key] < truss_internal_forces_negative_envelope[key]:
                truss_internal_forces_negative_envelope[key] = truss_internal_forces[key]
        
        for j in range(len(truss_node_forces)):
            if j > len(truss_node_forces_positive_envelope)-1:
                truss_node_forces_positive_envelope.append(truss_node_forces[j])
            elif truss_node_forces[j] > truss_node_forces_positive_envelope[j]:
                truss_node_forces_positive_envelope[j] = truss_node_forces[j]
            
            if j > len(truss_node_forces_negative_envelope)-1:
                truss_node_forces_negative_envelope.append(truss_node_forces[j])
            elif truss_node_forces[j] < truss_node_forces_negative_envelope[j]:
                truss_node_forces_negative_envelope[j] = truss_node_forces[j]
        
        for j in range(len(beam_sfd)):
            if j > len(beam_sfd_positive_envelope)-1:
                beam_sfd_positive_envelope.append(beam_sfd[j])
            elif beam_sfd[j] > beam_sfd_positive_envelope[j]:
                beam_sfd_positive_envelope[j] = beam_sfd[j]

            if j > len(beam_sfd_negative_envelope)-1:
                beam_sfd_negative_envelope.append(beam_sfd[j])
            elif beam_sfd[j] < beam_sfd_negative_envelope[j]:
                beam_sfd_negative_envelope[j] = beam_sfd[j]
        
        for j in range(len(beam_bmd)):
            if j > len(beam_bmd_positive_envelope)-1:
                beam_bmd_positive_envelope.append(beam_bmd[j])
            elif beam_bmd[j] > beam_bmd_positive_envelope[j]:
                beam_bmd_positive_envelope[j] = beam_bmd[j]
            
            if j > len(beam_bmd_negative_envelope)-1:
                beam_bmd_negative_envelope.append(beam_bmd[j])
            elif beam_bmd[j] < beam_bmd_negative_envelope[j]:
                beam_bmd_negative_envelope[j] = beam_bmd[j]
        
        for j in range(len(beam_curv)):
            if j > len(beam_curv_positive_envelope)-1:
                beam_curv_positive_envelope.append(beam_curv[j])
            elif beam_curv[j] > beam_curv_positive_envelope[j]:
                beam_curv_positive_envelope[j] = beam_curv[j]
            
            if j > len(beam_curv_negative_envelope)-1:
                beam_curv_negative_envelope.append(beam_curv[j])
            elif beam_curv[j] < beam_curv_negative_envelope[j]:
                beam_curv_negative_envelope[j] = beam_curv[j]
        
        print("Finished", str(i+1) +"/" +str(resolution), "iterations")
    
    force_and_property_envelopes =     {"truss_internal_forces_positive_envelope": truss_internal_forces_positive_envelope, 
                                        "truss_internal_forces_negative_envelope": truss_internal_forces_negative_envelope, 
                                        "truss_node_forces_positive_envelope": truss_node_forces_positive_envelope, 
                                        "truss_node_forces_negative_envelope": truss_node_forces_negative_envelope, 
                                        "beam_sfd_positive_envelope": beam_sfd_positive_envelope,
                                        "beam_sfd_negative_envelope": beam_sfd_negative_envelope,
                                        "beam_bmd_positive_envelope": beam_bmd_positive_envelope,
                                        "beam_bmd_negative_envelope": beam_bmd_negative_envelope,
                                        "beam_curv_positive_envelope": beam_curv_positive_envelope,
                                        "beam_curv_negative_envelope": beam_curv_negative_envelope}
    
    bridge.truss.update_min_envelope_FOS(force_and_property_envelopes)
    bridge.beam.update_min_envelope_FOS(force_and_property_envelopes)

    return force_and_property_envelopes

'''
Editable Bridge Properties
--------------------------
'''
bridge_height = 100

'''
Truss member cross sectional properties
'''
truss_cs_properties =  {"width": 75, "thickness": 1.27, "tab_width": 65, "tab_length": 10}

'''
Truss members with modified cross sectional properties (reinforced members)
(node_a_index, node_b_index) represents the member from node_a to node_b
'''
truss_modified_cs_members = {   (1,3): {"width": 100, "thickness": 1.27},                   
                                (3,5): {"width": 100, "thickness": 1.27},
                                (9,11): {"width": 100, "thickness": 1.27},
                                (11,13): {"width": 100, "thickness": 1.27},
                                (13,15): {"width": 100, "thickness": 1.27},
                                (15,17): {"width": 100, "thickness": 1.27},

                                (5,7): {"width": 88.97, "thickness": 1.9627}, 
                                (7,9): {"width": 88.97, "thickness": 1.9627}, 
                                (9,10): {"width": 67.5, "thickness": 2.1166},
                                (10,12): {"width": 67.5, "thickness": 2.1166},
                                (11,12): {"width": 67.5, "thickness": 2.1166},
                                (12,14): {"width": 75, "thickness": 2.54},
                                (14,15): {"width": 75, "thickness": 2.54},
                                (14,16): {"width": 75, "thickness": 2.54},
                                (16,17): {"width": 67.5, "thickness": 2.1166},
                                (16,18): {"width": 75, "thickness": 2.54},
                                (17,18): {"width": 61.5384, "thickness": 4.1275}}

'''
Truss joint nodes with modified glue tab lengths (reinforced tabs)
'''
truss_modified_tab_length_nodes = {4: 15, 5: 20, 6:20, 7: 25, 8: 20, 11: 15, 12: 20, 13: 20, 14: 20, 15: 25, 16: 25, 17: 25}

'''
Relative distances between the base of each truss "triangle" (adjacent pair of diagonal members) from the start to end of the truss
Allows for definition of asymmetrical and non periodic trusses
'''
truss_relative_base_xs = [40,232.5,232.5,130,195,195,110,115,40]

'''
Discrete beam segments and their cross sectional properties
(x_start, x_end) represents the beam segment from x_start to x_end
Cross sectional properties defined in a corresponding dictionary for each segment
'''
# beam_segments =    [((0, 1290),        {"top_thickness": 1.27,
#                                         "top_width": 100,
#                                         "web_thickness": 1.27,
#                                         "web_height": 75-1.27*2,
#                                         "web_spacing": 80-1.27*2,
#                                         "tab_width": 10,
#                                         "bottom_thickness": 1.27,
#                                         "bottom_width": 80})
                    # ]
# beam_diaphragms = [0, 30, 535, 565, 1045, 1075, 1060, 1090]

beam_segments =    [((0, 388.75),      {"top_thickness": 1.27,          #0-5
                                        "top_width": 100,
                                        "web_thickness": 1.27,
                                        "web_height": 100-1.27*2,
                                        "web_spacing": 75-1.27*2,
                                        "tab_width": 10,
                                        "bottom_thickness": 1.27,
                                        "bottom_width": 75}),

                    ((388.75, 732.5),  {"top_thickness": 1.9627,        #5-9
                                        "top_width": 88.97,
                                        "web_thickness": 1.27,
                                        "web_height": 100-1.27*3,
                                        "web_spacing": 75-1.27*2,
                                        "tab_width": 10,
                                        "bottom_thickness": 1.27,
                                        "bottom_width": 75}),

                    ((732.5, 830),     {"top_thickness": 1.27,          #9-10
                                        "top_width": 100,
                                        "web_thickness": 1.27,
                                        "web_height": 100-1.27*2,
                                        "web_spacing": 75-1.27*2,
                                        "tab_width": 10,
                                        "bottom_thickness": 1.27,
                                        "bottom_width": 75}),
                    
                    ((830, 1025),      {"top_thickness": 1.27,          #10-12
                                        "top_width": 100,
                                        "web_thickness": 1.27,
                                        "web_height": 100-1.27*2,
                                        "web_spacing": 75-1.27*2,
                                        "tab_width": 10,
                                        "bottom_thickness": 2.1166,
                                        "bottom_width": 67.5}),

                    ((1025, 1135),     {"top_thickness": 1.27,          #12-14
                                        "top_width": 100,
                                        "web_thickness": 1.27,
                                        "web_height": 100-1.27*2,
                                        "web_spacing": 75-1.27*2,
                                        "tab_width": 10,
                                        "bottom_thickness": 2.286,
                                        "bottom_width": 69.4444}),
                    
                    ((1135, 1290),     {"top_thickness": 1.27,          #14-18
                                        "top_width": 100,
                                        "web_thickness": 1.27,
                                        "web_height": 100-1.27*2,
                                        "web_spacing": 75-1.27*2,
                                        "tab_width": 10,
                                        "bottom_thickness": 2.54,
                                        "bottom_width": 75})
                    ]

'''
Static Main Code
----------------
'''
beam_diaphragms = []

'''
Place approximated vertical diaphragms at horizontal midpoint of each diagonal truss member
'''
x=0
for i in truss_relative_base_xs:
    beam_diaphragms.append(x + i/4)
    beam_diaphragms.append(x + 3*i/4)
    x += i

'''
Generate bridge with defined physical parameters

young_mod: 4000 MPa
poisson: 0.2
tension_ult: 30 MPa
comp_ult: 6 MPa
shear_ult: 4 MPa
glue_shear_ult: 2 MPa
'''
bridge = Bridge(True, 4000, 0.2, 30, 6, 4, 2, truss_cs_properties, truss_modified_cs_members, truss_modified_tab_length_nodes, beam_segments, beam_diaphragms, 1000)
bridge.gen_bridge_nonperiodic_truss(truss_relative_base_xs, bridge_height+10, False)

bridge.set_rxn_locs((15,0), (1075,0))
bridge.reset()

'''
Simulate either train loading or point loading by commentin
'''
# bridge_force_and_property_envelopes = run_train(bridge, 400, (52, 176, 164, 176, 164, 176), 10)
# bridge.display_bridge_FOS()
# bridge.draw_properties(bridge_force_and_property_envelopes)

bridge.add_load((565,bridge_height+10), Vector(0,-1,1500))
bridge.add_load((1265,bridge_height+10), Vector(0,-1,1500))
bridge.update_bridge()
bridge.update_bridge_FOS()
bridge.display_bridge_FOS()
bridge.draw_properties()

'''
Calculate deflection due to applied loads
'''
print("Truss Deflection:", bridge.truss.calc_displacement((530, 0), Vector(0, 1, 1)))
print("Beam Deflection:", bridge.beam.calc_deflection((530, 0)))
plt.show()

bridge.truss.draw_truss()