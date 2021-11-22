class Vector():
    def __init__(self, x_dir, y_dir, mag):
        self.x_dir = x_dir
        self.y_dir = y_dir
        self.mag = mag

    def get_opposite(self):
        return Vector(-self.x_dir, -self.y_dir, self.mag)

class TrussNode():
    def __init__(self, loc=(0,0)):
        self.loc = loc
        self.neighbours = []
        self.external_forces = []
        self.internal_forces = {}

    def new_neighbour_rel(self, rel):
        neighbour_loc = (self.loc[0] + rel[0], self.loc[1] + rel[1])
        self.neighbours.append(TrussNode(neighbour_loc))

        self.internal_forces[self.neighbours[-1]] = []

        return self.neighbours[-1]
    
    def add_neighbour(self, neighbour):
        self.neighbours.append(neighbour)

        self.internal_forces[self.neighbours[-1]] = []
    
    def add_external_force(self, force):
        self.external_forces.append(force)

    def calc_internal_forces(self):
        for neighbour in self.neighbours:
            if neighbour.internal_forces[self]:
                self.internal_forces[neighbour] = neighbour.internal_forces[self].get_opposite()
        
        



#simple truss
node_1 = TrussNode((0,0))
node_2 = node_1.new_neighbour_rel((6,0))
node_3 = node_2.new_neighbour_rel((6,0))
node_4 = node_3.new_neighbour_rel((6,0))

node_5 = TrussNode((3,4))
node_6 = node_5.new_neighbour_rel((1,0))
node_7 = node_6.new_neighbour_rel((1,0))

node_1.add_neighbour(node_5)
node_5.add_neighbour(node_2)
node_2.add_neighbour(node_6)
node_6.add_neighbour(node_3)
node_3.add_neighbour(node_7)
node_7.add_neighbour(node_4)