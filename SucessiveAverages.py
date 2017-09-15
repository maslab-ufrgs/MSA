'''
Created on 09/01/2015

@author: gabriel
'''
import string

# represents a node in the graph
class Node(object):
    def __init__(self, name):
        self.name = name     # name of the node
        self.dist = 1000000  # distance to this node from start node
        self.prev = None     # previous node to this node
        self.flag = 0        # access flag

# represents an edge in the graph
class Edge(object):
    def __init__(self, u, v, capacity, free_flow, alpha, beta):
        self.name = "%s-%s" % (u,v)
        
        self.start = u
        self.end = v
        
        self.capacity = capacity
        self.free_flow = free_flow
        self.alpha = alpha
        self.beta = beta
        
        self.flow = 0
        self.aux_flow = 0
                
        self.update_cost()
        
    # BPR value-delay function (VDF)
    def update_cost(self):
        self.cost = self.free_flow * (1 + self.beta * ((self.flow / self.capacity) ** self.alpha))

# read a text file and generate the graph according to declarations
def generateGraph(graph_file):
    N = []
    N_names = []
    E = []
    fname = open(graph_file, "r")
    line = fname.readline()
    line = line[:-1]
    while line:
        taglist = string.split(line)
        
        u = taglist[0]
        if u not in N_names:
            N.append(Node(u))
            N_names.append(u)
        v = taglist[1]
        if v not in N_names:
            N.append(Node(v))
            N_names.append(v)
        
        E.append(Edge(u, v, float(taglist[2]), float(taglist[3]), float(taglist[4]), float(taglist[5])))
        
        line = fname.readline()
        line = line[:-1]
    fname.close()
    return N, E

# reset graph's variables to default
def resetGraph(N, E):
    for node in N:
        node.dist = 1000000.0
        node.prev = None
        node.flag = 0

# returns the smallest node in N but not in S
def pickSmallestNode(N):
    minNode = None
    for node in N:
        if node.flag == 0:
            minNode = node
            break
    if minNode == None:
        return minNode
    for node in N:
        if node.flag == 0 and node.dist < minNode.dist:
            minNode = node
    return minNode

# returns the edges list of node u
def pickEdgesList(u, E):
    uv = []
    for edge in E:
        if edge.start == u.name:
            uv.append(edge)
    return uv

# Dijkstra's shortest path algorithm
def dijkstra(N, E, origin, destination, ignoredEdges):
    
    #reset the graph (so as to discard information from previous runs)
    resetGraph(N, E)
    
    # set origin node distance to zero, and get destination node
    dest = None
    for node in N:
        if node.name == origin:
            node.dist = 0
        if node.name == destination:
            dest = node
    
    u = pickSmallestNode(N)
    while u != None:
        u.flag = 1
        uv = pickEdgesList(u, E)
        n = None
        for edge in uv:
            
            # avoid ignored edges
            if edge in ignoredEdges:
                continue
            
            # take the node n
            for node in N:
                if node.name == edge.end:
                    n = node
                    break
            if n.dist > u.dist + edge.cost:
                n.dist = u.dist + edge.cost
                n.prev = u
        
        u = pickSmallestNode(N)
        # stop when destination is reached
        if u == dest:
            break
    
    # generate the final path
    S = []
    u = dest
    while u.prev != None:
        S.insert(0,u)
        u = u.prev
    S.insert(0,u)
    
    return S

# print vertices and edges
def printGraph(N, E):
    print('vertices:')
    for node in N:
        previous = node.prev
        if previous == None:
            print(node.name, node.dist, previous)
        else:
            print(node.name, node.dist, previous.name)
    print('edges:')
    for edge in E:
        print(edge.start, edge.end, edge.cost)

# get the directed edge from u to v
def getEdge(E, u, v):
    for edge in E:
        if edge.start == u and edge.end == v:
            return edge
    return None

# get the directed edge from u to v
def getNode(N, n):
    for node in N:
        if node.name == n:
            return node
    return None

# calculate path P's cost
def calcPathLength(P, N, E):
    if type(P[0]) is Edge:
        P = getPathAsNodes(P, N, E)
    length = 0
    prev = None
    for node in P:
        if prev != None:
            length += getEdge(E, prev.name, node.name).cost
        prev = node
    
    return length

# calculate path P's cost
def getPathAsEdges(P, E):
    path = []
    prev = None
    for node in P:
        if prev != None:
            path.append(getEdge(E, prev.name, node.name))
        prev = node
    
    return path

def getPathAsNodes(P, N, E):
    path = []
    path.append(getNode(N, P[0].start))
    for edge in P:
        path.append(getNode(N, edge.end))
    return path

# print the path S
def printPath(path, N, E):
    #S = N
    #if type(path[0]) is Edge:
    #    S = E
    strout = ''
    for e in path:
        if strout != '':
            strout += ' - '
        strout += e.name
        
    print "%g = %s" % (calcPathLength(path, N, E), strout)

def pathToStr(path, N, E):
    
    if type(path[0]) is Node:
        path = getPathAsEdges(path, E) 
    
    strout = ""
    
    for e in path:
        if strout != '':
            strout += ' - '
        strout += e.name
    
    return strout

def process_OD_matrix(od_file):
    
    od_file = open(od_file, "r")
    
    OD_matrix = {}
    
    for line in od_file:
        [o, d, demand] = line.split("\t")
        if o != d:
            OD_matrix[od_list_to_str(o, d)] = float(demand)
    
    od_file.close()
    
    return OD_matrix

def od_str_to_list(od):
    return od.split("|")

def od_list_to_str(o,d):
    return "%s|%s" % (o,d)

# method of successive averages
def run_MSA(its, N, E, OD_matrix):
    
    # initial value for phi
    phi = 1.0
    
    # a nested dictionary data structure to store, for each OD pair, 
    # its routes, and, for each route, its edges and flows 
    # an entry can be said a 4-uple: (OD, route string, route, flow) 
    od_routes_flow = {od : {} for od in OD_matrix}
    
    # iterations
    for n in xrange(1, its+1):
        
        # update phi
        phi = 1.0 / n
        
        # clear auxiliary flow of all links
        for e in E:
            e.aux_flow = 0
        
        # calculate auxiliary flow based on a all-or-nothing assignment
        min_routes = {}
        for od in OD_matrix:
            [o,d] = od_str_to_list(od)
            
            # compute shortest route
            route = getPathAsEdges(dijkstra(N, E, o, d, []), E)
            route_str = pathToStr(route, N, E)
            
            # store min route of this od pair
            min_routes[od] = [route_str, route]
            
            # if the min route is not in the od routes' list, add it
            if route_str not in od_routes_flow[od]:
                od_routes_flow[od][route_str] = [route, 0]
        
        # calculate current flow of all links
        for od in OD_matrix:
            for route in od_routes_flow[od]: 
                
                # route flow on previous iteration
                vna = od_routes_flow[od][route][1]
                
                # auxiliary route flow (0 if not the current best route)
                fa = 0
                if route == min_routes[od][0]:
                    fa = OD_matrix[od]
                
                # route flow of current iteration
                vna = max((1 - phi) * vna + phi * fa, 0)
                
                # update flows and costs
                od_routes_flow[od][route][1] = vna
                for e in od_routes_flow[od][route][0]:
                    e.aux_flow += vna
                
        for e in E:
            e.flow = e.aux_flow
            e.update_cost()
        
        # print (if desired) the main values calculated during current iteration
        if False:
            print "---- it %i ---------------------"%n
            print "phi %f"%phi
            for od in OD_matrix:
                print "%s (best route: %s)"%(od, min_routes[od][0])
                for route in od_routes_flow[od]:
                    print "\t%s"%route
                    fa = 0
                    if route == min_routes[od][0]:
                        fa = OD_matrix[od]
                    print "\tfa=%i"%fa
                    print "\tvna=%i"%od_routes_flow[od][route][1]
                    print ""
        
    # print the final assignment    
    evaluate_assignment(OD_matrix, od_routes_flow)

def evaluate_assignment(OD_matrix, od_routes_flow):
    # header
    print "od\troute\tflow\ttravel time\tdeviations"
    
    sum_tt = 0.0
    sum_deviations = 0
    delta_top = 0.0
    delta_bottom = 0.0
    
    for od in od_routes_flow:
        
        aux = []
        min_cost = 99999999999999
        
        # calculate some information of each route
        for route in od_routes_flow[od]:
            
            # calculate cost of the route
            cost = 0.0
            for e in od_routes_flow[od][route][0]:
                cost += e.cost
                sum_tt += e.cost * od_routes_flow[od][route][1]
            cost = round(cost * 100) / 100 #to handle imprecise double representation (IEEE 754, via stackoverflow.com/questions/15625556)
            
            # store minimum route cost of current OD pair
            if cost < min_cost:
                min_cost = cost
            
            # store the values in a temporary data structure to allow
            # the calculations of the "deviations from best" measure
            aux.append([od, route, od_routes_flow[od][route][1], cost])
            
        # read the temporary data structure and print the results
        for e in aux:
            
            # calculate the "deviations from best" measure
            deviations = 0
            if e[3] > min_cost:
                deviations = e[2]
                sum_deviations += deviations
            
            # update the top part of delta equation
            delta_top += e[2] * (e[3] - min_cost)
            
            # print
            print "%s\t%s\t%f\t%f\t%i" % (e[0], e[1], e[2], e[3], deviations)
        
        # update the bottom part of delta equation
        delta_bottom += OD_matrix[od]# * min_cost
        
    # print overall results
    print "Average travel time: %f min"%(sum_tt / sum([x for x in OD_matrix.values()]))
    print "Deviations: %i"%(sum_deviations)
    #print "%s: %.10f"%(u'\u03B4', (delta_top / delta_bottom))
    print "AEC: %.10f"%(delta_top / delta_bottom)
    
if __name__ == '__main__':
    
    net = "Sioux_Falls"
    
    # read graph from file
    N, E = generateGraph('graphs/%s.net'%net)
    
    # read list of OD_matrix-pairs
    OD_matrix = process_OD_matrix('graphs/%s.od'%net)
    
    # run MSA
    run_MSA(1000, N, E, OD_matrix)
