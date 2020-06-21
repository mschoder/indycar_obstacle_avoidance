

def extractLocalGraph(time2goMap: graph, planning_horizon: float, state: vector) -> graph:
    '''
    Extracts subgraph from pre-computed time-to-go map based on car's current location on the
    track and the planning horizon distance. Inserts a virtual node linking all points on the
    graph at the final station distance (s) as the terminus.
    '''

def deleteObstacleNodes(localGraph: graph, obstacleStateInfo: array) -> graph:
    '''
    For dynamic obstacles (moving vehicles), node removal is performed.
    For static obstacles (debris, stopped vehicles), edge removal is performed.
    Depends on specific implemetation of collision_check() TBD.
    '''

def determinePrimitives(feasibleGraph: graph, )
    '''
    Determine which action primitives are possible given the feasible graph
    '''

def getSubgraphs(feasibleGraph: graph, actionPrimitives: list) -> list(graph):
    '''
    Breaks apart feasible graph in to N subgraphs for determined set of primitives
    '''

def findShortestPath(subGraph: graph) -> list(nodes):
    '''
    Search implementation to find minimum-cost path to virtual segment end node
    '''

def getC2trajectory(path: list(nodes), veh_parameters) -> list(nodes):
    '''
    Modify shortest path with additional curvature continuity constraints
    '''