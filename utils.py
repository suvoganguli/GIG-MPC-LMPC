import numpy as np
from matplotlib.patches import Polygon
import pickle

def getLine(x1, y1, x2, y2):

    a = y2 - y1 # LHS
    b = -(x2 - x1)  # LHS
    c = x1 * y2 - x2 * y1  # RHS

    return a, b, c


def getPosIdx(E, N, path, posIdx0 = None):

    posIdx = {'number': -1}

    AR_Path = path.alongPathLines.AR
    BR_Path = path.alongPathLines.BR
    CR_Path = path.alongPathLines.CR
    AL_Path = path.alongPathLines.AL
    BL_Path = path.alongPathLines.BL
    CL_Path = path.alongPathLines.CL

    D1 = path.acrossPathLines.D1
    E1 = path.acrossPathLines.E1
    F1 = path.acrossPathLines.F1
    D2 = path.acrossPathLines.D2
    E2 = path.acrossPathLines.E2
    F2 = path.acrossPathLines.F2

    nSections = len(D1)

    if posIdx0 == None:
        kvec = range(nSections)
    else:
        kvec = np.arange(posIdx0['number'], nSections, 1)


    for k in kvec:
        inbox = insideBox(E, N, AR_Path[k], BR_Path[k], CR_Path[k], AL_Path[k], BL_Path[k], CL_Path[k],
                          D1[k], E1[k], F1[k], D2[k], E2[k], F2[k])

        if inbox == True:
            posIdx['number'] = k
            return posIdx

    return posIdx


def insideBox(x,y,AR,BR,CR,AL,BL,CL,D1,E1,F1,D2,E2,F2):

    chk1bool = False
    chk2bool = False
    chk3bool = False
    chk4bool = False

    chk1 = D1 * x + E1 * y - F1
    chk2 = D2 * x + E2 * y - F2
    chk3 = AR * x + BR * y - CR
    chk4 = AL * x + BL * y - CL

    myeps = 1e-3

    if chk1 >= -myeps:
        chk1bool = True
    if chk2 <= myeps:
        chk2bool = True
    if chk3 <= 0:
        chk3bool = True
    if chk4 >= 0:
        chk4bool = True

    if (chk1bool) and \
            (chk2bool) and \
            (chk3bool) and \
            (chk4bool):
        return True
    else:
        return False


def insideBox2(x, y, x_rect, y_rect):

    x1 = x_rect[0]
    x2 = x_rect[1]
    x3 = x_rect[2]
    x4 = x_rect[3]

    y1 = y_rect[0]
    y2 = y_rect[1]
    y3 = y_rect[2]
    y4 = y_rect[3]

    AR, BR, CR = getLine(x2, y2, x3, y3) # right-bottom to right-top
    AL, BL, CL = getLine(x1, y1, x4, y4) # left-bottom to left-top
    D1, E1, F1 = getLine(x2, y2, x1, y1) # right-bottom to left-bottom
    D2, E2, F2 = getLine(x2, y2, x3, y3) # right-top to left-top

    # print(AR, BR, CR)
    # print(AL, BL, CL)
    # print(D1, E1, F1)
    # print(D2, E2, F2)

    inside = insideBox(x, y, AR, BR, CR, AL, BL, CL, D1, E1, F1, D2, E2, F2)

    return inside


def insideBox3(x,y,x1,y1,x2,y2):

    AR, BR, CR = getLine(x2, y1, x2, y2) # right-bottom to right-top
    AL, BL, CL = getLine(x1, y1, x1, y2) # left-bottom to left-top
    D1, E1, F1 = getLine(x2, y1, x1, y1) # right-bottom to left-bottom
    D2, E2, F2 = getLine(x2, y2, x1, y2) # right-top to left-top

    inside = insideBox(x, y, AR, BR, CR, AL, BL, CL, D1, E1, F1, D2, E2, F2)

    return inside


def insideGap(x,y,AR,BR,CR,AL,BL,CL):

    chk1bool = False
    chk2bool = False

    chk1 = AR * x + BR * y - CR
    chk2 = AL * x + BL * y - CL

    myeps = 1e-3

    if chk1 >= -myeps:
        chk1bool = True
    if chk2 <= myeps:
        chk2bool = True

    if (chk1bool) and (chk2bool):
        return True
    else:
        return False

def changeAxis(x,y,theta):
    C = [[np.cos(theta), np.sin(theta)],[-np.sin(theta),np.cos(theta)]]
    X = np.array([x,y]).T
    Y = np.matmul(C,X)
    print(X)
    return Y[0], Y[1]


def intersect(laneLines, acrossLines, idx):
    a1 = laneLines.A_Lane[idx]
    b1 = laneLines.B_Lane[idx]
    c1 = laneLines.C_Lane[idx]

    a2 = acrossLines.D1[idx]
    b2 = acrossLines.E1[idx]
    c2 = acrossLines.F1[idx]

    determinant = a1*b2 - a2*b1

    if (determinant == 0):
        return []
    else:
        x = (c1*b2 - c2*b1) / determinant
        y = (a1*c2 - a2*c1) / determinant

    return x, y

def intersect2(a1,b1,c1,a2,b2,c2):

    determinant = a1*b2 - a2*b1

    if (determinant == 0):
        return []
    else:
        x = (c1*b2 - c2*b1) / determinant
        y = (a1*c2 - a2*c1) / determinant

    return x, y

def rotateRectangle(Ec, Nc, E, N, theta): # del_chi = - theta

    E = E - Ec
    N = N - Nc

    ERot = np.zeros(4)
    NRot = np.zeros(4)

    C = np.array([[np.cos(theta), +np.sin(theta)],[-np.sin(theta), np.cos(theta)]])

    for k in range(len(E)):
        p = np.array([E[k], N[k]])
        p = p[:,None]
        pRot = np.dot(C,p)
        ERot[k] = pRot[0]
        NRot[k] = pRot[1]

    ERot = E + Ec
    NRot = N + Nc

    return ERot, NRot

def getPatch(Ec,Nc,W,L,theta,fc):

    # create object with heading = 0 deg

    E1 = -W/2
    E2 = W/2
    E3 = E2
    E4 = E1
    N1 = -L/2
    N2 = N1
    N3 = L/2
    N4 = N3
    E = np.array([E1,E2,E3,E4])
    N = np.array([N1,N2,N3,N4])

    # Rotate object

    ERot = np.zeros(4)
    NRot = np.zeros(4)
    C = np.array([[np.cos(theta), +np.sin(theta)],[-np.sin(theta), np.cos(theta)]])
    for k in range(len(E)):
        p = np.array([E[k], N[k]])
        p = p[:,None]
        pRot = np.dot(C,p)
        ERot[k] = pRot[0] + Ec
        NRot[k] = pRot[1] + Nc

    vertices = np.array([ERot, NRot])
    polygon = Polygon(vertices.T, facecolor=fc, alpha=0.75)

    return polygon

def getEllipse(W, L):

    theta_sol = np.arctan(L / W)
    a = W / (2 * np.cos(theta_sol))
    b = L / (2 * np.sin(theta_sol))

    return a, b


def distance(p1,p2):
    x1 = p1[0]
    y1 = p1[1]
    x2 = p2[0]
    y2 = p2[1]

    d = np.sqrt((x1-x2)**2 + (y1-y2)**2)
    return d

def shiftRotate(vec1,dvec,Chi):

    dcm = np.array([
        [np.cos(Chi), np.sin(Chi)],
        [-np.sin(Chi), np.cos(Chi)]
    ])
    vec1 = vec1[:, None]
    dvec = dvec[:, None]
    vec2 = vec1 + np.matmul(dcm, dvec)

    return vec2

def rotate(vec1,Chi):

    dcm = np.array([
        [np.cos(Chi), np.sin(Chi)],
        [-np.sin(Chi), np.cos(Chi)]
    ])
    vec1 = np.squeeze(vec1)
    vec1 = vec1[:, None]

    vec2 = np.matmul(dcm, vec1)

    return np.squeeze(vec2)

def createGrid(gridSize, lengthSpace, widthSpace, heightSpace):

    # always set height at the middle of the space to run Laplacian Planner without error
    height = heightSpace / 2  # ft

    nE = widthSpace / gridSize
    nN = lengthSpace / gridSize
    nU = heightSpace / gridSize
    nU_low = nU

    class grid():
        def __init__(self):
            self.nE = nE
            self.nN = nN
            self.nU = nU
            self.nU_low = nU_low
            self.height = height
            self.gridSize = gridSize # ft
            self.lengthSpace = lengthSpace # ft
            self.widthSpace = widthSpace  # ft
            self.heightSpace = heightSpace  # ft
            pass

    return grid


def getColumns(inFile, delim=" ", header=True):     # delim="\t"
    """
    Get columns of data from inFile. The order of the rows is respected

    :param inFile: column file separated by delim
    :param header: if True the first line will be considered a header line
    :returns: a tuple of 2 dicts (cols, indexToName). cols dict has keys that
    are headings in the inFile, and values are a list of all the entries in that
    column. indexToName dict maps column index to names that are used as keys in
    the cols dict. The names are the same as the headings used in inFile. If
    header is False, then column indices (starting from 0) are used for the
    heading names (i.e. the keys in the cols dict)
    """
    cols = {}
    indexToName = {}
    for lineNum, line in enumerate(inFile):
        if lineNum == 0:
            headings = line.split(delim)
            i = 0
            for heading in headings:
                heading = heading.strip()
                if header:
                    cols[heading] = []
                    indexToName[i] = heading
                else:
                    # in this case the heading is actually just a cell
                    cols[i] = [heading]
                    indexToName[i] = i
                i += 1
        else:
            cells = line.split(delim)
            i = 0
            for cell in cells:
                cell = cell.strip()
                cols[indexToName[i]] += [cell]
                i += 1

    return cols, indexToName

# def savepkl(obj, file_pkl):
#     with open(file_pkl, 'wb') as f:
#         pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
#         #pickle.dump(obj, f, pickle.0)

def savepkl(obj, file_pkl):
    with open(file_pkl, 'wb') as f:
        pickle.dump(obj, f)

def loadpkl(file_pkl):
    with open(file_pkl, 'rb') as f:
        return pickle.load(f)

def makePathObj(pdata, path, obstacle):
    pathObj = { 'PathE': path.pathData.E,
                'PathN': path.pathData.N,
                'PathStartPoint': path.pathData.PathStartPoint,
                'PathEndPoint': path.pathData.PathEndPoint,
                'PathRightEndPointsE': path.pathData.PathRightEndPointsE,
                'PathRightEndPointsN': path.pathData.PathRightEndPointsN,
                'PathLeftEndPointsE': path.pathData.PathLeftEndPointsE,
                'PathLeftEndPointsN': path.pathData.PathLeftEndPointsN,
                'PathCenterEndPointsE': path.pathData.PathCenterEndPointsE,
                'PathCenterEndPointsN': path.pathData.PathCenterEndPointsN,
                'PathThetaEndpoints': path.pathData.Theta_endpoints,
                'PathDeltaYRoad':  pdata.delta_yRoad,
                'PathWidth': pdata.pathWidth,
                'ObstacleE': obstacle.E,
                'ObstacleN': obstacle.N,
                'ObstacleW': obstacle.w,
                'ObstacleL': obstacle.l,
                'ObstacleChi': obstacle.Chi,
            }

    return pathObj

def addCurrentPointToPath(path, startPoint, Chi):

    Theta = np.pi/2 - Chi

    path.pathData.E = np.append(startPoint[0], path.pathData.E)
    path.pathData.N = np.append(startPoint[1], path.pathData.N)
    path.pathData.PathStartPoint = startPoint

    RightEndPointE = startPoint[0] + path.pathWidth / 2 * np.sin(Theta)
    RightEndPointN = startPoint[1] - path.pathWidth / 2 * np.cos(Theta)
    LeftEndPointE = startPoint[0] - path.pathWidth / 2 * np.sin(Theta)
    LeftEndPointN = startPoint[1] + path.pathWidth / 2 * np.cos(Theta)

    path.pathData.PathRightEndPointsE = np.append(RightEndPointE, path.pathData.PathRightEndPointsE)
    path.pathData.PathRightEndPointsN = np.append(RightEndPointN, path.pathData.PathRightEndPointsN)
    path.pathData.PathLeftEndPointsE = np.append(LeftEndPointE, path.pathData.PathLeftEndPointsE)
    path.pathData.PathLeftEndPointsN = np.append(LeftEndPointN, path.pathData.PathLeftEndPointsN)

    path.pathData.PathCenterEndPointsE = np.append(startPoint[0], path.pathData.PathCenterEndPointsE)
    path.pathData.PathCenterEndPointsN = np.append(startPoint[1], path.pathData.PathCenterEndPointsN)

    path.pathData.Theta_endpoints = np.append(Theta, path.pathData.Theta_endpoints)

    return path


def obstacleDict_from_ClassInstance(obstacleClassInstance):
    obstacle = {
        'Present' : obstacleClassInstance.Present,
        'N': obstacleClassInstance.N,
        'E': obstacleClassInstance.E,
        'Chi': obstacleClassInstance.Chi,
        'N_corners': obstacleClassInstance.N_corners,
        'E_corners': obstacleClassInstance.E_corners,
        'w': obstacleClassInstance.w,
        'l': obstacleClassInstance.l,
        'sw': obstacleClassInstance.sw,
        'sl': obstacleClassInstance.sl,
        'sr': obstacleClassInstance.sr,
    }
    return obstacle

def obstacleClassInstance_from_Dict(obstacleDist):

    import obstacleData as od

    obstaclePresent = obstacleDist['Present']
    obstacleE       = obstacleDist['E']
    obstacleN       = obstacleDist['N']
    obstacleChi     = obstacleDist['Chi']
    obstacleWidth   = obstacleDist['w']
    obstacleLength  = obstacleDist['l']
    obstacleSafeWidth  = obstacleDist['sw']
    obstacleSafeLength = obstacleDist['sl']
    obstacleSafeRadius = obstacleDist['sr']

    obstacleClass = od.obstacleInfo(obstaclePresent, obstacleE, obstacleN, obstacleChi, obstacleWidth, obstacleLength,
                             obstacleSafeWidth, obstacleSafeLength, obstacleSafeRadius)
    obstacle = obstacleClass()
    return obstacle

def vehicleStop(T, x, mpciter, decelType, terminal_point, endPoint,
                lb_reachedGoal, lb_reachedNearGoal, zeroDistanceChange,
                t_slowDown_detected, tmeasure, V_cmd, lb_VTermSlowDown, lb_VdotValSlowDown, decel,
                t_slowDown, lb_VTerm, lb_VdotVal):
    # ----------------------------------------------------------------------------
    # Vehicle stopping

    breakLoop = False

    if decelType == 'Slow':

        # find detection time
        if (distance(terminal_point, endPoint) < lb_reachedNearGoal) and (t_slowDown_detected == False):
            t_slowDown = tmeasure
            t_slowDown_detected = True

        # slow down V_cmd near goal
        if t_slowDown_detected == True:

            V_cmd = V_cmd - decel * T

            lb_VTerm = lb_VTermSlowDown  # fps
            lb_VdotVal = lb_VdotValSlowDown  # fps2

            if (distance(terminal_point, endPoint) < lb_reachedGoal):
                print('Reached Goal')
                breakLoop = True

            if mpciter > 0:

                #print(distance(x[mpciter, 0:2], x[mpciter - 1, 0:2]))

                if distance(x[mpciter, 0:2], x[mpciter - 1, 0:2]) < zeroDistanceChange:
                    print('Stopped (distance moved close to zero)')
                    breakLoop = True

                if x[mpciter, 2] * x[mpciter - 1, 2] < 0:  # change in direction of terminal point
                    print('Stopped (V direction change)')
                    breakLoop = True


    elif decelType == 'Fast':
        if distance(terminal_point, endPoint) < lb_reachedNearGoal:
            print('Reached near goal')
            breakLoop = True

    return breakLoop, V_cmd, t_slowDown, t_slowDown_detected, lb_VTerm, lb_VdotVal