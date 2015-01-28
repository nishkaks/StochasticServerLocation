#-------------------------------------------------------------------------------
# Name:        ScenDecompMain.py
# Description: Main program for solving SSLP instances using scenario
#              decomposition algorithm for 0-1 stochastic programs
#
# Author:      Nishka K.S
#
# Created:     13/05/2014
# Parameter:   SMPS filename, example - "sslp_5_25_50"
# Usage:       ScenDecompMain.py sslp_5_25_50
#-------------------------------------------------------------------------------


from gurobipy import *
from ReadSTOFileSMPS import read_sto_file_smps
import time

# turn off logging
loggingOffOn = 0

# Used to fix the x variables in the problem based on the scenario suproblem
# solutions
def Fix_BoundsLB_HB(m,LOWER,UPPER):
    for i in range(nScen):
        varNum = 0
        for j in range(nServer):
            varNum = j + 1
            varTemp = m[i].getVarByName("x_"+str(varNum))
            varTemp.LB = LOWER[j]
            varTemp.UB = UPPER[j]
        m[i].update()


# Read file name from command line argument
if len(sys.argv) < 2:
    print 'Missing command line argument. Example Usage: ScenDecompMain.py sslp_5_25_50'
    quit()

file_name = sys.argv[1]

# -------------------
# READ INSTANCE DATA
# -------------------
# Get number of servers, clients and scenarios
nServer = int(file_name.split('_')[1])
nClient = int(file_name.split('_')[2])
nScen = int(file_name.split('_')[3])

# probability of each scenario - equal
ps = 1/nScen

# Read core file using gurobi built in function
# (changed .cor to .mps so that gurobi can read the mps file)
# This gives the Master problem structure
# Copy this and modify the h matrix for each scenario
tempMaster = read(file_name+".mps")

# custom utility to read the .sto file for SMPS format
# matrix h stores which of the clients materialize in each scenario
# h - one row per scenario, one column per client
h = read_sto_file_smps(file_name+".sto",nScen,nClient,nServer)

# --------------------------------------------
# BUILD COPIES OF THE MODEL FOR EACH SCENARIO
# --------------------------------------------

# scenario subproblems in the lower bounding scheme
# Relaxed with continous second stage variables
# integer cuts would be added to these models
numVar = (nServer * nClient) + nServer
SubprobScen = []
for i in range(nScen):
    SubprobScen.append(Model("SubprobScen"+str(i)))
    SubprobScen[i] = tempMaster.copy()
    SubprobScen[i].params.OutputFlag = loggingOffOn
    # Modify the h variable for each scenario
    constrNum = 1 + nServer
    for j in range(nClient):
        constrNum += 1
        constrttemp = SubprobScen[i].getConstrByName("c"+str(constrNum))
        constrttemp.setAttr("rhs", h[i][j])

    xVariable = SubprobScen[i].getVars()
    for k in range(nServer,numVar):
        xVariable[k].vtype = GRB.CONTINUOUS

    SubprobScen[i].update()


# Relaxed scenario subproblems with continuous second stage variables
# Used in upper bounding scheme, these would not have cuts added to them
MasterScenRelax = []
for i in range(nScen):
    MasterScenRelax.append(Model("MasterScenRelax"+str(i)))
    MasterScenRelax[i] = SubprobScen[i].copy()
    MasterScenRelax[i].params.OutputFlag = loggingOffOn

# Scenario subproblems with binary second stage variables
# used in upper bounding scheme
MasterScen = []
for i in range(nScen):
    MasterScen.append(Model("MasterScen"+str(i)))
    MasterScen[i] = SubprobScen[i].copy()
    MasterScen[i].params.OutputFlag = loggingOffOn
    xVariable = MasterScen[i].getVars()
    for k in range(nServer,numVar):
        xVariable[k].vtype = GRB.INTEGER


#SubprobScen[1].write("Scen1.lp")

# --------------------------------------------
# SCENARIO DECOMPOSITION ALGORITHM
# --------------------------------------------

# we use lambda = 0, hence there is no update lamdba step

upperbound = float('+inf')
lowerbound = float('-inf')
infinity = float('+inf')
integerCutNum = [0] * nScen
numIterations = 0

start_time = time.clock()

while(upperbound > lowerbound):
    numIterations +=1

    #Terminate if too long
    if (numIterations >=100):
        break

    # LOWER BOUNDING STEP:
    lowerbound = 0
    sHat = []
    for i in range(nScen):
        SubprobScen[i].optimize()
        if SubprobScen[i].status == GRB.status.INFEASIBLE:
            lowerbound = float('+inf')
        else:
            lowerbound += SubprobScen[i].objVal
            xVal = []
            varNum = 0
            for j in range(nServer):
                varNum = j+1
                varTemp = SubprobScen[i].getVarByName("x_"+str(varNum))
                xVal.append(varTemp.x)
            sHat.append(xVal)

    #Remove duplicate solutions to avoid duplicate cuts
    sHat = [list(t) for t in set(map(tuple, sHat))]

    if (lowerbound == infinity):
        break

    # UPPER BOUNDING STEP:
    for x in sHat:
        # start by solving the relaxed problem
        # Fix x values
        Fix_BoundsLB_HB(MasterScenRelax,x,x)
        ubar = 0
        for i in range(nScen):
            MasterScenRelax[i].optimize()
            ubar += MasterScenRelax[i].objVal
        # Solve integer subproblems only if it is promising
        if (ubar < upperbound):
            Fix_BoundsLB_HB(MasterScen,x,x)
            u = 0
            for i in range(nScen):
                MasterScen[i].optimize()
                u += MasterScen[i].objVal
            if (upperbound > u):
                upperbound = u
                xopt = x

    # Add integer cuts to scenario subproblems
    for i in range(nScen):
        xVariable = SubprobScen[i].getVars()
        integerCuts = {}
        for x in sHat:
            integerCuts[integerCutNum[i]] = SubprobScen[i].addConstr(
              quicksum(x[r] + (1- (2 * x[r])) * xVariable[r] for r in range(nServer)) >= 1,
                       'integerCuts_%s' % integerCutNum[i])
            integerCutNum[i] +=1
        SubprobScen[i].update()


#-----------------------------
# PRINT RESULTS
#-----------------------------
print "Total time taken      = ",time.clock() - start_time, "seconds"
print 'Optimal Value         = ',upperbound/nScen
print 'Number of Iterations  = ',numIterations

#SubprobScen[0].write("Scen1.lp")
