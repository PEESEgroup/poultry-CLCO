import pyomo.environ as pyo
import csv
from CLCO_Data import CLCO_Data
from sympy import symbols, evalf

FEEDSTOCK_SUPPLY = [23.25, 23.59, 28.56, 30.45, 52.72, 62.13, 76.13, 111.96, 125.89, 127.49, 151.68, 208.94,
                    272.69, 594.76, 758.82, 907.18, 907.45, 1285.48, 1356.12]
SUPPLY_X_COORD = [345, 696, 250, 572, 371, 312, 417,
                  350, 327, 555, 239, 558, 592, 303, 519, 190, 604, 334, 402]
SUPPLY_Y_COORD = [4695, 4535, 4682, 4680, 4764, 4747, 4874, 4738, 4722, 4937, 4765, 4583, 4773, 4682, 4618, 4740,
                  4955, 4799, 4762]

FACILITY_X_COORD = [i * 20 + 100 for i in range(int((700 - 100) / 20))]
FACILITY_Y_COORD = [i * 20 + 4440 for i in range(int((5000 - 4440) / 20))]

NPV = [[0 for i in range(6)] for j in range(6)]


# pre-compute the distances
def distance_calcs(A, max_dist):
    """
    Calculates the distance between a facility and a point producer
    :param A: dataset
    :param max_dist: maximum cutoff distance before transportation is assumed to be infeasible
    :return: a 3-D array containing all distances between every producer and the 2-D grid of potential facilities
    """
    return [[[999 if A.LOAD_TRANSIT_COST + A.OPEX['transit'] * (
            (SUPPLY_X_COORD[i] - x) ** 2 + (SUPPLY_Y_COORD[i] - y) ** 2) ** (1 / 2) > max_dist else
              A.LOAD_TRANSIT_COST + A.OPEX['transit'] * (
                      (SUPPLY_X_COORD[i] - x) ** 2 + (SUPPLY_Y_COORD[i] - y) ** 2) ** (1 / 2)
              for i in range(len(FEEDSTOCK_SUPPLY))] for x in FACILITY_X_COORD] for y in FACILITY_Y_COORD]


def FLP(t_cutoff_offset, t_dist, t_dist_no_logistics):
    """
    Model formulation for the facility location problem
    :param t_cutoff_offset: how many increments of 50km to set as the maximum transport distance from 450 km
    :param t_dist: list of transporation distances with cutoff.  [y][x][i] is the indexing scheme,
    with y the y position of the facility, x the x position of the facility, and i the producing facility number
    :param t_dist_no_logistics: list of transportation distances with no cutoff
    :return: N/A
    """
    # [y][x][i] is the indexing scheme for the transportation distances array
    transport_distances = t_dist

    M = pyo.ConcreteModel("FLP")

    M.num_facilities = pyo.Param(initialize=1, mutable=True)
    M.max_distance = pyo.Param(initialize=50, mutable=True)

    # set of locations for the facilities
    M.SupplyLocation = pyo.Set(initialize=[list(a) for a in zip(SUPPLY_X_COORD, SUPPLY_Y_COORD, FEEDSTOCK_SUPPLY)])
    M.Facility_X = pyo.Set(initialize=FACILITY_X_COORD)
    M.Facility_Y = pyo.Set(initialize=FACILITY_Y_COORD)

    # set of supply amounts
    M.Facility = pyo.Var(M.Facility_X, M.Facility_Y, within=pyo.Binary, initialize=0)
    M.Facility_Cost = pyo.Var(M.Facility_X, M.Facility_Y, within=pyo.NonNegativeReals, initialize=1)
    facility_upper = 10000
    M.Facility_Capacity = pyo.Var(M.Facility_X, M.Facility_Y, within=pyo.NonNegativeReals, bounds=(0, facility_upper),
                                  initialize=1)
    M.Amount_shipped = pyo.Var(M.Facility_X, M.Facility_Y, M.SupplyLocation, within=pyo.NonNegativeReals, initialize=1)
    M.Total_Costs = pyo.Var(within=pyo.NonNegativeReals, initialize=1)
    M.TransportCost = pyo.Var(M.Facility_X, M.Facility_Y, within=pyo.NonNegativeReals, initialize=0)

    M.const = pyo.ConstraintList()

    # limit the number of facilities in the run to t_cutoff_offset
    # not limited in this run
    M.const.add(expr=sum(M.Facility[x, y] for x in M.Facility_X for y in M.Facility_Y) == M.num_facilities)

    for k in M.SupplyLocation:
        # the amount of material leaving each node is equal to the supply at the node
        # amount shipped from supply node k to facility location l
        M.const.add(
            expr=sum(M.Amount_shipped[x, y, k[0], k[1], k[2]] for x in M.Facility_X for y in M.Facility_Y) == k[2])

    for x in M.Facility_X:
        for y in M.Facility_Y:
            # constraints based on the geometry of new york to speed up problem solve time
            if x < 390 and y > 4800:
                M.const.add(expr=M.Facility[x, y] == 0)
            if x < 500 and y < 4650:
                M.const.add(expr=M.Facility[x, y] == 0)
            if x > 620:
                M.const.add(expr=M.Facility[x, y] == 0)
            if y < 4530:
                M.const.add(expr=M.Facility[x, y] == 0)

            # calculate the capacity for each facility
            # originally sum(M.Amount_shipped[x, y, k] for k in M.SupplyLocation)
            # <= M.Facility_Capacity[x, y] * M.Facility[x, y], but applied glovers linearization
            M.const.add(expr=0 <= sum(M.Amount_shipped[x, y, k] for k in M.SupplyLocation))
            M.const.add(
                expr=sum(M.Amount_shipped[x, y, k] for k in M.SupplyLocation) <= facility_upper * M.Facility[x, y])
            M.const.add(expr=M.Facility_Capacity[x, y] - facility_upper * (1 - M.Facility[x, y]) <=
                             sum(M.Amount_shipped[x, y, k] for k in M.SupplyLocation))
            M.const.add(expr=sum(M.Amount_shipped[x, y, k] for k in M.SupplyLocation) <= M.Facility_Capacity[x, y])

            # impose a minimum limit on the size of each facility
            M.const.add(expr=M.Facility_Capacity[x, y] >= M.Facility[x, y] * 250)

            # implementing piecewise linear approximation with a dummy CAPEX equation that mimics pyrolysis
            q = symbols("q")
            thermochem = [i * 500 for i in range(21)]
            capex_cost = 60000 * (q ** .6)
            capex_est = [capex_cost.evalf(subs={q: x}) for x in thermochem]

            # Piecewise linear approximation to capital cost terms
            M.add_component(str(x) + "," + str(y), pyo.Piecewise(M.Facility_Cost[x, y], M.Facility_Capacity[x, y],
                                                                 pw_pts=thermochem,
                                                                 pw_constr_type='EQ',
                                                                 f_rule=capex_est,
                                                                 pw_repn='SOS2'))

            M.const.add(
                expr=M.TransportCost[x, y] == sum(t_dist_no_logistics[int((y - 4440) / 20)][int((x - 100) / 20)][
                                                      int(FEEDSTOCK_SUPPLY.index(k[2]))] * M.Amount_shipped[
                                                      x, y, k[0], k[1], k[2]] for k in M.SupplyLocation))

    # total costs depends upon the distance transported as well as the facility cost
    M.const.add(expr=M.Total_Costs == sum(M.Facility_Cost[x, y] for x in M.Facility_X for y in M.Facility_Y) +
                     sum(87 * transport_distances[int((y - 4440) / 20)][int((x - 100) / 20)][
                         int(FEEDSTOCK_SUPPLY.index(k[2]))] * M.Amount_shipped[x, y, k[0], k[1], k[2]]
                         for x in M.Facility_X for y in M.Facility_Y for k in
                         M.SupplyLocation))  # 87 to account for discounting over the 10 year time horizon as opposed to 120.

    M.obj = pyo.Objective(expr=M.Total_Costs, sense=pyo.minimize)

    # iterates through the different number of facilities that can exist at this transportation distance cutoff
    for i in range(6):
        # update and solve the model
        M.num_facilities = i + 1
        instance = M.create_instance()
        opt = pyo.SolverFactory('gurobi')
        if not i == 0:
            opt.options[
                'mipgap'] = .04 + t_cutoff_offset / 100 + i / 50 + i * i / 125  # for gurobi so as to not burn the computer down
        print(opt.solve(instance, tee=True))  # keepfiles = True

        # store the NPV of the location for later
        NPV[t_cutoff_offset][i] = pyo.value(instance.Total_Costs)
        print(NPV)

        # save information about the specifics of the model solution
        print_model(i + 1, instance, 450 - 50 * t_cutoff_offset)


def print_model(scenario, model, transit_dist):
    """
    Prints the relevant Total_Costs and amount of manure transported between any two locations
    :param scenario: The number of facilities
    :param model: the solved model
    :param transit_dist: the maximum transportation distance allowed to each facility
    :return: a csv file containing all pertinant information
    """
    with open("data/FLP/range" + str(transit_dist) + "facilities" + str(scenario) + 'transport.csv', 'w',
              encoding='UTF8',
              newline='') as f:
        write = csv.writer(f)

        # print out all the data that isn't part of the SOS constraints
        for v in model.component_objects(pyo.Var, active=True):
            if "SOS2" not in str(v):
                row = ["Variable", v]
                write.writerow(row)
                for index in v:
                    row = [v, index, pyo.value(v[index])]
                    write.writerow(row)


if __name__ == '__main__':
    # facility location problem is traditionally scenario number 4
    A = CLCO_Data(4)

    for i in range(6):
        '''
        There are 6 scenarios, corresponding to maximum transportation prices of $45/ton, $40/ton, $35/ton, $30/ton, 
        $25/ton, $20/ton, which corresponds to transportation circles of 450km, 400km, etc.
        '''
        FLP(i, distance_calcs(A, 45 - 5 * i), distance_calcs(A, 10000))

    print(NPV)
