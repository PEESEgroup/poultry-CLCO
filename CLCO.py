import os.path
import glob
import pyomo.environ as pyo
from pyomo.opt import TerminationCondition
import csv
from CLCO_Data import CLCO_Data
from sympy import symbols
import matplotlib.pyplot as plt
import math


def utopian(m, lca_midpoint, lca_type):
    """
    Finds the utopia and nadir points for the pareto front with NPV and any LCA subcategory
    :param m: the model
    :param lca_midpoint: the LCA midpoint type
    :param lca_type: the LCA type: ALCA, CLCA
    :return: two lists containing the utopia and nadir points
    """
    # step 1: normalize the objective functions
    # solve the model to get optimal objective function values
    utopia = []
    nadir = []

    # only one function can be active at a time
    m.Obj2.deactivate()
    m.Obj.activate()
    model = m
    opt = pyo.SolverFactory('gurobi')
    print(opt.solve(model))  # keepfiles = True

    if scenario < 4000:
        utopia.append(pyo.value(model.npv[0]))
        temp = pyo.value(model.total_LCA_midpoints[0, lca_type, lca_midpoint])
    else:
        utopia.append(pyo.value(model.total_LCA_midpoints[0, lca_type, 'climate change']))
        temp = pyo.value(model.total_LCA_midpoints[0, lca_type, 'eutrophication: freshwater'])

    # only one function can be active at a time
    m.Obj.deactivate()
    m.Obj2.activate()
    model = m
    opt = pyo.SolverFactory('gurobi')
    print(opt.solve(model))  # keepfiles = True

    # get the nadir and utopia points
    if scenario < 4000:
        utopia.append(pyo.value(model.total_LCA_midpoints[0, lca_type, lca_midpoint]))
        nadir.append(pyo.value(model.npv[0]))
    else:
        utopia.append(pyo.value(model.total_LCA_midpoints[0, lca_type, 'eutrophication: freshwater']))
        nadir.append(pyo.value(model.total_LCA_midpoints[0, lca_type, 'climate change']))

    nadir.append(pyo.value(temp))

    # add constraints on the anchor points
    if scenario < 4000:
        m.const.add(expr=m.npv[0] >= nadir[0])
        m.const.add(expr=m.total_LCA_midpoints[0, lca_type, lca_midpoint] <= nadir[1])
    else:
        m.const.add(expr=m.total_LCA_midpoints[0, lca_type, 'climate change'] <= nadir[0])
        m.const.add(expr=m.total_LCA_midpoints[0, lca_type, 'eutrophication: freshwater'] <= nadir[1])

    print("scalar numerator", (nadir[0] - utopia[0]), "nadir npv", nadir[0], "utopia npv", utopia[0])
    print("scalar denominator", (nadir[1] - utopia[1]), "nadir lca_midpoint", nadir[1], "utopia lca_midpoint",
          utopia[1])

    if (nadir[1] - utopia[1]) <= .000001:
        raise ValueError("No Pareto Front can be generated - no change in denominator")
    scalar = abs((nadir[0] - utopia[0]) / (nadir[1] - utopia[1]))
    print("scalar", scalar)

    # deactivate remaining objective function
    m.Obj2.deactivate()
    return utopia, nadir


def aws(M, divs, midpoint, utopia, nadir, scenario, lca_type):
    """
    the main control loop for conducting adaptive weight sums to identify the pareto front
    :param M: the model being used
    :param divs: the initial number of divisions
    :param midpoint: the LCA ReCiPe lca_midpoint used
    :param utopia: the utopia point for the pareto front
    :param nadir: the nadir point for the pareto front
    :param scenario: the scenario number
    :param lca_type: the type of LCA being conducted (ALCA/CLCA)
    :return: the npv values and the lca_midpoint values that comprise the pareto front
    """
    # following Adaptive weighted-sum method for bi-objective optimization:Pareto front generation
    print("\n\n AWS!!!!")
    npv_new = []
    gwp_new = []

    # introduce the scaling factors
    scalar = abs((nadir[0] - utopia[0]) / (nadir[1] - utopia[1]))

    npv_vals, gwp_vals, models = ws(M, divs, lca_type, midpoint, scalar)

    # append the anchor point to the pareto front
    if len(npv_vals) > 0:
        print("new pareto point appended - anchor", npv_vals[0], gwp_vals[0])
        npv_new.append(npv_vals[0])
        gwp_new.append(gwp_vals[0])
        print_model(scenario, models[0], int(pyo.value(models[0].npv[0])), "TEA")
        print_model(scenario, models[0], int(pyo.value(models[0].npv[0])), "LCA", midpoint=midpoint)
    else:
        # return if there are no values from weighted sums, there's no point in further refinement
        return [], []

    # calculate the euclidean distances between the points on the pareto front
    dist = pareto_point_distance(gwp_vals, npv_vals)

    # avoid those damn divide by zero errors
    if len(dist) == 0:
        return [], []

    # refinement along sections of the pareto front
    delta, n = ws_bound_refinement(dist, nadir, utopia)

    # go through the points on the pareto front
    for i in range(len(n)):
        if n[i] > 1.6:
            # determine the offset distances
            try:
                theta = math.atan(-((gwp_vals[i] - gwp_vals[i + 1]) / (npv_vals[i] - npv_vals[i + 1])))
            except ZeroDivisionError:
                # if we divide by zero, there is no need to subdivide this region any further
                return [], []

            delta1 = delta * math.cos(theta)
            delta2 = delta * math.sin(theta)

            print("\ndelta1", delta1)
            print("delta2", delta2)

            # solve the new submodel with additional constraints by calling the aws method
            submodel = M.clone()
            print("left point", npv_vals[i], gwp_vals[i])
            print("right point", npv_vals[i + 1], gwp_vals[i + 1])
            print("new lower npv bound", npv_vals[i] + delta1)
            print("new upper gwp bound", gwp_vals[i + 1] + delta2)
            if scenario < 4000:
                submodel.const.add(expr=sum(submodel.npv[l] for l in submodel.Location) >= npv_vals[i] + delta1)
                submodel.const.add(
                    expr=sum(submodel.total_LCA_midpoints[l, lca_type, midpoint] for l in submodel.Location) * scalar <=
                         gwp_vals[i + 1] + delta2)
            else:
                submodel.const.add(expr=sum(submodel.total_LCA_midpoints[l, lca_type, 'climate change'] for l in submodel.Location) >= npv_vals[i] + delta1)
                submodel.const.add(
                    expr=sum(submodel.total_LCA_midpoints[l, lca_type, 'eutrophication: freshwater'] for l in submodel.Location) * scalar <=
                         gwp_vals[i + 1] + delta2)

            # get results back from the submodel
            x, y = aws(submodel, math.ceil(n[i]), midpoint, utopia, nadir, scenario, lca_type)

            # append the returned values to our lists of points on the pareto front
            for j in range(len(x)):
                npv_new.append(x[j])
                gwp_new.append(y[j])

        # if the distance between this new point and the previous is less than half of delta, don't add the new point
        pf_dist = (((npv_new[len(npv_new) - 1] - npv_vals[i + 1]) ** 2) +
                   ((gwp_new[len(gwp_new) - 1] - gwp_vals[i + 1]) ** 2)) ** (1 / 2)
        print("old point", npv_new[len(npv_new) - 1], gwp_new[len(gwp_new) - 1])
        print("new point", npv_vals[i + 1], gwp_vals[i + 1])
        print("distance along pareto front", pf_dist)
        if pf_dist > delta / 2:
            print("new pareto point appended!")
            npv_new.append(npv_vals[i + 1])
            gwp_new.append(gwp_vals[i + 1])
            print_model(scenario, models[i + 1], int(pyo.value(models[i + 1].npv[0])), "TEA")
            print_model(scenario, models[i + 1], int(pyo.value(models[i + 1].npv[0])), "LCA", midpoint=midpoint)

    # because the pareto front is monotonic, we can sort the points without losing order (is it even necessary to sort?)
    print("pareto front points", npv_new, gwp_new)
    npv_new.sort()
    gwp_new.sort()
    print("pareto front points sorted", npv_new, gwp_new)

    # return from the algorithm
    print("\n\n\n\n\n return from an iteration!!!!!")
    print(npv_new, gwp_new)
    return npv_new, gwp_new


def pareto_point_distance(gwp_vals, npv_vals):
    """
    calculates the distance between the pareto points
    :param gwp_vals: the list of lca_midpoint values
    :param npv_vals: the list of npv values
    :return: the list continaing the distance between the points in the list
    """
    dist = []
    for i in range(len(npv_vals) - 1):
        dist.append((((npv_vals[i] - npv_vals[i + 1]) ** 2) + ((gwp_vals[i] - gwp_vals[i + 1]) ** 2)) ** (1 / 2))
    print("dist", dist)
    return dist


def ws_bound_refinement(dist, nadir, utopia):
    """
    Updates the lower and upper bounds on regions of the aws curve
    :param dist: the list containing the distance between points
    :param nadir: the nadir point for the pareto curve
    :param utopia: the utopia point for the pareto curve
    :return: the offset parameter delta and the number of refinements to conduct in each segment
    """
    avg_length = sum(i for i in dist) / len(dist)
    C = 1.5
    delta = abs(nadir[0] - utopia[0]) / 8
    if delta < 1000:
        delta = 1000
    print("delta", delta)
    n = []
    for i in range(len(dist)):
        # to cut down on the number of iterations
        if avg_length > delta:
            temp = C * dist[i] / avg_length
        else:
            temp = C * dist[i] / delta
        if temp > 6:
            n.append(6)
        else:
            n.append(temp)
    print("n", n)
    return delta, n


def ws(M, divisions, lca_type, midpoint, scalar):
    """
    Conducts weighted sums on a region of a pareto front
    :param M: the model being optimized
    :param divisions: the number of divisions to be made on the pareto front
    :param lca_type: the LCA type (ALCA, CLCA)
    :param midpoint: the ReCiPe lca_midpoint used on the pareto front
    :param scalar: the scalar for the objective function so that the distance function works accurately
    :return: a list of npv values, lca_midpoint values, and copies of the models for points alongside the pareto front
    """
    # lists of objective values
    models = []
    npv_vals = []
    gwp_vals = []

    # do normal weighted sums
    for i in range(divisions + 1):
        alpha = 1 / divisions * i
        print("\n\nalpha", alpha)
        M.alpha = alpha
        model = M
        opt = pyo.SolverFactory('gurobi')
        opt.options['TimeLimit'] = 600
        try:
            results = opt.solve(model, tee=True)

            # check for infeasible solutions - don't add infeasible solutions to the pareto front
            if not results.solver.termination_condition == TerminationCondition.optimal:
                print("infeasible solution reached")
            else:
                npv_vals.append(pyo.value(model.npv[0]))
                gwp_vals.append(pyo.value(model.total_LCA_midpoints[0, lca_type, midpoint]) * scalar)
                models.append(model.clone())
        except ValueError:
            print("model did not find a solution within the time limit")
    print("\n\n\nvalues from weighted sums iteration")
    print("npv values", npv_vals)
    print("gwp values", gwp_vals)
    for i in range(len(models)):
        print("model npv", i, pyo.value(models[i].npv[0]))

    return npv_vals, gwp_vals, models


def pareto_front(M, midpoint, scenario, A, lca_type):
    """
    Calculates the pareto front
    :param M: the model used in the pareto front
    :param midpoint: the lca_midpoint used for the objective function
    :param scenario:  the scenario number
    :param A: parameters from the CLCO_Data
    :param lca_type: the LCA type (ALCA/CLCA)
    :return:
    """
    if scenario < 4000:
        M.Obj = pyo.Objective(expr=sum(M.npv[l] for l in M.Location), sense=pyo.maximize)
        M.Obj2 = pyo.Objective(expr=sum(M.total_LCA_midpoints[l, lca_type, midpoint] for l in M.Location),
                               sense=pyo.minimize)
    else:
        M.Obj = pyo.Objective(expr=sum(M.total_LCA_midpoints[l, lca_type, "climate change"] for l in M.Location),
                           sense=pyo.minimize)
        M.Obj2 = pyo.Objective(expr=sum(M.total_LCA_midpoints[l, lca_type, 'eutrophication: freshwater'] for l in M.Location),
                               sense=pyo.minimize)
    try:
        utopia, nadir = utopian(M, midpoint, lca_type)
    except ValueError as err:
        print(err.args)
        print("no further attempt to find the pareto front")
        return

    print("\n\n finished computing utopia and nadir points")

    # use the new objective function with the new weights
    if scenario < 4000:
        M.combined = pyo.Objective(
            expr=M.alpha * sum(M.npv[l] for l in M.Location) -
                 (1 - M.alpha) * abs((nadir[0] - utopia[0]) / (nadir[1] - utopia[1])) * sum(
                M.total_LCA_midpoints[l, lca_type, midpoint] for l in M.Location), sense=pyo.maximize)
    else:
        M.combined = pyo.Objective(
            expr=M.alpha * sum(M.total_LCA_midpoints[l, lca_type, "climate change"] for l in M.Location) +
                 (1 - M.alpha) * abs((nadir[0] - utopia[0]) / (nadir[1] - utopia[1])) * sum(
                M.total_LCA_midpoints[l, lca_type, 'eutrophication: freshwater'] for l in M.Location), sense=pyo.minimize)

    print("returned to control method")
    # gather the points on the pareto front
    x, y = aws(M, 4, midpoint, utopia, nadir, scenario, lca_type)

    # rescale the y points back to their original values
    scalar = abs((nadir[0] - utopia[0]) / (nadir[1] - utopia[1]))
    y_rescaled = [sum(i / (scalar * A.FEEDSTOCK_SUPPLY[l] * A.TIME_PERIODS) for l in M.Location) for i in y]
    x_rescaled = [sum(i / (A.FEEDSTOCK_SUPPLY[l] * A.TIME_PERIODS) for l in M.Location) for i in x]
    print("rescaled points")

    # plot the points
    plt.clf()
    plt.plot(x_rescaled, y_rescaled, 'ob-')
    plt.title("Pareto Front (Adaptive Weighted Sums)")
    if scenario < 4000:
        plt.xlabel("NPV ($USD) per ton manure")
        plt.ylabel(str(midpoint) + " impact per ton manure")
    else:
        plt.xlabel("climate change impact (kg CO2-eq/ton manure)")
        plt.ylabel("freshwater eutrophication impact (kg P-eq/ton manure)")
    plt.savefig(save_plot(scenario, midpoint=midpoint), dpi=300)
    print("saved fig")
    # plt.show()
    return 1


def initialize_model(scenario, j, midpoint, lca_type):
    """
    Initializes the CLCO model
    :param scenario: the scenario number
    :param j: the county number
    :param midpoint: the lca_midpoint type for MOO
    :param lca_type: the LCA type for MOO
    :return:
    """
    A = CLCO_Data(scenario)

    #### MODEL
    M = pyo.ConcreteModel(scenario)

    add_sets(A, M, j)

    add_variables(M)

    add_constraints(A, M, scenario)

    # solving the model
    if 1000 < scenario < 3000 or scenario in [4501, 4502, 4503, 4511, 4512, 4513]:
        return pareto_front(M, midpoint, scenario, A, lca_type)
    elif 2999 < scenario < 9999:
        if int((scenario / 100) % 10) == 0:
            M.Obj = pyo.Objective(expr=sum(M.npv[l] for l in M.Location), sense=pyo.maximize)
        if int((scenario / 100) % 10) == 1:
            M.Obj = pyo.Objective(
                expr=sum(M.npv[l] - M.total_LCA_midpoints[l, "CLCA", "climate change"] / 5 for l in M.Location),
                sense=pyo.maximize)
        if int((scenario / 100) % 10) == 2:
            M.Obj = pyo.Objective(
                expr=sum(M.npv[l] - M.total_LCA_midpoints[l, "CLCA", "climate change"] for l in M.Location),
                sense=pyo.maximize)

    elif scenario == 50:
        M.Obj = pyo.Objective(expr=sum(M.total_LCA_midpoints[l, lca_type, "climate change"] for l in M.Location),
                              sense=pyo.minimize)
    elif scenario in [425]:
        M.Obj = pyo.Objective(
            expr=sum(M.npv[l] - M.total_LCA_midpoints[l, "CLCA", "climate change"] / 5 for l in M.Location),
            sense=pyo.maximize)
    else:
        M.Obj = pyo.Objective(expr=sum(M.npv[l] for l in M.Location), sense=pyo.maximize)

    model = M
    opt = pyo.SolverFactory('gurobi')

    if scenario in [5]:
        opt.options["mipgap"] = .001

    print(opt.solve(model, tee=True))

    model.npv.pprint()

    # print data to csv
    print_model(scenario, model, j, "TEA")
    print_model(scenario, model, j, "LCA")


def add_constraints(A, M, scenario):
    """
    Adds constraints to the model
    :param A: parameters for the model
    :param M: the model
    :param scenario: the scenario number
    :return:
    """
    # constraint list definition
    M.const = pyo.ConstraintList()

    for l in M.Location:
        # only one stage can be selected for each location for AD
        M.const.add(expr=sum(M.decision_ad_stage[l, stage] for stage in M.ADStages) == 1)

        # time dependent constraints
        for t in M.Time:
            facility_constraints(A, M, l, scenario, t)
            revenue_constraints(A, M, l, scenario, t)
            opex_constraints(A, M, l, t)

        # time independent constraints
        capex_constraints(A, M, l)
        lca_constraints(A, M, l)

        # NPV discounting for each location
        M.const.add(expr=M.npv[l] == - sum(M.process_capex[l, tech] + M.storage_capex[l, tech] for tech in M.Technology)
                         + sum(M.opex_revenues[l, t, tech, revenue_type]
                               for t in M.Time for tech in M.Technology for revenue_type in M.OPEXSubRevenues) -
                         sum(M.opex_costs[l, t, tech, cost_type]
                             for cost_type in M.OPEXSubCosts for t in M.Time for tech in M.Technology))


def lca_constraints(A, M, l):
    """
    Adds constraints on ALCA/CLCA to model
    :param A: data
    :param M: model
    :param l: location
    :return: N/A
    """
    # LCA calculations
    for cat in M.LCAMidpointCat:
        for lca_type in M.LCATypes:
            # mid point for a point in time
            M.const.add(expr=M.LCA_midpoints[l, lca_type, 'natural gas', cat] ==
                             A.IMPACT[lca_type, 'natural gas', cat] * sum(M.purchased_fuel[l, t, tech]
                                                                          for tech in M.Technology for t in M.Time))
            M.const.add(expr=M.LCA_midpoints[l, lca_type, 'grid electricity', cat] ==
                             A.IMPACT[lca_type, 'grid electricity', cat] * sum(M.purchased_power[l, t, tech]
                                                                               for tech in M.Technology for t in
                                                                               M.Time))
            M.const.add(expr=M.LCA_midpoints[l, lca_type, 'diesel', cat] == A.IMPACT[lca_type, 'diesel', cat]
                             * (sum(M.inputs[l, t, tech, 'bio-oil diesel'] for tech in M.Technology for t in M.Time) +
                             sum(M.inputs[l, t, tech, 'transportation'] * A.INTRA_COUNTY_TRANSPORT_DISTANCE[l]
                                 * A.DIESEL_USE for t in M.Time for tech in M.Technology)))
            M.const.add(expr=M.LCA_midpoints[l, lca_type, 'water', cat] == A.IMPACT[lca_type, 'water', cat] *
                             sum(M.inputs[l, t, tech, 'water'] for tech in M.Technology for t in M.Time))
            M.const.add(expr=M.LCA_midpoints[l, lca_type, 'biochar-chp', cat] ==
                             sum(A.IMPACT[lca_type, 'biochar-chp', temp, cat] * M.biochar_from_pyrolysis[
                                 l, t, feed, temp, 'CHP'] for temp in M.PyrolysisTemperatures
                                 for feed in M.PyrolysisFeedstocks for t in M.Time))
            M.const.add(expr=M.LCA_midpoints[l, lca_type, 'biochar-land', cat] ==
                             sum(A.IMPACT[lca_type, 'biochar-land', temp, cat] * M.biochar_from_pyrolysis[
                                 l, t, feed, temp, 'land'] for temp in M.PyrolysisTemperatures
                                 for feed in M.PyrolysisFeedstocks for t in M.Time))
            M.const.add(expr=M.LCA_midpoints[l, lca_type, 'biochar-disposal', cat] ==
                             A.IMPACT[lca_type, 'biochar-disposal', cat] *
                             sum(M.biochar_from_pyrolysis[l, t, feed, temp, 'disposal']
                                 for temp in M.PyrolysisTemperatures for feed in M.PyrolysisFeedstocks for t in M.Time))
            M.const.add(expr=M.LCA_midpoints[l, lca_type, 'pyro-bio-oil-chp', cat] ==
                             sum(A.IMPACT[lca_type, 'pyro-bio-oil-chp', temp, cat] * M.biooil_from_pyrolysis[
                                 l, t, feed, temp, 'CHP'] for temp in M.PyrolysisTemperatures for feed
                                 in M.PyrolysisFeedstocks for t in M.Time))
            M.const.add(expr=M.LCA_midpoints[l, lca_type, 'syngas-chp', cat] ==
                             sum(A.IMPACT[lca_type, 'syngas-chp', temp, cat] * M.syngas_from_pyrolysis[
                                 l, t, feed, temp, 'CHP'] for temp in M.PyrolysisTemperatures
                                 for feed in M.PyrolysisFeedstocks for t in M.Time))
            M.const.add(expr=M.LCA_midpoints[l, lca_type, 'syngas-disposal', cat] ==
                             sum(A.IMPACT[lca_type, 'syngas-disposal', temp, cat] * M.syngas_from_pyrolysis[
                                 l, t, feed, temp, 'disposal'] for temp in M.PyrolysisTemperatures
                                 for feed in M.PyrolysisFeedstocks for t in M.Time))
            M.const.add(expr=M.LCA_midpoints[l, lca_type, 'pyro-ap-disposal', cat] ==
                             A.IMPACT[lca_type, 'pyro-ap-disposal', cat] *
                             sum(M.ap_from_pyrolysis[l, t, feed, temp, 'disposal']
                                 for temp in M.PyrolysisTemperatures for feed in M.PyrolysisFeedstocks for t in M.Time))
            M.const.add(expr=M.LCA_midpoints[l, lca_type, 'htl-hydrochar-land', cat] ==
                             A.IMPACT[lca_type, 'htl-hydrochar-land', cat]
                             * sum(M.hydrochar_from_htl[l, t, feed, temp, 'land']
                                   for temp in M.HTLTemperatures for feed in M.HTLFeedstocks for t in M.Time))
            M.const.add(expr=M.LCA_midpoints[l, lca_type, 'htl-hydrochar-chp', cat] ==
                             A.IMPACT[lca_type, 'htl-hydrochar-chp', cat]
                             * sum(M.hydrochar_from_htl[l, t, feed, temp, 'CHP']
                                   for temp in M.HTLTemperatures for feed in M.HTLFeedstocks for t in M.Time))
            M.const.add(expr=M.LCA_midpoints[l, lca_type, 'htl-hydrochar-disposal', cat] ==
                             A.IMPACT[lca_type, 'htl-hydrochar-disposal', cat] *
                             sum(M.hydrochar_from_htl[l, t, feed, temp, 'disposal'] for temp in M.HTLTemperatures
                                 for feed in M.HTLFeedstocks for t in M.Time))
            M.const.add(expr=M.LCA_midpoints[l, lca_type, 'htl-bio-oil-chp', cat] ==
                             A.IMPACT[lca_type, 'htl-bio-oil-chp', cat] * sum(M.biooil_from_htl[l, t, feed, temp, 'CHP']
                                                                              for temp in M.HTLTemperatures for feed in
                                                                              M.HTLFeedstocks for t in M.Time))
            M.const.add(expr=M.LCA_midpoints[l, lca_type, 'htl-gp-disposal', cat] ==
                             A.IMPACT[lca_type, 'htl-gp-disposal', cat] *
                             sum(M.gp_from_htl[l, t, feed, temp, 'disposal']
                                 for temp in M.HTLTemperatures for feed in M.HTLFeedstocks for t in M.Time))
            M.const.add(expr=M.LCA_midpoints[l, lca_type, 'htl-ap-disposal', cat] ==
                             A.IMPACT[lca_type, 'htl-ap-disposal', cat] *
                             sum(M.ap_from_htl[l, t, feed, temp, 'disposal']
                                 for temp in M.HTLTemperatures for feed in M.HTLFeedstocks for t in M.Time))
            M.const.add(expr=M.LCA_midpoints[l, lca_type, 'htc-hydrochar-land', cat] ==
                             A.IMPACT[lca_type, 'htc-hydrochar-land', cat] * sum(
                M.hydrochar_from_htc[l, t, feed, temp, 'land']
                for temp in M.HTCTemperatures for feed in M.HTCFeedstocks for t in M.Time))
            M.const.add(expr=M.LCA_midpoints[l, lca_type, 'htc-hydrochar-chp', cat] ==
                             sum(A.IMPACT[lca_type, 'htc-hydrochar-chp', temp, cat] * M.hydrochar_from_htc[
                                 l, t, feed, temp, 'CHP']
                                 for temp in M.HTCTemperatures for feed in M.HTCFeedstocks for t in M.Time))
            M.const.add(expr=M.LCA_midpoints[l, lca_type, 'htc-hydrochar-disposal', cat] ==
                             A.IMPACT[lca_type, 'htc-hydrochar-disposal', cat] *
                             sum(M.hydrochar_from_htc[l, t, feed, temp, 'disposal']
                                 for temp in M.HTCTemperatures for feed in M.HTCFeedstocks for t in M.Time))
            M.const.add(expr=M.LCA_midpoints[l, lca_type, 'htc-gp-disposal', cat] ==
                             sum(A.IMPACT[lca_type, 'htc-gp-disposal', temp, cat] *
                                 M.gp_from_htc[l, t, feed, temp, 'disposal']
                                 for temp in M.HTCTemperatures for feed in M.HTCFeedstocks for t in M.Time))
            M.const.add(expr=M.LCA_midpoints[l, lca_type, 'htc-ap-disposal', cat] ==
                             A.IMPACT[lca_type, 'htc-ap-disposal', cat] *
                             sum(M.ap_from_htc[l, t, feed, temp, 'disposal']
                                 for temp in M.HTCTemperatures for feed in M.HTCFeedstocks for t in M.Time))
            M.const.add(expr=M.LCA_midpoints[l, lca_type, 'digestate-land', cat] ==
                             A.IMPACT[lca_type, 'digestate-land', cat] * sum(M.digestate_from_ad[l, t, stage, 'land']
                                                                             for stage in M.ADStages for t in M.Time))
            M.const.add(expr=M.LCA_midpoints[l, lca_type, 'digestate-disposal', cat] ==
                             A.IMPACT[lca_type, 'digestate-disposal', cat] *
                             sum(M.digestate_from_ad[l, t, stage, 'disposal'] for stage in M.ADStages for t in M.Time))
            M.const.add(expr=M.LCA_midpoints[l, lca_type, 'biogas-disposal', cat] ==
                             A.IMPACT[lca_type, 'biogas-disposal', cat] * sum(M.biogas_from_ad[l, t, stage, 'disposal']
                                                                              for stage in M.ADStages for t in M.Time))
            M.const.add(expr=M.LCA_midpoints[l, lca_type, 'biogas-chp', cat] ==
                             A.IMPACT[lca_type, 'biogas-chp', cat] * sum(M.biogas_from_ad[l, t, stage, 'CHP']
                                                                         for stage in M.ADStages for t in M.Time))
            M.const.add(expr=M.LCA_midpoints[l, lca_type, 'manure-land', cat] ==
                             A.IMPACT[lca_type, 'manure-land', cat] *
                             sum(M.feedstock_from_storage[l, t] for t in M.Time))
            M.const.add(expr=M.LCA_midpoints[l, lca_type, 'facility construction', cat] ==
                             A.IMPACT[lca_type, 'facility construction', "ad", cat] * M.process_capacity[l, 'AD'] +
                             A.IMPACT[lca_type, 'facility construction', "chem", cat] * (M.process_capex[l, 'Pyrolysis']
                                                                                         + M.process_capex[l, 'HTL'] +
                                                                                         M.process_capex[l, 'HTC']))
            # CHP facility construction LCA impacts are spread out across the various CHP categories
            M.const.add(expr=M.LCA_midpoints[l, lca_type, 'storage-facility-solids', cat] ==
                             A.IMPACT[lca_type, 'storage-facility-solids', cat] * (M.feedstock_storage_capacity[l] +
                                                                                   M.pyrolysis_storage_capacity[
                                                                                       l, 'Biochar'] +
                                                                                   M.htl_storage_capacity[
                                                                                       l, 'Hydrochar'] +
                                                                                   M.htc_storage_capacity[
                                                                                       l, 'Hydrochar'] +
                                                                                   M.ad_storage_capacity[
                                                                                       l, 'digestate']))
            M.const.add(expr=M.LCA_midpoints[l, lca_type, 'storage-facility-liquids', cat] ==
                             A.IMPACT[lca_type, 'storage-facility-liquids', cat] *
                             (M.pyrolysis_storage_capacity[l, 'AP'] +
                              M.pyrolysis_storage_capacity[l, 'Syngas'] +
                              M.htl_storage_capacity[l, 'AP'] +
                              M.htc_storage_capacity[l, 'AP'] +
                              M.ad_storage_capacity[l, 'biogas']))
            M.const.add(expr=M.LCA_midpoints[l, lca_type, "biochar market", cat] == 0)
            '''A.IMPACT[lca_type,
                'biochar market', cat] *
                             sum(M.biochar_from_pyrolysis[l, t, feedstock, temp, 'market']
                                  for l in M.Location for t in M.Time for feedstock in M.PyrolysisFeedstocks
                                  for temp in M.PyrolysisTemperatures))'''  # for when biochar can be sold on the market - still needs to be implemented
            M.const.add(expr=M.LCA_midpoints[l, lca_type, "bio-oil market", cat] ==
                             A.IMPACT[lca_type, 'bio-oil market', cat] *
                             (sum(M.biooil_from_pyrolysis[l, t, feedstock, temp, 'market']
                                  for l in M.Location for t in M.Time for feedstock in M.PyrolysisFeedstocks
                                  for temp in M.PyrolysisTemperatures)
                              + sum(M.biooil_from_htl[l, t, feedstock, temp, 'market']
                                    for l in M.Location for t in M.Time for feedstock in M.HTLFeedstocks
                                    for temp in M.HTLTemperatures)))
            M.const.add(expr=M.LCA_midpoints[l, lca_type, "hydrochar market", cat] ==
                             A.IMPACT[lca_type, 'hydrochar market', cat] *
                             (sum(M.hydrochar_from_htc[l, t, feedstock, temp, 'market']
                                  for l in M.Location for t in M.Time for feedstock in M.HTCFeedstocks
                                  for temp in M.HTCTemperatures)
                              + sum(M.hydrochar_from_htl[l, t, feedstock, temp, 'market']
                                    for l in M.Location for t in M.Time for feedstock in M.HTLFeedstocks
                                    for temp in M.HTLTemperatures)))

            M.const.add(expr=M.total_LCA_midpoints[l, lca_type, cat] ==
                             sum(M.LCA_midpoints[l, lca_type, origin, cat] for origin in M.ALCAInputs))

            # LCA specific constraints
            if lca_type == "ALCA":
                M.const.add(expr=M.LCA_midpoints[l, lca_type, "avoided electricity", cat] == 0)
                M.const.add(expr=M.LCA_midpoints[l, lca_type, 'N fertilizer', cat] == 0)
                M.const.add(expr=M.LCA_midpoints[l, lca_type, 'P fertilizer', cat] == 0)
                M.const.add(expr=M.LCA_midpoints[l, lca_type, 'K fertilizer', cat] == 0)
            elif lca_type == "CLCA":
                M.const.add(expr=M.LCA_midpoints[l, lca_type, "avoided electricity", cat] ==
                                 -A.IMPACT[lca_type, 'grid electricity', cat] *
                                 sum(M.chp_market[l, t, 'electricity'] for t in M.Time))
                M.const.add(expr=M.LCA_midpoints[l, lca_type, 'N fertilizer', cat] ==
                                 -A.IMPACT[lca_type, 'N fertilizer', cat] *
                                 sum(M.avoided_fertilizers[l, t, tech, 'N']
                                     for l in M.Location for t in M.Time for tech in M.Technology))
                M.const.add(expr=M.LCA_midpoints[l, lca_type, 'P fertilizer', cat] ==
                                 -A.IMPACT[lca_type, 'P fertilizer', cat] *
                                 sum(M.avoided_fertilizers[l, t, tech, 'P']
                                     for l in M.Location for t in M.Time for tech in M.Technology))
                M.const.add(expr=M.LCA_midpoints[l, lca_type, 'K fertilizer', cat] ==
                                 -A.IMPACT[lca_type, 'K fertilizer', cat] *
                                 sum(M.avoided_fertilizers[l, t, tech, 'K']
                                     for l in M.Location for t in M.Time for tech in M.Technology))


def capex_constraints(A, M, l):
    """
    Constraints for calculating the CAPEX for the model
    :param A: parameter information
    :param M: the model
    :param l: the location
    :return: N/A
    """
    # CAPEX for storage and process capacity
    # implementing piecewise linear approximation
    q = symbols("q")
    lpa_xvals = []
    [lpa_xvals.append(x) for x in A.ORIGINAL_FEEDSTOCK_SUPPLY if x not in lpa_xvals]
    lpa_xvals.append(0.00000001)
    for i in range(16):
        lpa_xvals.append(i * 1000 + 4500)
    for i in range(19):
        lpa_xvals.append(i * 10000 + 20000)
    lpa_xvals.sort()
    pyro_process = A.CAPEX['Pyrolysis', 'process', 'coefficient'] * (q * 1000 / A.HOURS_PER_PERIOD) ** A.CAPEX[
        'Pyrolysis', 'process', 'exponent']
    pyro_storage = A.CAPEX['Pyrolysis', 'storage', 'coefficient'] * (q / A.HOURS_PER_PERIOD) ** A.CAPEX[
        'Pyrolysis', 'storage', 'exponent']
    htl_process = A.CAPEX['HTL', 'process', 'coefficient'] * (q / A.DRY_BIOMASS_REF) ** A.CAPEX[
        'HTL', 'process', 'exponent']
    htl_storage = A.CAPEX['HTL', 'storage', 'coefficient'] * (q / A.HOURS_PER_PERIOD) ** A.CAPEX[
        'HTL', 'storage', 'exponent']
    htc_process = A.CAPEX['HTC', 'process', 'coefficient'] * (q / A.DRY_BIOMASS_REF) ** A.CAPEX[
        'HTC', 'process', 'exponent']
    htc_storage = A.CAPEX['HTC', 'storage', 'coefficient'] * (q / A.HOURS_PER_PERIOD) ** A.CAPEX[
        'HTC', 'storage', 'exponent']
    ad_process = A.CAPEX['AD', 'process', 'coefficient'] * q ** A.CAPEX['AD', 'process', 'exponent']
    ad_storage = A.CAPEX['AD', 'storage', 'coefficient'] * q ** A.CAPEX['AD', 'storage', 'exponent']
    chp_process = A.CAPEX['CHP', 'process', 'coefficient'] * q ** A.CAPEX['CHP', 'process', 'exponent']
    solid_storage = A.CAPEX['Solid', 'storage', 'coefficient'] * q ** A.CAPEX['Solid', 'storage', 'exponent']
    pyro_proc_lpa = [pyro_process.evalf(subs={q: x}) for x in lpa_xvals]
    pyro_stor_lpa = [pyro_storage.evalf(subs={q: x}) for x in lpa_xvals]
    solids_products = [solid_storage.evalf(subs={q: x}) for x in lpa_xvals]
    htl_proc_lpa = [htl_process.evalf(subs={q: x}) for x in lpa_xvals]
    htl_stor_lpa = [htl_storage.evalf(subs={q: x}) for x in lpa_xvals]
    htc_proc_lpa = [htc_process.evalf(subs={q: x}) for x in lpa_xvals]
    htc_stor_lpa = [htc_storage.evalf(subs={q: x}) for x in lpa_xvals]
    ad_proc_lpa = [ad_process.evalf(subs={q: x}) for x in lpa_xvals]
    ad_stor_lpa = [ad_storage.evalf(subs={q: x}) for x in lpa_xvals]
    chp_proc_lpa = [chp_process.evalf(subs={q: x}) for x in lpa_xvals]

    # Piecewise linear approximation to capital cost terms
    for tech in M.Technology:
        if tech == 'Pyrolysis':
            M.con5 = pyo.Piecewise(M.process_capex[l, 'Pyrolysis'], M.process_capacity[l, 'Pyrolysis'],
                                   pw_pts=lpa_xvals,
                                   pw_constr_type='EQ',
                                   f_rule=pyro_proc_lpa,
                                   pw_repn='SOS2')
            M.con6 = pyo.Piecewise(M.pyro_storage_cost[l, 'Biochar'], M.pyrolysis_storage_capacity[l, 'Biochar'],
                                   pw_pts=lpa_xvals,
                                   pw_constr_type='EQ',
                                   f_rule=solids_products,
                                   pw_repn='SOS2')
            M.con7 = pyo.Piecewise(M.pyro_storage_cost[l, 'Biooil'], M.pyrolysis_storage_capacity[l, 'Biooil'],
                                   pw_pts=lpa_xvals,
                                   pw_constr_type='EQ',
                                   f_rule=pyro_stor_lpa,
                                   pw_repn='SOS2')
            M.con8 = pyo.Piecewise(M.pyro_storage_cost[l, 'AP'], M.pyrolysis_storage_capacity[l, 'AP'],
                                   pw_pts=lpa_xvals,
                                   pw_constr_type='EQ',
                                   f_rule=pyro_stor_lpa,
                                   pw_repn='SOS2')
            M.con9 = pyo.Piecewise(M.pyro_storage_cost[l, 'Syngas'], M.pyrolysis_storage_capacity[l, 'Syngas'],
                                   pw_pts=lpa_xvals,
                                   pw_constr_type='EQ',
                                   f_rule=pyro_stor_lpa,
                                   pw_repn='SOS2')
            M.const.add(expr=M.storage_capex[l, 'Pyrolysis'] == sum(
                M.pyro_storage_cost[l, prod] for prod in M.PyrolysisProducts))

        elif tech == 'HTL':
            M.con10 = pyo.Piecewise(M.process_capex[l, 'HTL'], M.process_capacity[l, 'HTL'],
                                    pw_pts=lpa_xvals,
                                    pw_constr_type='EQ',
                                    f_rule=htl_proc_lpa,
                                    pw_repn='SOS2')
            M.con11 = pyo.Piecewise(M.htl_storage_cost[l, 'Hydrochar'], M.htl_storage_capacity[l, 'Hydrochar'],
                                    pw_pts=lpa_xvals,
                                    pw_constr_type='EQ',
                                    f_rule=solids_products,
                                    pw_repn='SOS2')
            M.con12 = pyo.Piecewise(M.htl_storage_cost[l, 'Biooil'], M.htl_storage_capacity[l, 'Biooil'],
                                    pw_pts=lpa_xvals,
                                    pw_constr_type='EQ',
                                    f_rule=htl_stor_lpa,
                                    pw_repn='SOS2')
            M.con13 = pyo.Piecewise(M.htl_storage_cost[l, 'AP'], M.htl_storage_capacity[l, 'AP'],
                                    pw_pts=lpa_xvals,
                                    pw_constr_type='EQ',
                                    f_rule=htl_stor_lpa,
                                    pw_repn='SOS2')
            M.con14 = pyo.Piecewise(M.htl_storage_cost[l, 'GP'], M.htl_storage_capacity[l, 'GP'],
                                    pw_pts=lpa_xvals,
                                    pw_constr_type='EQ',
                                    f_rule=htl_stor_lpa,
                                    pw_repn='SOS2')
            M.const.add(expr=M.storage_capex[l, 'HTL'] == sum(M.htl_storage_cost[l, prod] for prod in M.HTLProducts))
        elif tech == 'HTC':
            M.con15 = pyo.Piecewise(M.process_capex[l, 'HTC'], M.process_capacity[l, 'HTC'],
                                    pw_pts=lpa_xvals,
                                    pw_constr_type='EQ',
                                    f_rule=htc_proc_lpa,
                                    pw_repn='SOS2')
            M.con16 = pyo.Piecewise(M.htc_storage_cost[l, 'Hydrochar'], M.htc_storage_capacity[l, 'Hydrochar'],
                                    pw_pts=lpa_xvals,
                                    pw_constr_type='EQ',
                                    f_rule=solids_products,
                                    pw_repn='SOS2')
            M.con17 = pyo.Piecewise(M.htc_storage_cost[l, 'AP'], M.htc_storage_capacity[l, 'AP'],
                                    pw_pts=lpa_xvals,
                                    pw_constr_type='EQ',
                                    f_rule=htc_stor_lpa,
                                    pw_repn='SOS2')
            M.con18 = pyo.Piecewise(M.htc_storage_cost[l, 'GP'], M.htc_storage_capacity[l, 'GP'],
                                    pw_pts=lpa_xvals,
                                    pw_constr_type='EQ',
                                    f_rule=htc_stor_lpa,
                                    pw_repn='SOS2')
            M.const.add(expr=M.storage_capex[l, 'HTC'] == sum(M.htc_storage_cost[l, prod] for prod in M.HTCProducts))
        elif tech == 'AD':
            M.con1 = pyo.Piecewise(M.process_capex[l, 'AD'], M.process_capacity[l, 'AD'],
                                   pw_pts=lpa_xvals,
                                   pw_constr_type='EQ',
                                   f_rule=ad_proc_lpa,
                                   pw_repn='SOS2')
            M.con2 = pyo.Piecewise(M.ad_storage_cost[l, 'biogas'], M.ad_storage_capacity[l, 'biogas'],
                                   pw_pts=lpa_xvals,
                                   pw_constr_type='EQ',
                                   f_rule=ad_stor_lpa,
                                   pw_repn='SOS2')
            M.con4 = pyo.Piecewise(M.ad_storage_cost[l, 'digestate'], M.ad_storage_capacity[l, 'digestate'],
                                   pw_pts=lpa_xvals,
                                   pw_constr_type='EQ',
                                   f_rule=solids_products,
                                   pw_repn='SOS2')

            M.const.add(expr=M.storage_capex[l, 'AD'] == sum(M.ad_storage_cost[l, prod] for prod in M.ADProducts))

        elif tech == 'CHP':
            M.con3 = pyo.Piecewise(M.process_capex[l, 'CHP'], M.process_capacity[l, 'CHP'],
                                   pw_pts=lpa_xvals,
                                   pw_constr_type='EQ',
                                   f_rule=chp_proc_lpa,
                                   pw_repn='SOS2')  # process capacity of CHP is in kWh
            M.const.add(expr=M.storage_capex[l, 'CHP'] == 0)

        elif tech == 'Feedstock':
            M.const.add(expr=M.process_capex[l, 'Feedstock'] == A.CAPEX['Feedstock', 'process', 'coefficient'] *
                             M.process_capacity[l, 'Feedstock'])
            M.const.add(expr=M.storage_capex[l, 'Feedstock'] == A.CAPEX['Feedstock', 'storage', 'coefficient'] *
                             M.feedstock_storage_capacity[l])


def opex_constraints(A, M, l, t):
    """
    Constraints on the OPEX costs for the model
    :param A: data
    :param M: model
    :param l: location
    :param t: time
    :return: N/A
    """
    for tech in M.Technology:
        # heat inputs for each technology
        if tech == "AD":
            M.const.add(expr=M.inputs[l, t, 'AD', 'heat'] == sum(A.OPEX['AD', 'Heat'] * M.ad_capacity[l, stage] for stage in M.ADStages))
        elif tech == "CHP":
            M.const.add(expr=M.inputs[l, t, 'CHP', 'heat'] == A.OPEX['CHP', 'Heat'] *
                             A.CHP_HEAT_EFFICIENCY * M.chp_in[l, t])
        elif tech == "Feedstock":
            M.const.add(expr=M.inputs[l, t, 'Feedstock', 'heat'] == 0)
        elif tech == "Pyrolysis":
            M.const.add(expr=M.inputs[l, t, 'Pyrolysis', 'heat'] == sum(
                A.OPEX['Pyrolysis', 'Heat', temp] * M.pyrolysis_out[l, t, feed, pyro_prod, temp]
                for temp in M.PyrolysisTemperatures for feed in M.PyrolysisFeedstocks
                for pyro_prod in M.PyrolysisProducts))
        elif tech == "HTL":
            M.const.add(expr=M.inputs[l, t, 'HTL', 'heat'] ==
                             sum(M.htl_out[l, t, feed, prod, temp] * A.OPEX['HTL', 'Heat']
                                 for feed in M.HTLFeedstocks for temp in M.HTLTemperatures for prod in M.HTLProducts))
        elif tech == "HTC":
            M.const.add(expr=M.inputs[l, t, 'HTC', 'heat'] == sum(
                A.OPEX['HTC', 'Heat', temp] * M.htc_out[l, t, feed, prod, temp]
                for feed in M.HTCFeedstocks for temp in M.HTCTemperatures for prod in M.HTCProducts))

        # electricity inputs for each technology
        if tech == "AD":
            M.const.add(expr=M.inputs[l, t, 'AD', 'electricity'] == sum(A.OPEX['AD', 'Electricity'] * M.ad_capacity[l, stage] for stage in M.ADStages))
        elif tech == "CHP":
            M.const.add(expr=M.inputs[l, t, 'CHP', 'electricity'] == A.OPEX['CHP', 'Electricity'] *
                             A.CHP_ELECTRICITY_EFFICIENCY * M.chp_in[l, t])
        elif tech == "Pyrolysis":
            M.const.add(expr=M.inputs[l, t, 'Pyrolysis', 'electricity'] == sum(
                A.OPEX['Pyrolysis', 'Electricity'] * M.pyrolysis_out[l, t, feed, pyro_prod, temp] for feed in
                M.PyrolysisFeedstocks for temp in M.PyrolysisTemperatures for pyro_prod in M.PyrolysisProducts))
        elif tech == "HTL":
            M.const.add(
                expr=M.inputs[l, t, 'HTL', 'electricity'] ==
                     sum(A.OPEX['HTL', 'Electricity'] * M.htl_out[l, t, feed, prod, temp] for feed in
                         M.HTLFeedstocks for temp in M.HTLTemperatures for prod in M.HTLProducts))
        elif tech == "HTC":
            M.const.add(expr=M.inputs[l, t, 'HTC', 'electricity'] == sum(A.OPEX['HTC', 'Electricity'] *
                                                                         M.htc_out[l, t, feed, prod, temp]
                                                                         for feed in M.HTCFeedstocks for temp in
                                                                         M.HTCTemperatures for prod in M.HTCProducts))
        else:
            M.const.add(expr=M.inputs[l, t, tech, 'electricity'] == 0)

        # water inputs
        if tech == "HTC":
            M.const.add(expr=M.inputs[l, t, 'HTC', 'water'] == A.HTC_WATER * M.htc_in[l, t, 'feedstock'])
        elif tech == "HTL":
            M.const.add(expr=M.inputs[l, t, 'HTL', 'water'] == A.HTL_WATER * M.htl_in[l, t, 'feedstock'])
        else:
            M.const.add(expr=M.inputs[l, t, tech, 'water'] == 0)

        # diesel inputs for biodiesel blending
        if tech == "HTL":
            M.const.add(expr=M.inputs[l, t, 'HTL', 'bio-oil diesel'] == A.DIESEL_BIOOIL_RATIO * A.TON_DIESEL_TO_GAL *
                             sum(M.biooil_from_htl[l, t, 'feedstock', temp, 'CHP'] for temp in M.HTLTemperatures))
        elif tech == "Pyrolysis":
            M.const.add(
                expr=M.inputs[l, t, 'Pyrolysis', 'bio-oil diesel'] == A.DIESEL_BIOOIL_RATIO * A.TON_DIESEL_TO_GAL *
                     sum(M.biooil_from_pyrolysis[l, t, 'feedstock', temp, 'CHP'] for temp in M.PyrolysisTemperatures))
        else:
            M.const.add(expr=M.inputs[l, t, tech, 'bio-oil diesel'] == 0)

        # transportation inputs ton
        if tech == "Feedstock":
            M.const.add(
                expr=M.inputs[l, t, tech, "transportation"] == M.feedstock_to_storage[l, t] + M.feedstock_from_storage[
                    l, t])
        elif tech == "Pyrolysis":
            M.const.add(
                expr=M.inputs[l, t, tech, "transportation"] == (M.pyrolysis_in[l, t, 'feedstock'] + sum(
                    M.biochar_from_pyrolysis[l, t, feed, temp, 'land'] for feed in
                    M.PyrolysisFeedstocks for temp in M.PyrolysisTemperatures)))
        elif tech == "HTL":
            M.const.add(expr=M.inputs[l, t, tech, "transportation"] == (M.htl_in[l, t, 'feedstock'] + sum(
                M.hydrochar_from_htl[l, t, feed, temp, 'land'] for feed in M.HTLFeedstocks for
                temp in M.HTLTemperatures)))
        elif tech == "HTC":
            M.const.add(expr=M.inputs[l, t, tech, "transportation"] == (
                    M.htc_in[l, t, 'feedstock'] + sum(M.hydrochar_from_htc[l, t, feed, temp, 'land']
                                                      for feed in M.HTCFeedstocks for temp in
                                                      M.HTCTemperatures)))
        elif tech == "AD":
            M.const.add(expr=M.inputs[l, t, tech, "transportation"] == M.ad_in[l, t, 'feedstock'] +
                             sum(M.digestate_from_ad[l, t, temp, 'land'] for temp in M.ADStages))
        else:
            M.const.add(expr=M.inputs[l, t, tech, "transportation"] == 0)

        # inputs to heat and electricity from the CHP plant
        M.const.add(
            expr=M.inputs[l, t, tech, 'heat'] == M.chp_out[l, t, tech, 'heat'] + A.CHP_HEAT_EFFICIENCY *
                 M.purchased_fuel[l, t, tech])  # units of MJ
        M.const.add(expr=M.inputs[l, t, tech, 'electricity'] == M.chp_out[l, t, tech, 'electricity'] +
                         M.purchased_power[l, t, tech])  # units of kWh
        # OPEX COSTS
        # water costs
        M.const.add(
            expr=M.opex_costs[l, t, tech, 'water'] == A.OPEX['Freshwater'] * M.inputs[l, t, tech, 'water']
                 / ((1 + A.MONTHLY_DISCOUNT_RATE) ** int(t)))
        # heat costs
        M.const.add(
            expr=M.opex_costs[l, t, tech, 'heat'] == (A.OPEX['Fuel'] * M.purchased_fuel[l, t, tech]) / (
                    (1 + A.MONTHLY_DISCOUNT_RATE) ** int(t)))
        # electricity costs
        M.const.add(
            expr=M.opex_costs[l, t, tech, 'electricity'] == (A.OPEX['Electricity'] * M.purchased_power[
                l, t, tech]) / ((1 + A.MONTHLY_DISCOUNT_RATE) ** int(t)))

        # plant labor costs
        if tech == "AD" or tech == "CHP" or tech == "Feedstock":
            M.const.add(expr=M.opex_costs[l, t, tech, 'labor'] == 0 / (
                    (1 + A.MONTHLY_DISCOUNT_RATE) ** int(t)))  # are already incorporated in operational costs
        else:
            M.const.add(expr=M.opex_costs[l, t, tech, 'labor'] ==
                             (A.OPEX['Labor Cost'] * (M.process_capacity[l, tech] /
                                                      A.HOURS_PER_PERIOD ** A.OPEX['Labor Exponent'])) / (
                                     (1 + A.MONTHLY_DISCOUNT_RATE) ** int(t)))
        # transportation costs
        M.const.add(expr=M.opex_costs[l, t, tech, 'transportation'] ==
                         (A.LOAD_TRANSIT_COST + A.OPEX['transit'] * A.INTRA_COUNTY_TRANSPORT_DISTANCE[l]) *
                         M.inputs[l, t, tech, "transportation"] /
                         (1 + A.MONTHLY_DISCOUNT_RATE) ** int(t))

        # Disposal operational expenses depend on the amount of feedstock being disposed of and the method of disposal
        if tech == "AD":
            M.const.add(expr=M.opex_costs[l, t, tech, 'disposal'] ==
                             (A.OPEX['Landfill'] * sum(
                                 M.digestate_from_ad[l, t, stage, 'disposal'] for stage in M.ADStages) +
                              sum(A.OPEX['Atmosphere'] * M.biogas_from_ad[l, t, stage, 'disposal'] for stage in
                                  M.ADStages)) /
                             ((1 + A.MONTHLY_DISCOUNT_RATE) ** int(t)))
        elif tech == "Pyrolysis":
            M.const.add(expr=M.opex_costs[l, t, tech, 'disposal'] ==
                             (A.OPEX['Landfill'] * sum(M.biochar_from_pyrolysis[l, t, feed, temp, 'disposal']
                                                       for feed in M.PyrolysisFeedstocks for temp in
                                                       M.PyrolysisTemperatures) +
                              A.OPEX['Wastewater'] * sum(M.ap_from_pyrolysis[l, t, feed, temp, 'disposal']
                                                         for feed in M.PyrolysisFeedstocks for temp in
                                                         M.PyrolysisTemperatures) +
                              A.OPEX['Atmosphere'] * sum(M.syngas_from_pyrolysis[l, t, feed, temp, 'disposal']
                                                         for feed in M.PyrolysisFeedstocks for temp in
                                                         M.PyrolysisTemperatures)) /
                             ((1 + A.MONTHLY_DISCOUNT_RATE) ** int(t)))
        elif tech == "HTL":
            M.const.add(expr=M.opex_costs[l, t, 'HTL', 'disposal'] ==
                             (A.OPEX['Landfill'] * sum(M.hydrochar_from_htl[l, t, feed, temp, 'disposal']
                                                       for feed in M.HTLFeedstocks for temp in M.HTLTemperatures) +
                              A.OPEX['Wastewater'] * sum(M.ap_from_htl[l, t, feed, temp, 'disposal']
                                                         for feed in M.HTLFeedstocks for temp in M.HTLTemperatures) +
                              A.OPEX['Atmosphere'] * sum(M.gp_from_htl[l, t, feed, temp, 'disposal']
                                                         for feed in M.HTLFeedstocks for temp in M.HTLTemperatures)) /
                             ((1 + A.MONTHLY_DISCOUNT_RATE) ** int(t)))
        elif tech == "HTC":
            M.const.add(expr=M.opex_costs[l, t, 'HTC', 'disposal'] ==
                             (A.OPEX['Landfill'] * sum(M.hydrochar_from_htc[l, t, feed, temp, 'disposal']
                                                       for feed in M.HTCFeedstocks for temp in M.HTCTemperatures) +
                              A.OPEX['Wastewater'] * sum(M.ap_from_htc[l, t, feed, temp, 'disposal']
                                                         for feed in M.HTCFeedstocks for temp in M.HTCTemperatures) +
                              A.OPEX['Atmosphere'] * sum(M.gp_from_htc[l, t, feed, temp, 'disposal']
                                                         for feed in M.HTCFeedstocks for temp in M.HTCTemperatures)) /
                             ((1 + A.MONTHLY_DISCOUNT_RATE) ** int(t)))
        else:  # by definition nothing from chp and feedstock can be disposed
            M.const.add(expr=M.opex_costs[l, t, tech, 'disposal'] == 0 / ((1 + A.MONTHLY_DISCOUNT_RATE) ** int(t)))

        # diesel costs
        M.const.add(
            expr=M.opex_costs[l, t, tech, 'diesel'] == A.DIESEL_PRICE * M.inputs[l, t, tech, 'bio-oil diesel'] / (
                    (1 + A.MONTHLY_DISCOUNT_RATE) ** int(t)))

        # opex plant costs
        M.const.add(expr=M.opex_costs[l, t, tech, 'TPC'] == (A.OPEX_TPC * M.process_capex[l, tech]) / (
                (1 + A.MONTHLY_DISCOUNT_RATE) ** int(t)))


def revenue_constraints(A, M, l, scenario, t):
    """
    Adds revenue constraints to the model
    :param A: data
    :param M: model
    :param l: location
    :param scenario: scenario
    :param t: time
    :return: N/A
    """
    # calculating the amount of avoided fertilizers
    for fertilizer in M.AvoidedFertilizers:
        M.const.add(expr=M.avoided_fertilizers[l, t, 'Pyrolysis', fertilizer] == sum(
            A.NUTRIENTS['Pyrolysis', feed, temp, fertilizer] * M.biochar_from_pyrolysis[
                l, t, feed, temp, 'land'] for feed in M.PyrolysisFeedstocks for temp in
            M.PyrolysisTemperatures))
        M.const.add(expr=M.avoided_fertilizers[l, t, 'HTL', fertilizer] == sum(
            A.NUTRIENTS['HTL', feed, temp, fertilizer] * M.hydrochar_from_htl[l, t, feed, temp, 'land']
            for feed in M.HTLFeedstocks for temp in M.HTLTemperatures))
        M.const.add(expr=M.avoided_fertilizers[l, t, 'HTC', fertilizer] == sum(
            A.NUTRIENTS['HTC', feed, temp, fertilizer] * M.hydrochar_from_htc[l, t, feed, temp, 'land']
            for feed in M.HTCFeedstocks for temp in M.HTCTemperatures))
        M.const.add(expr=M.avoided_fertilizers[l, t, 'AD', fertilizer] == sum(
            A.NUTRIENTS['AD', fertilizer] * M.digestate_from_ad[l, t, temp, 'land']
            for temp in M.ADStages))
        M.const.add(expr=M.avoided_fertilizers[l, t, 'CHP', fertilizer] == 0)
        M.const.add(expr=M.avoided_fertilizers[l, t, 'Feedstock', fertilizer] == A.NUTRIENTS[
            'Feedstock', fertilizer] * M.feedstock_from_storage[l, t])

    # feedstock can only be applied to land in the dedicated month of the year
    if not int(t) % (int(A.TIME_PERIODS / A.YEARS)) == int(A.LAND_APPLICATION_MONTH):
        # no avoided fertilizers allowed when it is not application month
        for fert in M.AvoidedFertilizers:
            for tech in M.Technology:
                M.const.add(expr=M.avoided_fertilizers[l, t, tech, fert] == 0)

    # Selling products on the markets
    # adding the avoided fertilizer and coal revenues to opex
    for tech in M.Technology:
        M.const.add(expr=M.opex_revenues[l, t, tech, 'avoided fertilizer'] == sum(
            A.REVENUE[fert] * M.avoided_fertilizers[l, t, tech, fert] for fert in M.AvoidedFertilizers) / (
                                 (1 + A.MONTHLY_DISCOUNT_RATE) ** int(t)))

        # selling hydrochar as avoided coal, bio oil, and electricity on the markets.
        if tech == "Pyrolysis":
            M.const.add(expr=M.opex_revenues[l, t, tech, 'avoided coal'] == 0 / (
                    (1 + A.MONTHLY_DISCOUNT_RATE) ** int(t)))
            M.const.add(expr=M.opex_revenues[l, t, tech, 'bio oil'] == (A.REVENUE['Biooil'] *
                                                                        sum(A.HHV[
                                                                                'Pyrolysis', feedstock, 'Biooil', temp] *
                                                                            M.biooil_from_pyrolysis[
                                                                                l, t, feedstock, temp, 'market']
                                                                            for feedstock in
                                                                            M.PyrolysisFeedstocks for temp in
                                                                            M.PyrolysisTemperatures)) / (
                                     (1 + A.MONTHLY_DISCOUNT_RATE) ** int(t)))
            M.const.add(expr=M.opex_revenues[l, t, tech, 'electricity'] == 0 / (
                    (1 + A.MONTHLY_DISCOUNT_RATE) ** int(t)))
        elif tech == "HTL":
            M.const.add(expr=M.opex_revenues[l, t, tech, 'avoided coal'] == (A.REVENUE['hydrochar'] *
                                                                             sum(A.HHV[
                                                                                     tech, feedstock, 'Hydrochar', temp] *
                                                                                 M.hydrochar_from_htl[
                                                                                     l, t, feedstock, temp, 'market']
                                                                                 for feedstock in
                                                                                 M.HTLFeedstocks for temp in
                                                                                 M.HTLTemperatures) /
                                                                             ((
                                                                                      1 + A.MONTHLY_DISCOUNT_RATE) ** int(
                                                                                 t))))
            M.const.add(expr=M.opex_revenues[l, t, tech, 'bio oil'] == (A.REVENUE['Biooil'] *
                                                                        sum(A.HHV[
                                                                                'HTL', feedstock, 'Biooil', temp] *
                                                                            M.biooil_from_htl[
                                                                                l, t, feedstock, temp, 'market']
                                                                            for feedstock in M.HTLFeedstocks for
                                                                            temp in M.HTLTemperatures)) / (
                                     (1 + A.MONTHLY_DISCOUNT_RATE) ** int(t)))
            M.const.add(expr=M.opex_revenues[l, t, tech, 'electricity'] == 0 / (
                    (1 + A.MONTHLY_DISCOUNT_RATE) ** int(t)))
        elif tech == "HTC":
            M.const.add(expr=M.opex_revenues[l, t, tech, 'avoided coal'] == (A.REVENUE['hydrochar'] *
                                                                             sum(A.HHV[
                                                                                     tech, feedstock, 'Hydrochar', temp] *
                                                                                 M.hydrochar_from_htc[
                                                                                     l, t, feedstock, temp, 'market']
                                                                                 for feedstock in
                                                                                 M.HTCFeedstocks for temp in
                                                                                 M.HTCTemperatures) /
                                                                             ((
                                                                                      1 + A.MONTHLY_DISCOUNT_RATE) ** int(
                                                                                 t))))
            M.const.add(
                expr=M.opex_revenues[l, t, tech, 'bio oil'] == 0 / ((1 + A.MONTHLY_DISCOUNT_RATE) ** int(t)))
            M.const.add(expr=M.opex_revenues[l, t, tech, 'electricity'] == 0 / (
                    (1 + A.MONTHLY_DISCOUNT_RATE) ** int(t)))
        elif tech == 'CHP':
            M.const.add(expr=M.opex_revenues[l, t, tech, 'avoided coal'] == 0 / (
                    (1 + A.MONTHLY_DISCOUNT_RATE) ** int(t)))
            M.const.add(
                expr=M.opex_revenues[l, t, tech, 'bio oil'] == 0 / ((1 + A.MONTHLY_DISCOUNT_RATE) ** int(t)))
            M.const.add(
                expr=M.opex_revenues[l, t, tech, 'electricity'] == (A.REVENUE['electricity'] * M.chp_market[
                    l, t, 'electricity']) / ((1 + A.MONTHLY_DISCOUNT_RATE) ** int(t)))
        elif tech == 'AD':
            M.const.add(expr=M.opex_revenues[l, t, tech, 'avoided coal'] == 0 / (
                    (1 + A.MONTHLY_DISCOUNT_RATE) ** int(t)))
            M.const.add(
                expr=M.opex_revenues[l, t, tech, 'bio oil'] == 0 / ((1 + A.MONTHLY_DISCOUNT_RATE) ** int(t)))
            M.const.add(expr=M.opex_revenues[l, t, tech, 'electricity'] == 0 / (
                    (1 + A.MONTHLY_DISCOUNT_RATE) ** int(t)))
        else:
            M.const.add(expr=M.opex_revenues[l, t, tech, 'avoided coal'] == 0 / (
                    (1 + A.MONTHLY_DISCOUNT_RATE) ** int(t)))
            M.const.add(
                expr=M.opex_revenues[l, t, tech, 'bio oil'] == 0 / ((1 + A.MONTHLY_DISCOUNT_RATE) ** int(t)))
            M.const.add(expr=M.opex_revenues[l, t, tech, 'electricity'] == 0 / (
                    (1 + A.MONTHLY_DISCOUNT_RATE) ** int(t)))

        # other revenue sources
        if scenario in [5, 10103]:
            if tech == "Feedstock":
                M.const.add(expr=M.opex_revenues[l, t, tech, "incentive 1"] == 0)
                # carbon credits can go here
                # M.const.add(expr=M.opex_revenues[l, t, tech, "incentive 1"] == -100/1000*sum(M.total_LCA_midpoints[l, "ALCA", "climate change"]/A.TIME_PERIODS for l in M.Location))
            else:
                M.const.add(expr=M.opex_revenues[l, t, tech, "incentive 1"] == 0)
        else:
            M.const.add(expr=M.opex_revenues[l, t, tech, "incentive 1"] == 0)
        M.const.add(expr=M.opex_revenues[l, t, tech, "incentive 2"] == 0)
        M.const.add(expr=M.opex_revenues[l, t, tech, "potting media"] == 0)


def facility_constraints(A, M, l, scenario, t):
    """
    Adds facility constraints to the model
    :param A: data
    :param M: model
    :param l: location
    :param scenario: scenario
    :param t: time
    :return: N/A
    """
    # FEEDSTOCK distribution
    if scenario in [1, 421, 431, 441]:
        M.const.add(expr=M.ad_in[l, t, 'feedstock'] == A.FEEDSTOCK_SUPPLY[l])
        M.const.add(expr=M.htc_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.pyrolysis_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.htl_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.feedstock_to_storage[l, t] == 0)
    elif scenario in [2, 422, 432, 442, 10003] or (1000 < scenario < 1019) or (2000 < scenario < 2019):
        M.const.add(expr=M.feedstock_to_storage[l, t] == A.FEEDSTOCK_SUPPLY[l])
        M.const.add(expr=M.ad_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.htc_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.pyrolysis_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.htl_in[l, t, 'feedstock'] == 0)
    elif scenario in [5, 425, 435, 445, 50, 51, 52] or (1100 < scenario < 1119) or (2100 < scenario < 2119) \
            or 2999 < scenario < 3003 or 3099 < scenario < 3103 or 3199 < scenario < 3203:
        M.const.add(expr=M.pyrolysis_in[l, t, 'feedstock'] == A.FEEDSTOCK_SUPPLY[l])
        M.const.add(expr=M.ad_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.htl_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.htc_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.feedstock_to_storage[l, t] == 0)
    elif scenario in [6, 426, 436, 446] or (1200 < scenario < 1219) or (2200 < scenario < 2219):
        M.const.add(expr=M.htl_in[l, t, 'feedstock'] == A.FEEDSTOCK_SUPPLY[l])
        M.const.add(expr=M.pyrolysis_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.htc_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.ad_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.feedstock_to_storage[l, t] == 0)
    elif scenario in [7, 427, 437, 447] or (1300 < scenario < 1319) or (2300 < scenario < 2319):
        M.const.add(expr=M.htc_in[l, t, 'feedstock'] == A.FEEDSTOCK_SUPPLY[l])
        M.const.add(expr=M.pyrolysis_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.htl_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.ad_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.feedstock_to_storage[l, t] == 0)
    elif scenario in [8, 428, 438, 448] or (1400 < scenario < 1419) or (2400 < scenario < 2419):
        M.const.add(expr=M.htc_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.pyrolysis_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.htl_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.ad_in[l, t, 'feedstock'] == A.FEEDSTOCK_SUPPLY[l])
        M.const.add(expr=M.feedstock_to_storage[l, t] == 0)
    elif scenario in [9, 429, 439, 449]:
        M.const.add(expr=M.htc_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.pyrolysis_in[l, t, 'feedstock'] == A.FEEDSTOCK_SUPPLY[l] / 2)
        M.const.add(expr=M.htl_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.ad_in[l, t, 'feedstock'] == A.FEEDSTOCK_SUPPLY[l] / 2)
        M.const.add(expr=M.feedstock_to_storage[l, t] == 0)
    elif scenario in [10103]:
        if A.FEEDSTOCK_SUPPLY[l] > 300: #this accounts for 92.4% of all manure in the state
            M.const.add(expr=M.pyrolysis_in[l, t, 'feedstock'] == A.FEEDSTOCK_SUPPLY[l])
            M.const.add(expr=M.feedstock_to_storage[l, t] == 0)
            M.const.add(expr=sum(M.biochar_from_pyrolysis[l, t, feed, temp, 'CHP'] for feed in M.PyrolysisFeedstocks for temp in M.PyrolysisTemperatures) == 0) #all biochar applied to land
            M.const.add(expr=M.process_capacity[l, 'CHP'] == 1.08731587 * A.FEEDSTOCK_SUPPLY[l])  # to ensure that capex decisions aren't influenced by changes in energy supply of the feedstock
        else:
            M.const.add(expr=M.feedstock_to_storage[l, t] == A.FEEDSTOCK_SUPPLY[l])
            M.const.add(expr=M.pyrolysis_in[l, t, 'feedstock'] == 0)

        # ensure that no other feedstocks are used
        M.const.add(expr=M.ad_in[l, t, 'COD'] == 0)
        M.const.add(expr=M.ad_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.htl_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.htc_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.decision_pyrolysis_temperature[l, t, 'feedstock', 500] == 1)
    elif scenario in [10203]:
        if A.FEEDSTOCK_SUPPLY[l] > 3.2: #this accounts for 99.96% of all manure in the state
            M.const.add(expr=M.pyrolysis_in[l, t, 'feedstock'] == A.FEEDSTOCK_SUPPLY[l])
            M.const.add(expr=M.feedstock_to_storage[l, t] == 0)
            M.const.add(expr=sum(
                M.biochar_from_pyrolysis[l, t, feed, temp, 'CHP'] for feed in M.PyrolysisFeedstocks for temp in
                M.PyrolysisTemperatures) == 0)
            M.const.add(expr=M.process_capacity[l, 'CHP'] == 1.08731587*A.FEEDSTOCK_SUPPLY[l]) #to ensure that capex decisions aren't influenced by changes in energy supply of the feedstock
        else:
            M.const.add(expr=M.feedstock_to_storage[l, t] == A.FEEDSTOCK_SUPPLY[l])
            M.const.add(expr=M.pyrolysis_in[l, t, 'feedstock'] == 0)


        # ensure that no other feedstocks are used
        M.const.add(expr=M.ad_in[l, t, 'COD'] == 0)
        M.const.add(expr=M.ad_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.htl_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.htc_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.decision_pyrolysis_temperature[l, t, 'feedstock', 500] == 1)
    elif scenario in [10]:
        M.const.add(expr=M.htc_in[l, t, 'feedstock'] == A.FEEDSTOCK_SUPPLY[l])
        M.const.add(expr=M.pyrolysis_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.htl_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.ad_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.feedstock_to_storage[l, t] == 0)
        M.const.add(expr=sum(M.hydrochar_from_htc[l, t, feed, temp, 'CHP'] for feed in M.HTCFeedstocks for temp in M.HTCTemperatures) == 0)
        M.const.add(expr=sum(M.hydrochar_from_htc[l, t, feed, temp, 'land'] for feed in M.HTCFeedstocks for temp in M.HTCTemperatures) == 0)
    elif scenario in [11]:
        M.const.add(expr=M.htc_in[l, t, 'feedstock'] == A.FEEDSTOCK_SUPPLY[l])
        M.const.add(expr=M.pyrolysis_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.htl_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.ad_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.feedstock_to_storage[l, t] == 0)
        M.const.add(expr=sum(M.hydrochar_from_htc[l, t, feed, temp, 'CHP'] for feed in M.HTCFeedstocks for temp in M.HTCTemperatures) == 0)
        M.const.add(expr=sum(M.hydrochar_from_htc[l, t, feed, temp, 'market'] for feed in M.HTCFeedstocks for temp in M.HTCTemperatures) == 0)
        for feed in M.HTCFeedstocks:
            M.const.add(expr=M.decision_htc_temperature[l, t, feed, 200] == 1)
    elif scenario in [12]:
        M.const.add(expr=M.htl_in[l, t, 'feedstock'] == A.FEEDSTOCK_SUPPLY[l])
        M.const.add(expr=M.pyrolysis_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.htc_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.ad_in[l, t, 'feedstock'] == 0)
        M.const.add(expr=M.feedstock_to_storage[l, t] == 0)
        M.const.add(expr=sum(M.hydrochar_from_htl[l, t, feed, temp, 'land'] for feed in M.HTLFeedstocks for temp in M.HTLTemperatures) == 0)
        M.const.add(expr=sum(M.hydrochar_from_htl[l, t, feed, temp, 'market'] for feed in M.HTLFeedstocks for temp in M.HTLTemperatures) == 0)
    else:
        M.const.add(expr=M.htl_in[l, t, 'feedstock'] + M.htc_in[l, t, 'feedstock'] + M.ad_in[l, t, 'feedstock'] +
                         M.feedstock_to_storage[l, t] + M.pyrolysis_in[l, t, 'feedstock'] == A.FEEDSTOCK_SUPPLY[l])

    # DIRECT LAND APPLICATION
    M.const.add(expr=M.feedstock_to_storage[l, t] <= M.process_capacity[l, 'Feedstock'])
    M.const.add(expr=M.feedstock_from_storage[l, t] <= M.feedstock_storage[l, t])
    # direct land cost storage size updates
    if t > 0:
        M.const.add(expr=M.feedstock_storage[l, t] == M.feedstock_storage[
            l, t - 1] + M.feedstock_to_storage[l, t - 1] - M.feedstock_from_storage[l, t - 1])
    else:
        M.const.add(expr=M.feedstock_storage[l, t] == 0)
        M.const.add(expr=M.feedstock_from_storage[l, t] == 0)
    if t == A.TIME_PERIODS - 1:
        # raw feedstock has no other disposal pathway other than direct land application, so it must be allowed to
        # store a certain amount for next year
        M.const.add(
            expr=M.feedstock_storage[l, t] <= (A.TIME_PERIODS / A.YEARS - A.LAND_APPLICATION_MONTH - 1) *
                 A.FEEDSTOCK_SUPPLY[l] + 1)  # constrain likely to cause infeasibility issues
    # ensuring that the feedstock storage does not overflow
    M.const.add(expr=M.feedstock_storage[l, t] <= M.feedstock_storage_capacity[l])

    # PYROLYSIS
    M.const.add(expr=sum(M.pyrolysis_in[l, t, feed] for feed in M.PyrolysisFeedstocks) <=
                     M.process_capacity[l, 'Pyrolysis'])
    ##pyrolysis yield conversion
    # only one temperature can be selected for each location at each time period
    for feed in M.PyrolysisFeedstocks:
        M.const.add(expr=sum(M.decision_pyrolysis_temperature[l, t, feed, temp] for temp in
                             M.PyrolysisTemperatures) == 1)
        # getting pyrolysis yields from inputs
        for temp in M.PyrolysisTemperatures:
            for pyro_prod in M.PyrolysisProducts:
                # the amount of yield coming out of the reaction is dependent on yield conversion factors and the amount of material entering
                M.const.add(
                    expr=M.pyrolysis_out[l, t, feed, pyro_prod, temp] / A.PYRO_YIELD[feed, pyro_prod, temp] <=
                         M.pyrolysis_in[l, t, feed])
                M.const.add(
                    expr=M.pyrolysis_out[l, t, feed, pyro_prod, temp] / A.PYRO_YIELD[feed, pyro_prod, temp] <=
                         A.FEEDSTOCK_SUPPLY[l] * M.decision_pyrolysis_temperature[l, t, feed, temp])
                M.const.add(expr=M.pyrolysis_out[l, t, feed, pyro_prod, temp] / A.PYRO_YIELD[
                    feed, pyro_prod, temp] >= 0)
                M.const.add(expr=M.pyrolysis_in[l, t, feed] - A.FEEDSTOCK_SUPPLY[l] * (
                        1 - M.decision_pyrolysis_temperature[l, t, feed, temp]) <= M.pyrolysis_out[
                                     l, t, feed, pyro_prod, temp] / A.PYRO_YIELD[feed, pyro_prod, temp])

                M.const.add(
                    expr=M.pyrolysis_from_storage[l, t, feed, pyro_prod, temp] <=
                         M.pyrolysis_storage[l, t, feed, pyro_prod, temp])
                M.const.add(
                    expr=M.pyrolysis_to_storage[l, t, feed, pyro_prod, temp] <=
                         M.pyrolysis_storage_capacity[l, pyro_prod])
                # storage is empty in the first time period, and is the boundary conditions
                if t > 0:
                    M.const.add(expr=M.pyrolysis_storage[l, t, feed, pyro_prod, temp] ==
                                     M.pyrolysis_storage[l, t - 1, feed, pyro_prod, temp] +
                                     M.pyrolysis_to_storage[l, t - 1, feed, pyro_prod, temp] -
                                     M.pyrolysis_from_storage[l, t - 1, feed, pyro_prod, temp])
                else:
                    M.const.add(expr=M.pyrolysis_storage[l, t, feed, pyro_prod, temp] == 0)
                    M.const.add(expr=M.pyrolysis_from_storage[l, t, feed, pyro_prod, temp] == 0)

                if t == A.TIME_PERIODS - 1:
                    M.const.add(expr=M.pyrolysis_storage[l, t, feed, pyro_prod, temp] == 0)

            # moving products to stage 3 of the model
            M.const.add(
                expr=M.pyrolysis_out[l, t, feed, 'Biochar', temp] + M.pyrolysis_from_storage[
                    l, t, feed, 'Biochar', temp] == sum(
                    M.biochar_from_pyrolysis[l, t, feed, temp, loc] for loc in
                    M.PyroBiocharLocations))
            M.const.add(
                expr=M.pyrolysis_out[l, t, feed, 'Biooil', temp] + M.pyrolysis_from_storage[
                    l, t, feed, 'Biooil', temp] == sum(
                    M.biooil_from_pyrolysis[l, t, feed, temp, loc] for loc in
                    M.PyroBiooilLocations))
            M.const.add(expr=M.pyrolysis_out[l, t, feed, 'Syngas', temp] + M.pyrolysis_from_storage[
                l, t, feed, 'Syngas', temp] == sum(
                M.syngas_from_pyrolysis[l, t, feed, temp, loc] for loc in
                M.PyroSyngasLocations))
            M.const.add(
                expr=M.pyrolysis_out[l, t, feed, 'AP', temp] + M.pyrolysis_from_storage[
                    l, t, feed, 'AP', temp] == sum(
                    M.ap_from_pyrolysis[l, t, feed, temp, loc] for loc in
                    M.PyroAPLocations))

            # linking the edges of the graph for biochar
            M.const.add(expr=M.pyrolysis_to_storage[l, t, feed, 'Biochar', temp] ==
                             M.biochar_from_pyrolysis[l, t, feed, temp, 'storage'])
            M.const.add(expr=M.pyrolysis_to_chp[l, t, feed, 'Biochar', temp] ==
                             A.HHV['Pyrolysis', feed, 'Biochar', temp] * M.biochar_from_pyrolysis[
                                 l, t, feed, temp, 'CHP'])

            # linking the edges of the graph for biooil
            M.const.add(expr=M.pyrolysis_to_storage[
                                 l, t, feed, 'Biooil', temp] == 0)  # biooil cannot be stored
            M.const.add(expr=M.pyrolysis_to_chp[
                                 l, t, feed, 'Biooil', temp] == (
                                     A.HHV['Pyrolysis', feed, 'Biooil', temp] + 3 * 42800) *
                             M.biooil_from_pyrolysis[l, t, feed, temp, 'CHP'])

            # linking the edges of the graph for syngas
            M.const.add(expr=M.pyrolysis_to_storage[l, t, feed, 'Syngas', temp] ==
                             M.syngas_from_pyrolysis[
                                 l, t, feed, temp, 'storage'])
            M.const.add(expr=M.pyrolysis_to_chp[l, t, feed, 'Syngas', temp]
                             == A.HHV['Pyrolysis', feed, 'Syngas', temp] * M.syngas_from_pyrolysis[
                                 l, t, feed, temp, 'CHP'])

            # linking the edges of the graph for aqueous phase
            M.const.add(expr=M.pyrolysis_to_storage[l, t, feed, 'AP', temp] == 0)
            M.const.add(expr=M.pyrolysis_to_chp[l, t, feed, 'AP', temp] == 0)
    # ensuring that the pyrolysis storage does not overflow
    for prod in M.PyrolysisProducts:
        M.const.add(expr=sum(
            M.pyrolysis_storage[l, t, feed, prod, temp] for feed in M.PyrolysisFeedstocks for temp in
            M.PyrolysisTemperatures) <= M.pyrolysis_storage_capacity[l, prod])
    ## HTL
    M.const.add(expr=sum(M.htl_in[l, t, feed] for feed in M.HTLFeedstocks) <= M.process_capacity[l, 'HTL'])
    ##htl yield conversion
    # only one temperature can be selected for each location at each time period
    for feed in M.HTLFeedstocks:
        M.const.add(
            expr=sum(M.decision_htl_temperature[l, t, feed, temp] for temp in M.HTLTemperatures) == 1)
        # getting pyrolysis yields from inputs
        for prod in M.HTLProducts:
            for temp in M.HTLTemperatures:
                # the amount of yield coming out of the reaction is dependent on yield conversion factors and the amount of material entering
                M.const.add(
                    expr=M.htl_out[l, t, feed, prod, temp] / A.HTL_YIELD[feed, prod, temp] <=
                         M.htl_in[l, t, feed])
                M.const.add(
                    expr=M.htl_out[l, t, feed, prod, temp] / A.HTL_YIELD[feed, prod, temp] <=
                         A.FEEDSTOCK_SUPPLY[l] * M.decision_htl_temperature[l, t, feed, temp])
                M.const.add(expr=M.htl_out[l, t, feed, prod, temp] / A.HTL_YIELD[
                    feed, prod, temp] >= 0)
                M.const.add(expr=M.htl_in[l, t, feed] - A.FEEDSTOCK_SUPPLY[l] * (
                        1 - M.decision_htl_temperature[l, t, feed, temp]) <= M.htl_out[
                                     l, t, feed, prod, temp] / A.HTL_YIELD[feed, prod, temp])
                M.const.add(
                    expr=M.htl_from_storage[l, t, feed, prod, temp] <= M.htl_storage[
                        l, t, feed, prod, temp])
                M.const.add(
                    expr=M.htl_to_storage[l, t, feed, prod, temp] <= M.htl_storage_capacity[
                        l, prod])

                # storage is empty in the first time period, and is the boundary conditions
                if (t > 0):
                    M.const.add(expr=M.htl_storage[l, t, feed, prod, temp] == M.htl_storage[
                        l, t - 1, feed, prod, temp] + M.htl_to_storage[l, t - 1, feed, prod, temp] -
                                     M.htl_from_storage[l, t - 1, feed, prod, temp])
                else:
                    M.const.add(expr=M.htl_storage[l, t, feed, prod, temp] == 0)
                    M.const.add(expr=M.htl_from_storage[l, t, feed, prod, temp] == 0)

                if t == A.TIME_PERIODS - 1:
                    M.const.add(expr=M.htl_storage[l, t, feed, prod, temp] == 0)
    for feed in M.HTLFeedstocks:
        for temp in M.HTLTemperatures:
            # moving products to stage 3 of the model
            M.const.add(expr=M.htl_out[l, t, feed, 'Hydrochar', temp] + M.htl_from_storage[
                l, t, feed, 'Hydrochar', temp] == sum(
                M.hydrochar_from_htl[l, t, feed, temp, loc] for loc in M.HTLHydrocharLocations))
            M.const.add(expr=M.htl_out[l, t, feed, 'Biooil', temp] == sum(
                M.biooil_from_htl[l, t, feed, temp, loc] for loc in M.HTLBiooilLocations))
            M.const.add(
                expr=M.htl_out[l, t, feed, 'GP', temp] + M.htl_from_storage[l, t, feed, 'GP', temp] == sum(
                    M.gp_from_htl[l, t, feed, temp, loc] for loc in M.HTLGPLocations))
            M.const.add(
                expr=M.htl_out[l, t, feed, 'AP', temp] + M.htl_from_storage[l, t, feed, 'AP', temp] == sum(
                    M.ap_from_htl[l, t, feed, temp, loc] for loc in M.HTLAPLocations))

            # linking the edges of the graph for hydrochar
            M.const.add(expr=M.htl_to_storage[l, t, feed, 'Hydrochar', temp] == M.hydrochar_from_htl[
                l, t, feed, temp, 'storage'])
            M.const.add(
                expr=M.htl_to_chp[l, t, feed, 'Hydrochar', temp] == A.HHV['HTL', feed, 'Hydrochar', temp] *
                     M.hydrochar_from_htl[l, t, feed, temp, 'CHP'])

            # linking the edges of the graph for biooil
            M.const.add(expr=M.htl_to_storage[l, t, feed, 'Biooil', temp] == 0)  # biooil cannot be stored
            M.const.add(
                expr=M.htl_to_chp[l, t, feed, 'Biooil', temp] == (
                        A.HHV['HTL', feed, 'Biooil', temp] + 3 * 42800) *
                     M.biooil_from_htl[l, t, feed, temp, 'CHP'])

            # linking the edges of the graph for GP
            M.const.add(expr=M.htl_to_storage[l, t, feed, 'GP', temp] == 0)
            M.const.add(expr=M.htl_to_chp[l, t, feed, 'GP', temp] == 0)

            # linking the edges of the graph for AP
            M.const.add(expr=M.htl_to_storage[l, t, feed, 'AP', temp] == M.ap_from_htl[
                l, t, feed, temp, 'storage'])
            M.const.add(expr=M.htl_to_chp[l, t, feed, 'AP', temp] == 0)
    # ensuring that the HTL storage does not overflow
    for prod in M.HTLProducts:
        M.const.add(expr=sum(M.htl_storage[l, t, feed, prod, temp] for feed in M.HTLFeedstocks
                             for temp in M.HTLTemperatures) <= M.htl_storage_capacity[l, prod])
    ## HTC
    # products entering htc must be under capacity
    M.const.add(expr=sum(M.htc_in[l, t, feed] for feed in M.HTCFeedstocks) <= M.process_capacity[l, 'HTC'])
    ##htc yield conversion
    # only one temperature can be selected for each location at each time period
    for feed in M.HTCFeedstocks:
        M.const.add(
            expr=sum(M.decision_htc_temperature[l, t, feed, temp] for temp in M.HTCTemperatures) == 1)
        # getting HTC yields from inputs
        for temp in M.HTCTemperatures:
            for prod in M.HTCProducts:
                # the amount of yield coming out of the reaction is dependent on yield conversion factors and the amount of material entering
                M.const.add(
                    expr=M.htc_out[l, t, feed, prod, temp] / A.HTC_YIELD[feed, prod, temp] <=
                         M.htc_in[l, t, feed])
                M.const.add(
                    expr=M.htc_out[l, t, feed, prod, temp] / A.HTC_YIELD[feed, prod, temp] <=
                         A.FEEDSTOCK_SUPPLY[l] * M.decision_htc_temperature[l, t, feed, temp])
                M.const.add(expr=M.htc_out[l, t, feed, prod, temp] / A.HTC_YIELD[
                    feed, prod, temp] >= 0)
                M.const.add(expr=M.htc_in[l, t, feed] - A.FEEDSTOCK_SUPPLY[l] * (
                        1 - M.decision_htc_temperature[l, t, feed, temp]) <= M.htc_out[
                                     l, t, feed, prod, temp] / A.HTC_YIELD[feed, prod, temp])

                M.const.add(
                    expr=M.htc_from_storage[l, t, feed, prod, temp] <= M.htc_storage[
                        l, t, feed, prod, temp])
                M.const.add(
                    expr=M.htc_to_storage[l, t, feed, prod, temp] <= M.htc_storage_capacity[
                        l, prod])

                # storage is empty in the first time period, and is the boundary conditions
                if t > 0:
                    M.const.add(expr=M.htc_storage[l, t, feed, prod, temp] == M.htc_storage[
                        l, t - 1, feed, prod, temp] + M.htc_to_storage[l, t - 1, feed, prod, temp] -
                                     M.htc_from_storage[l, t - 1, feed, prod, temp])
                else:
                    M.const.add(expr=M.htc_storage[l, t, feed, prod, temp] == 0)
                    M.const.add(expr=M.htc_from_storage[l, t, feed, prod, temp] == 0)

                if t == A.TIME_PERIODS - 1:
                    M.const.add(expr=M.htc_storage[l, t, feed, prod, temp] == 0)

            # moving products to stage 3 of the model
            M.const.add(expr=M.htc_out[l, t, feed, 'Hydrochar', temp] + M.htc_from_storage[
                l, t, feed, 'Hydrochar', temp] == sum(
                M.hydrochar_from_htc[l, t, feed, temp, loc] for loc in M.HTCHydrocharLocations))
            M.const.add(
                expr=M.htc_out[l, t, feed, 'GP', temp] + M.htc_from_storage[l, t, feed, 'GP', temp] == sum(
                    M.gp_from_htc[l, t, feed, temp, loc] for loc in M.HTCGPLocations))
            M.const.add(
                expr=M.htc_out[l, t, feed, 'AP', temp] + M.htc_from_storage[l, t, feed, 'AP', temp] == sum(
                    M.ap_from_htc[l, t, feed, temp, loc] for loc in M.HTCAPLocations))

            # linking the edges of the graph for hydrochar
            M.const.add(expr=M.htc_to_storage[l, t, feed, 'Hydrochar', temp] == M.hydrochar_from_htc[
                l, t, feed, temp, 'storage'])
            M.const.add(
                expr=M.htc_to_chp[l, t, feed, 'Hydrochar', temp] == A.HHV['HTC', feed, 'Hydrochar', temp] *
                     M.hydrochar_from_htc[l, t, feed, temp, 'CHP'])

            # linking the edges of the graph for gaseous phase
            M.const.add(expr=M.htc_to_storage[l, t, feed, 'GP', temp] == 0)
            M.const.add(expr=M.htc_to_chp[l, t, feed, 'GP', temp] == 0)

            # linking the edges of the graph for AP
            M.const.add(expr=M.htc_to_storage[l, t, feed, 'AP', temp] == M.ap_from_htc[
                l, t, feed, temp, 'storage'])
            M.const.add(expr=M.htc_to_chp[l, t, feed, 'AP', temp] == 0)

            # ensuring that the HTC storage does not overflow
        for prod in M.HTCProducts:
            M.const.add(
                expr=sum(M.htc_storage[l, t, feed, prod, temp] for feed in M.HTCFeedstocks for temp in
                         M.HTCTemperatures) <= M.htc_storage_capacity[l, prod])
    ## AD
    # products entering ad must be under capacity
    for stage in M.ADStages:
        for feed in M.ADFeedstocks:
            # glovers linearization for ad_in_glovers = ad_in*decision_ad_stage
            M.const.add(
                expr=M.ad_in[l, t, feed] - A.FEEDSTOCK_SUPPLY[l] * (1 - M.decision_ad_stage[l, stage]) <=
                     M.ad_in_glovers[l, t, feed, stage])
            M.const.add(expr=M.ad_in_glovers[l, t, feed, stage] <= A.FEEDSTOCK_SUPPLY[l] * M.decision_ad_stage[
                l, stage])
            M.const.add(expr=M.ad_in_glovers[l, t, feed, stage] <= M.ad_in[l, t, feed])
            M.const.add(expr=M.ad_in_glovers[l, t, feed, stage] >= 0)

        M.const.add(
            expr=sum(M.ad_in_glovers[l, t, feed, stage] * A.LOADING[feed, stage] for feed in M.ADFeedstocks)
                 == M.ad_capacity[l, stage])
    M.const.add(expr=sum(M.ad_capacity[l, stage] for stage in M.ADStages) == M.process_capacity[l, 'AD'])

    if scenario in [5, 425, 435, 445, 50, 51, 52] or (1100 < scenario < 1119) or (2100 < scenario < 2119) \
            or 2999 < scenario < 3003 or 3099 < scenario < 3103 or 3199 < scenario < 3203:
        # scenarios in which we want pyrolysis only, no AD
        M.const.add(expr=M.ad_in[l, t, 'COD'] == 0)
    else:
        M.const.add(expr=M.ad_in[l, t, 'COD'] == sum(A.COD['Pyrolysis', feed, 'AP', temp] * M.ap_from_pyrolysis[
            l, t, feed, temp, 'AD'] for feed in M.PyrolysisFeedstocks
                                                     for temp in M.PyrolysisTemperatures))  # units: tons COD
    ##ad yield conversion
    # getting AD yields from inputs
    for prod in M.ADProducts:
        for stage in M.ADStages:
            # the amount of yield coming out of the reaction is dependent on yield conversion factors and the amount of material entering
            if t == 0:
                M.const.add(expr=M.ad_out[l, t, prod, stage] == 0)  # no yield from first month of ad
            else:
                M.const.add(expr=M.ad_out[l, t, prod, stage] == sum(M.ad_in_glovers[l, t - 1, feed, stage] *
                                                                    A.AD_YIELD[feed, prod, stage] for feed
                                                                    in M.ADFeedstocks))
                # yields from ad take time to materialize, biogas yield is in Nm^3, digestate is in kg

            M.const.add(expr=M.ad_from_storage[l, t, prod, stage] <= M.ad_storage[l, t, prod, stage])
            M.const.add(expr=M.ad_to_storage[l, t, prod, stage] <= M.ad_storage_capacity[l, prod])

            # storage is empty in the first time period, and is the boundary conditions
            if t > 0:
                M.const.add(expr=M.ad_storage[l, t, prod, stage] == M.ad_storage[
                    l, t - 1, prod, stage] + M.ad_to_storage[l, t - 1, prod, stage] -
                                 M.ad_from_storage[l, t - 1, prod, stage])
            else:
                M.const.add(expr=M.ad_storage[l, t, prod, stage] == 0)
                M.const.add(expr=M.ad_from_storage[l, t, prod, stage] == 0)

            if t == A.TIME_PERIODS - 1:
                M.const.add(expr=M.ad_storage[l, t, prod, stage] == 0)
    for stage in M.ADStages:
        # moving products to stage 3 of the model
        M.const.add(expr=M.ad_out[l, t, 'digestate', stage] + M.ad_from_storage[
            l, t, 'digestate', stage] == sum(
            M.digestate_from_ad[l, t, stage, loc] for loc in M.ADDigestateLocations))
        M.const.add(expr=M.ad_out[l, t, 'biogas', stage] + M.ad_from_storage[
            l, t, 'biogas', stage] == sum(
            M.biogas_from_ad[l, t, stage, loc] for loc in M.ADBiogasLocations))

        # linking the edges of the graph for digestate
        # if scenario is 1, then digestate can only be disposed of
        if scenario in [1, 421, 431, 441]:
            M.const.add(expr=M.digestate_from_ad[l, t, stage, 'storage'] == 0)
            M.const.add(expr=M.digestate_from_ad[l, t, stage, 'land'] == 0)

        M.const.add(expr=M.ad_to_storage[l, t, 'digestate', stage] == M.digestate_from_ad[
            l, t, stage, 'storage'])
        M.const.add(expr=M.ad_to_chp[l, t, 'digestate', stage] == 0)

        M.const.add(expr=M.ad_to_storage[l, t, 'biogas', stage] == M.biogas_from_ad[
            l, t, stage, 'storage'])
        M.const.add(
            expr=M.ad_to_chp[l, t, 'biogas', stage] == A.HHV['methane'] * M.biogas_from_ad[
                l, t, stage, 'CHP'])  # units MJ

        # ensuring the AD storage does not overflow
        for prod in M.ADProducts:
            M.const.add(expr=sum(M.ad_storage[l, t, prod, stage] for stage in
                                 M.ADStages) <= M.ad_storage_capacity[l, prod])
    # CHP CONSTRAINTS
    M.const.add(expr=M.chp_in[l, t] * A.MJ_TO_KW / A.HOURS_PER_PERIOD <= M.process_capacity[l, 'CHP'])
    # different available products in different scenarios
    M.const.add(expr=M.chp_in[l, t] == sum(
        M.pyrolysis_to_chp[l, t, feedstock, pyro_prod, temp]
        for temp in M.PyrolysisTemperatures for pyro_prod in M.PyrolysisProducts for feedstock in
        M.PyrolysisFeedstocks) +
                     sum(M.htl_to_chp[l, t, feedstock, htl_prod, temp]
                         for feedstock in M.HTLFeedstocks for htl_prod
                         in M.HTLProducts for temp in M.HTLTemperatures) +
                     sum(M.htc_to_chp[l, t, feedstock, htc_prod, temp]
                         for feedstock in M.HTCFeedstocks for htc_prod
                         in M.HTCProducts for temp in M.HTCTemperatures) +
                     sum(M.ad_to_chp[l, t, ad_prod, temp]
                         for ad_prod in M.ADProducts for temp in
                         M.ADStages))  # units in MJ
    # CHP yields heat and power depends on the types of feedstocks used
    M.const.add(
        expr=sum(M.chp_out[l, t, tech, 'heat'] for tech in M.Technology) + M.chp_market[l, t, 'heat'] ==
             A.CHP_HEAT_EFFICIENCY * M.chp_in[l, t])  # units of J
    M.const.add(expr=sum(M.chp_out[l, t, tech, 'electricity'] for tech in M.Technology) +
                     M.chp_market[l, t, 'electricity'] == A.CHP_ELECTRICITY_EFFICIENCY * M.chp_in[
                         l, t])  # units of kWh out, MJ in


def add_variables(M):
    """
    Adds variables to the model
    :param M: the model
    :return:
    """
    ### VARIABLES
    ## LAND APPLICATION VARIABLES
    M.feedstock_storage_capacity = pyo.Var(M.Location, initialize=0, within=pyo.NonNegativeReals)
    M.feedstock_to_storage = pyo.Var(M.Location, M.Time, initialize=0, within=pyo.NonNegativeReals)
    M.feedstock_storage = pyo.Var(M.Location, M.Time, initialize=0, within=pyo.NonNegativeReals)
    M.feedstock_from_storage = pyo.Var(M.Location, M.Time, initialize=0, within=pyo.NonNegativeReals)

    ## PYROLYSIS VARIABLES
    M.decision_pyrolysis_temperature = pyo.Var(M.Location, M.Time, M.PyrolysisFeedstocks, M.PyrolysisTemperatures,
                                               initialize=0, within=pyo.Binary)
    M.pyrolysis_out = pyo.Var(M.Location, M.Time, M.PyrolysisFeedstocks, M.PyrolysisProducts,
                              M.PyrolysisTemperatures, initialize=0, within=pyo.NonNegativeReals)
    M.pyrolysis_in = pyo.Var(M.Location, M.Time, M.PyrolysisFeedstocks, initialize=0,
                             within=pyo.NonNegativeReals)
    M.pyrolysis_to_storage = pyo.Var(M.Location, M.Time, M.PyrolysisFeedstocks, M.PyrolysisProducts,
                                     M.PyrolysisTemperatures, initialize=0, within=pyo.NonNegativeReals)
    M.pyrolysis_to_chp = pyo.Var(M.Location, M.Time, M.PyrolysisFeedstocks, M.PyrolysisProducts,
                                 M.PyrolysisTemperatures, initialize=0, within=pyo.NonNegativeReals)
    M.pyrolysis_from_storage = pyo.Var(M.Location, M.Time, M.PyrolysisFeedstocks, M.PyrolysisProducts,
                                       M.PyrolysisTemperatures, initialize=0, within=pyo.NonNegativeReals)
    M.pyrolysis_storage = pyo.Var(M.Location, M.Time, M.PyrolysisFeedstocks, M.PyrolysisProducts,
                                  M.PyrolysisTemperatures, initialize=0, within=pyo.NonNegativeReals)
    M.pyrolysis_storage_capacity = pyo.Var(M.Location, M.PyrolysisProducts, initialize=0,
                                           within=pyo.NonNegativeReals, bounds=(0, 200000))
    M.biochar_from_pyrolysis = pyo.Var(M.Location, M.Time, M.PyrolysisFeedstocks, M.PyrolysisTemperatures,
                                       M.PyroBiocharLocations, initialize=0, within=pyo.NonNegativeReals)
    M.biooil_from_pyrolysis = pyo.Var(M.Location, M.Time, M.PyrolysisFeedstocks, M.PyrolysisTemperatures,
                                      M.PyroBiooilLocations, initialize=0, within=pyo.NonNegativeReals)
    M.syngas_from_pyrolysis = pyo.Var(M.Location, M.Time, M.PyrolysisFeedstocks, M.PyrolysisTemperatures,
                                      M.PyroSyngasLocations, initialize=0, within=pyo.NonNegativeReals)
    M.ap_from_pyrolysis = pyo.Var(M.Location, M.Time, M.PyrolysisFeedstocks, M.PyrolysisTemperatures,
                                  M.PyroAPLocations, initialize=0, within=pyo.NonNegativeReals)
    M.pyro_storage_cost = pyo.Var(M.Location, M.PyrolysisProducts, initialize=0, within=pyo.NonNegativeReals)

    ## HTL VARIABLES
    M.decision_htl_temperature = pyo.Var(M.Location, M.Time, M.HTLFeedstocks, M.HTLTemperatures,
                                         initialize=0, within=pyo.Binary)
    M.htl_out = pyo.Var(M.Location, M.Time, M.HTLFeedstocks, M.HTLProducts, M.HTLTemperatures,
                        initialize=0, within=pyo.NonNegativeReals)
    M.htl_in = pyo.Var(M.Location, M.Time, M.HTLFeedstocks, initialize=0, within=pyo.NonNegativeReals)
    M.htl_to_storage = pyo.Var(M.Location, M.Time, M.HTLFeedstocks, M.HTLProducts,
                               M.HTLTemperatures, initialize=0, within=pyo.NonNegativeReals)
    M.htl_to_chp = pyo.Var(M.Location, M.Time, M.HTLFeedstocks, M.HTLProducts,
                           M.HTLTemperatures, initialize=0, within=pyo.NonNegativeReals)
    M.htl_from_storage = pyo.Var(M.Location, M.Time, M.HTLFeedstocks, M.HTLProducts,
                                 M.HTLTemperatures, initialize=0, within=pyo.NonNegativeReals)
    M.htl_storage = pyo.Var(M.Location, M.Time, M.HTLFeedstocks, M.HTLProducts,
                            M.HTLTemperatures, initialize=0, within=pyo.NonNegativeReals)
    M.htl_storage_capacity = pyo.Var(M.Location, M.HTLProducts, initialize=0, within=pyo.NonNegativeReals,
                                     bounds=(0, 200000))
    M.hydrochar_from_htl = pyo.Var(M.Location, M.Time, M.HTLFeedstocks, M.HTLTemperatures,
                                   M.HTLHydrocharLocations, initialize=0, within=pyo.NonNegativeReals)
    M.biooil_from_htl = pyo.Var(M.Location, M.Time, M.HTLFeedstocks, M.HTLTemperatures,
                                M.HTLBiooilLocations, initialize=0, within=pyo.NonNegativeReals)
    M.gp_from_htl = pyo.Var(M.Location, M.Time, M.HTLFeedstocks, M.HTLTemperatures,
                            M.HTLGPLocations, initialize=0, within=pyo.NonNegativeReals)
    M.ap_from_htl = pyo.Var(M.Location, M.Time, M.HTLFeedstocks, M.HTLTemperatures,
                            M.HTLAPLocations, initialize=0, within=pyo.NonNegativeReals)
    M.htl_storage_cost = pyo.Var(M.Location, M.HTLProducts, initialize=0, within=pyo.NonNegativeReals)
    ## HTC VARIABLES
    M.decision_htc_temperature = pyo.Var(M.Location, M.Time, M.HTCFeedstocks, M.HTCTemperatures,
                                         initialize=0, within=pyo.Binary)
    M.htc_out = pyo.Var(M.Location, M.Time, M.HTCFeedstocks, M.HTCProducts, M.HTCTemperatures,
                        initialize=0, within=pyo.NonNegativeReals)
    M.htc_in = pyo.Var(M.Location, M.Time, M.HTCFeedstocks, initialize=0, within=pyo.NonNegativeReals)
    M.htc_to_storage = pyo.Var(M.Location, M.Time, M.HTCFeedstocks, M.HTCProducts,
                               M.HTCTemperatures, initialize=0, within=pyo.NonNegativeReals)
    M.htc_to_chp = pyo.Var(M.Location, M.Time, M.HTCFeedstocks, M.HTCProducts,
                           M.HTCTemperatures, initialize=0, within=pyo.NonNegativeReals)
    M.htc_from_storage = pyo.Var(M.Location, M.Time, M.HTCFeedstocks, M.HTCProducts,
                                 M.HTCTemperatures, initialize=0, within=pyo.NonNegativeReals)
    M.htc_storage = pyo.Var(M.Location, M.Time, M.HTCFeedstocks, M.HTCProducts,
                            M.HTCTemperatures, initialize=0, within=pyo.NonNegativeReals)
    M.htc_storage_capacity = pyo.Var(M.Location, M.HTCProducts, initialize=0, within=pyo.NonNegativeReals,
                                     bounds=(0, 200000))
    M.hydrochar_from_htc = pyo.Var(M.Location, M.Time, M.HTCFeedstocks, M.HTCTemperatures,
                                   M.HTCHydrocharLocations, initialize=0, within=pyo.NonNegativeReals)
    M.gp_from_htc = pyo.Var(M.Location, M.Time, M.HTCFeedstocks, M.HTCTemperatures,
                            M.HTCGPLocations, initialize=0, within=pyo.NonNegativeReals)
    M.ap_from_htc = pyo.Var(M.Location, M.Time, M.HTCFeedstocks, M.HTCTemperatures,
                            M.HTCAPLocations, initialize=0, within=pyo.NonNegativeReals)
    M.htc_storage_cost = pyo.Var(M.Location, M.HTCProducts, initialize=0, within=pyo.NonNegativeReals)

    ## AD VARIABLES
    M.ad_in_glovers = pyo.Var(M.Location, M.Time, M.ADFeedstocks, M.ADStages, initialize=0, within=pyo.NonNegativeReals)
    M.decision_ad_stage = pyo.Var(M.Location, M.ADStages, initialize=0,
                                  within=pyo.Binary)  # can only choose one reactor size
    M.ad_capacity = pyo.Var(M.Location, M.ADStages, initialize=0, within=pyo.NonNegativeReals)
    M.ad_out = pyo.Var(M.Location, M.Time, M.ADProducts, M.ADStages,
                       initialize=0, within=pyo.NonNegativeReals)
    M.ad_in = pyo.Var(M.Location, M.Time, M.ADFeedstocks, initialize=0, within=pyo.NonNegativeReals)
    M.ad_to_storage = pyo.Var(M.Location, M.Time, M.ADProducts,
                              M.ADStages, initialize=0, within=pyo.NonNegativeReals)
    M.ad_to_chp = pyo.Var(M.Location, M.Time, M.ADProducts,
                          M.ADStages, initialize=0, within=pyo.NonNegativeReals)
    M.ad_from_storage = pyo.Var(M.Location, M.Time, M.ADProducts,
                                M.ADStages, initialize=0, within=pyo.NonNegativeReals)
    M.ad_storage = pyo.Var(M.Location, M.Time, M.ADProducts,
                           M.ADStages, initialize=0, within=pyo.NonNegativeReals)
    M.ad_storage_capacity = pyo.Var(M.Location, M.ADProducts, initialize=0, within=pyo.NonNegativeReals,
                                    bounds=(0, 200000))
    M.ad_storage_cost = pyo.Var(M.Location, M.ADProducts, initialize=0, within=pyo.NonNegativeReals)
    M.digestate_from_ad = pyo.Var(M.Location, M.Time, M.ADStages,
                                  M.ADDigestateLocations, initialize=0, within=pyo.NonNegativeReals)
    M.biogas_from_ad = pyo.Var(M.Location, M.Time, M.ADStages,
                               M.ADBiogasLocations, initialize=0, within=pyo.NonNegativeReals)
    ##CHP
    M.chp_in = pyo.Var(M.Location, M.Time, initialize=0, within=pyo.NonNegativeReals)
    M.chp_out = pyo.Var(M.Location, M.Time, M.Technology, M.CHPProducts, initialize=0, within=pyo.NonNegativeReals)
    M.chp_market = pyo.Var(M.Location, M.Time, M.CHPProducts, initialize=0, within=pyo.NonNegativeReals)

    ##FACTOR VARIABLES
    M.opex_costs = pyo.Var(M.Location, M.Time, M.Technology, M.OPEXSubCosts, initialize=0, within=pyo.NonNegativeReals)
    M.opex_revenues = pyo.Var(M.Location, M.Time, M.Technology, M.OPEXSubRevenues, initialize=0,
                              within=pyo.NonNegativeReals)
    M.inputs = pyo.Var(M.Location, M.Time, M.Technology, M.InputProducts, initialize=0, within=pyo.NonNegativeReals)
    M.process_capex = pyo.Var(M.Location, M.Technology, initialize=0, within=pyo.NonNegativeReals)
    M.storage_capex = pyo.Var(M.Location, M.Technology, initialize=0, within=pyo.NonNegativeReals)
    M.npv = pyo.Var(M.Location, initialize=0, within=pyo.Reals, bounds=(-1e12, 1e12))

    M.process_capacity = pyo.Var(M.Location, M.Technology, initialize=0,
                                 within=pyo.NonNegativeReals, bounds=(0, 200000))  # process capacity for each process
    M.avoided_fertilizers = pyo.Var(M.Location, M.Time, M.Technology, M.AvoidedFertilizers, initialize=0,
                                    within=pyo.NonNegativeReals)

    ## LCA VARIABLES
    M.LCA_midpoints = pyo.Var(M.Location, M.LCATypes, M.ALCAInputs, M.LCAMidpointCat, within=pyo.Reals, initialize=0)
    M.total_LCA_midpoints = pyo.Var(M.Location, M.LCATypes, M.LCAMidpointCat, within=pyo.Reals, initialize=0)

    ##Purchased goods
    M.purchased_power = pyo.Var(M.Location, M.Time, M.Technology, initialize=0, within=pyo.NonNegativeReals)
    M.purchased_fuel = pyo.Var(M.Location, M.Time, M.Technology, initialize=0, within=pyo.NonNegativeReals)


def add_sets(A, M, j):
    """
    Defines sets for the model
    :param A: CLCO parameters
    :param M: the model itself
    :param j: the location number
    :return:
    """
    ### SETS
    # technology is the set of technologies in stage 2 of the model
    M.Location = pyo.Set(initialize=[j])

    M.Technology = pyo.Set(initialize=['Pyrolysis', 'AD', 'HTL', 'HTC', 'CHP', 'Feedstock'])

    M.Time = pyo.Set(initialize=[j for j in range(A.TIME_PERIODS)])

    M.AvoidedFertilizers = pyo.Set(initialize=['N', 'P', 'K'])

    M.LCATypes = pyo.Set(initialize=['ALCA', 'CLCA'])

    M.PyrolysisProducts = pyo.Set(initialize=['Biochar', 'Syngas', 'Biooil', 'AP'])
    M.PyrolysisTemperatures = pyo.Set(initialize=[400, 450, 500, 550, 600, 700, 800])
    M.PyrolysisFeedstocks = pyo.Set(initialize=['feedstock'])
    M.PyroBiocharLocations = pyo.Set(initialize=['storage', 'CHP', 'land', 'disposal'])
    M.PyroBiooilLocations = pyo.Set(initialize=['CHP', 'market'])
    M.PyroSyngasLocations = pyo.Set(initialize=['storage', 'CHP', 'disposal'])
    M.PyroAPLocations = pyo.Set(initialize=['AD', 'disposal'])

    M.HTLProducts = pyo.Set(initialize=['Hydrochar', 'GP', 'AP', 'Biooil'])
    M.HTLTemperatures = pyo.Set(initialize=[350])
    M.HTLFeedstocks = pyo.Set(initialize=['feedstock'])
    M.HTLHydrocharLocations = pyo.Set(initialize=['storage', 'land', 'CHP', 'market', 'disposal'])
    M.HTLBiooilLocations = pyo.Set(initialize=['market', 'CHP'])
    M.HTLGPLocations = pyo.Set(initialize=['disposal'])
    M.HTLAPLocations = pyo.Set(initialize=['storage', 'disposal'])

    M.HTCProducts = pyo.Set(initialize=['Hydrochar', 'GP', 'AP'])
    M.HTCTemperatures = pyo.Set(initialize=[180, 200, 220, 250])
    M.HTCFeedstocks = pyo.Set(initialize=['feedstock'])
    M.HTCHydrocharLocations = pyo.Set(initialize=['storage', 'CHP', 'land', 'market', 'disposal'])
    M.HTCGPLocations = pyo.Set(initialize=['disposal'])
    M.HTCAPLocations = pyo.Set(initialize=['storage', 'disposal'])

    M.ADProducts = pyo.Set(initialize=['digestate', 'biogas'])
    M.ADStages = pyo.Set(initialize=[1.5, 3, 4.5])
    M.ADFeedstocks = pyo.Set(initialize=['feedstock', 'COD'])
    M.ADBiogasLocations = pyo.Set(initialize=['storage', 'CHP', 'disposal'])
    M.ADDigestateLocations = pyo.Set(initialize=['storage', 'land', 'disposal'])

    M.OPEXSubCosts = pyo.Set(
        initialize=['heat', 'electricity', 'disposal', 'transportation', 'water',
                    'labor', 'diesel', 'TPC'])
    M.OPEXSubRevenues = pyo.Set(
        initialize=['avoided fertilizer', 'bio oil', 'avoided coal', 'electricity', 'potting media', 'incentive 1',
                    'incentive 2'])
    M.InputProducts = pyo.Set(initialize=['heat', 'electricity', 'water', 'bio-oil diesel', 'transportation'])
    M.CHPProducts = pyo.Set(initialize=['heat', 'electricity'])

    # LCA categories
    M.ALCAInputs = pyo.Set(initialize=['natural gas', 'grid electricity', 'diesel', 'water', 'biochar-chp',
                                       'biochar-land', 'biochar-disposal', 'pyro-bio-oil-chp', 'syngas-chp',
                                       'syngas-disposal', 'pyro-ap-disposal', 'htl-hydrochar-land',
                                       'htl-hydrochar-chp', 'htl-hydrochar-disposal', 'htl-bio-oil-chp',
                                       'htl-gp-disposal', 'htl-ap-disposal', 'htc-hydrochar-land', 'htc-hydrochar-chp',
                                       'htc-hydrochar-disposal', 'htc-gp-disposal', 'htc-ap-disposal', 'digestate-land',
                                       'digestate-disposal', 'biogas-disposal', 'biogas-chp', 'manure-land',
                                       'facility construction', 'N fertilizer', 'P fertilizer', 'K fertilizer',
                                       'storage-facility-solids', 'storage-facility-liquids', "biochar market",
                                       "bio-oil market", "hydrochar market", "avoided electricity"])
    M.LCAMidpointCat = pyo.Set(
        initialize=['acidification', 'climate change', 'ecotoxicity: freshwater', 'ecotoxicity: marine',
                    'ecotoxicity: terrestrial',
                    'energy resources', 'eutrophication: freshwater', 'eutrophication: marine',
                    'human toxicity: carcinogenic',
                    'human toxicity: non -carcinogenic', 'ionising radiation', 'land use', 'material resources',
                    'ozone depletion', 'particulate matter formation', 'photochemical oxidant formation: human health',
                    'photochemical oxidant formation: terrestrial ecosystems', 'water use'])

    M.alpha = pyo.Param(default=0, mutable=True)


def save_plot(scenario, FLP=False, midpoint=""):
    """
    Saves the pareto plot for quick investigations
    :param scenario: the scenario number
    :param FLP: whether or not this is being used to solve the facility location problem
    :param midpoint: the lca_midpoint for hte Pareto front
    :return:
    """
    if FLP:
        filename = "data/FLP/data/" + str(scenario) + 'aws'
    else:
        filename = "data/" + str(scenario) + "/" + midpoint.replace(" ", "").replace(":", "-") + '_aws'
        print("filename", filename)
    return filename


def print_model(scenario, model, location_num, model_type="TEA", FLP=False, midpoint=""):
    """
    Prints the relevant data in the model out to a .csv file
    :param scenario: the sceanrio number
    :param model: the model being used
    :param location_num: the location number
    :param model_type: the model type TEA/LCA
    :param FLP: boolean for if facility location problem is being conducted
    :param midpoint: the LCA midpoints being used
    :return:
    """
    # construct the filenames
    if FLP:
        filepath = "data/FLP/" + model_type + "/data/"
    else:
        filepath = "data/" + str(scenario) + "/" + model_type

    # check to see if path exists, and if it doesn't, create it
    if not os.path.exists(filepath):
        os.makedirs(filepath)

    filename = filepath + "/location" + str(location_num) + midpoint.replace(" ", "").replace(":", "-") + '_data.csv'

    with open(filename, 'w', encoding='UTF8', newline='') as f:
        write = csv.writer(f)

        # identify the necessary categories
        labels = []
        if model_type == "TEA":
            labels.append("opex_costs")
            labels.append("opex_revenues")
            labels.append("process_capex")
            labels.append("storage_capex")
        elif model_type == "LCA":
            labels.append("LCA_midpoints")

        # print out all the data
        for v in model.component_objects(pyo.Var, active=True):
            if str(v) in labels:
                row = ["Variable", v]
                write.writerow(row)
                for index in v:
                    row = [v, index, pyo.value(v[index])]
                    write.writerow(row)


if __name__ == '__main__':
    '''
    SCENARIO LIST:
    1: AD + CHP w/ disposal of digestate, NPV max, county level
    2: Direct Land Application, NPV max, county level
    3: Optimal County Level Results, NPV max, county level
    4: FLP, NPV max, county level
    5: Pyrolysis + CHP, NPV max, county level
    6: HTL + CHP, NPV max, county level
    7: HTC + CHP, NPV max, county level
    8: AD + CHP, NPV max, county level
    9: AD + Pyrolysis + CHP, NPV max, county level
    10: HTC w/ selling hydrochar, NPV max, county level
    11: HTC w/ DLA, NPV max, county level
    12: HTL w/ combustion of hydrochar, NPV max, county level
    50: Pyrolysis + CHP GWP min, county level
    51: Pyroylsis + CHP Onondaga county Pareto Front min GWP max NPV
    52: Pyroylsis + CHP Jefferson county Pareto Front min GWP max NPV
    421: AD + CHP w/ disposal of digestate, NPV max, two optimal facilities
    422: Direct Land Application, NPV max, two optimal facilities
    423: Optimal County Level Results, NPV max, two optimal facilities
    424: FLP, NPV max, two optimal facilities
    425: Pyrolysis + CHP, NPV max, two optimal facilities
    426: HTL + CHP, NPV max, two optimal facilities
    427: HTC + CHP, NPV max, two optimal facilities
    428: AD + CHP, NPV max, two optimal facilities
    429: AD + Pyrolysis, NPV max, two optimal facilities
    431: AD + CHP w/ disposal of digestate, NPV max, three optimal facilities
    432: Direct Land Application, NPV max, three optimal facilities
    433: Optimal County Level Results, NPV max, three optimal facilities
    434: FLP, NPV max, three optimal facilities
    435: Pyrolysis + CHP, NPV max, three optimal facilities
    436: HTL + CHP, NPV max, three optimal facilities
    437: HTC + CHP, NPV max, three optimal facilities
    438: AD + CHP, NPV max, three optimal facilities
    439: AD + Pyrolysis, NPV max, three optimal facilities
    441: AD + CHP w/ disposal of digestate, NPV max, four optimal facilities
    442: Direct Land Application, NPV max, four optimal facilities
    443: Optimal County Level Results, NPV max, four optimal facilities
    444: FLP, NPV max, four optimal facilities
    445: Pyrolysis + CHP, NPV max, four optimal facilities
    446: HTL + CHP, NPV max, four optimal facilities
    447: HTC + CHP, NPV max, four optimal facilities
    448: AD + CHP, NPV max, four optimal facilities
    449: AD + Pyrolysis, NPV max, four optimal facilities
    1001:, CLCA,  Direct Land Application, Pareto Front NPV max GWP min, largest facility from 2 facility FLP
    1002:, CLCA,  Direct Land Application, Pareto Front NPV max GWP min, Onondaga county
    1003:, CLCA,  Direct Land Application, Pareto Front NPV max GWP min, Jefferson county
    1011:, CLCA,  Direct Land Application, Pareto Front NPV max freshwater eutrophication min, largest facility from 2 facility FLP
    1012:, CLCA,  Direct Land Application, Pareto Front NPV max freshwater eutrophication min, Onondaga county
    1013:, CLCA,  Direct Land Application, Pareto Front NPV max freshwater eutrophication min, Jefferson county
    1101:, CLCA,  Pyrolysis, Pareto Front NPV max GWP min, largest facility from 2 facility FLP
    1102:, CLCA,  Pyrolysis, Pareto Front NPV max GWP min, Onondaga county
    1103:, CLCA,  Pyrolysis, Pareto Front NPV max GWP min, Jefferson county
    1111:, CLCA,  Pyrolysis, Pareto Front NPV max freshwater eutrophication min, largest facility from 2 facility FLP
    1112:, CLCA,  Pyrolysis, Pareto Front NPV max freshwater eutrophication min, Onondaga county
    1113:, CLCA,  Pyrolysis, Pareto Front NPV max freshwater eutrophication min, Jefferson county
    1201:, CLCA,  HTL, Pareto Front NPV max GWP min, largest facility from 2 facility FLP
    1202:, CLCA,  HTL, Pareto Front NPV max GWP min, Onondaga county
    1203:, CLCA,  HTL, Pareto Front NPV max GWP min, Jefferson county
    1211:, CLCA,  HTL, Pareto Front NPV max freshwater eutrophication min, largest facility from 2 facility FLP
    1212:, CLCA,  HTL, Pareto Front NPV max freshwater eutrophication min, Onondaga county
    1213:, CLCA,  HTL, Pareto Front NPV max freshwater eutrophication min, Jefferson county
    1301:, CLCA,  HTC, Pareto Front NPV max GWP min, largest facility from 2 facility FLP
    1302:, CLCA,  HTC, Pareto Front NPV max GWP min, Onondaga county
    1303:, CLCA,  HTC, Pareto Front NPV max GWP min, Jefferson county
    1311:, CLCA,  HTC, Pareto Front NPV max freshwater eutrophication min, largest facility from 2 facility FLP
    1312:, CLCA,  HTC, Pareto Front NPV max freshwater eutrophication min, Onondaga county
    1313:, CLCA,  HTC, Pareto Front NPV max freshwater eutrophication min, Jefferson county
    1401:, CLCA,  AD, Pareto Front NPV max GWP min, largest facility from 2 facility FLP
    1402:, CLCA,  AD, Pareto Front NPV max GWP min, Onondaga county
    1403:, CLCA,  AD, Pareto Front NPV max GWP min, Jefferson county
    1411:, CLCA,  AD, Pareto Front NPV max freshwater eutrophication min, largest facility from 2 facility FLP
    1412:, CLCA,  AD, Pareto Front NPV max freshwater eutrophication min, Onondaga county
    1413:, CLCA,  AD, Pareto Front NPV max freshwater eutrophication min, Jefferson county
    1501:, CLCA,  All, Pareto Front NPV max GWP min, largest facility from 2 facility FLP
    1502:, CLCA,  All, Pareto Front NPV max GWP min, Onondaga county
    1503:, CLCA,  All, Pareto Front NPV max GWP min, Jefferson county
    1511:, CLCA,  All, Pareto Front NPV max freshwater eutrophication min, largest facility from 2 facility FLP
    1512:, CLCA,  All, Pareto Front NPV max freshwater eutrophication min, Onondaga county
    1513:, CLCA,  All, Pareto Front NPV max freshwater eutrophication min, Jefferson county
    2001:, ALCA,  Direct Land Application, Pareto Front NPV max GWP min, largest facility from 2 facility FLP
    2002:, ALCA,  Direct Land Application, Pareto Front NPV max GWP min, Onondaga county
    2003:, ALCA,  Direct Land Application, Pareto Front NPV max GWP min, Jefferson county
    2011:, ALCA,  Direct Land Application, Pareto Front NPV max freshwater eutrophication min, largest facility from 2 facility FLP
    2012:, ALCA,  Direct Land Application, Pareto Front NPV max freshwater eutrophication min, Onondaga county
    2013:, ALCA,  Direct Land Application, Pareto Front NPV max freshwater eutrophication min, Jefferson county
    2101:, ALCA,  Pyrolysis, Pareto Front NPV max GWP min, largest facility from 2 facility FLP
    2102:, ALCA,  Pyrolysis, Pareto Front NPV max GWP min, Onondaga county
    2103:, ALCA,  Pyrolysis, Pareto Front NPV max GWP min, Jefferson county
    2111:, ALCA,  Pyrolysis, Pareto Front NPV max freshwater eutrophication min, largest facility from 2 facility FLP
    2112:, ALCA,  Pyrolysis, Pareto Front NPV max freshwater eutrophication min, Onondaga county
    2113:, ALCA,  Pyrolysis, Pareto Front NPV max freshwater eutrophication min, Jefferson county
    2201:, ALCA,  HTL, Pareto Front NPV max GWP min, largest facility from 2 facility FLP
    2202:, ALCA,  HTL, Pareto Front NPV max GWP min, Onondaga county
    2203:, ALCA,  HTL, Pareto Front NPV max GWP min, Jefferson county
    2211:, ALCA,  HTL, Pareto Front NPV max freshwater eutrophication min, largest facility from 2 facility FLP
    2212:, ALCA,  HTL, Pareto Front NPV max freshwater eutrophication min, Onondaga county
    2213:, ALCA,  HTL, Pareto Front NPV max freshwater eutrophication min, Jefferson county
    2301:, ALCA,  HTC, Pareto Front NPV max GWP min, largest facility from 2 facility FLP
    2302:, ALCA,  HTC, Pareto Front NPV max GWP min, Onondaga county
    2303:, ALCA,  HTC, Pareto Front NPV max GWP min, Jefferson county
    2311:, ALCA,  HTC, Pareto Front NPV max freshwater eutrophication min, largest facility from 2 facility FLP
    2312:, ALCA,  HTC, Pareto Front NPV max freshwater eutrophication min, Onondaga county
    2313:, ALCA,  HTC, Pareto Front NPV max freshwater eutrophication min, Jefferson county
    2401:, ALCA,  AD, Pareto Front NPV max GWP min, largest facility from 2 facility FLP
    2402:, ALCA,  AD, Pareto Front NPV max GWP min, Onondaga county
    2403:, ALCA,  AD, Pareto Front NPV max GWP min, Jefferson county
    2411:, ALCA,  AD, Pareto Front NPV max freshwater eutrophication min, largest facility from 2 facility FLP
    2412:, ALCA,  AD, Pareto Front NPV max freshwater eutrophication min, Onondaga county
    2413:, ALCA,  AD, Pareto Front NPV max freshwater eutrophication min, Jefferson county
    2501:, ALCA,  All, Pareto Front NPV max GWP min, largest facility from 2 facility FLP
    2502:, ALCA,  All, Pareto Front NPV max GWP min, Onondaga county
    2503:, ALCA,  All, Pareto Front NPV max GWP min, Jefferson county
    2511:, ALCA,  All, Pareto Front NPV max freshwater eutrophication min, largest facility from 2 facility FLP
    2512:, ALCA,  All, Pareto Front NPV max freshwater eutrophication min, Onondaga county
    2513:, ALCA,  All, Pareto Front NPV max freshwater eutrophication min, Jefferson county
    3000: 1 pyrolysis plant for the whole state - max NPV
    3001: 2 pyrolysis plants for the whole state - max NPV
    3002: pyrolysis plant in every county, calculate biochar break-even prices - max NPV
    3003: plant in every county, calculate biochar break-even prices - max NPV
    3100: 1 pyrolysis plant for the whole state - max NPV and min GWP
    3101: 2 pyrolysis plants for the whole state - max NPV and min GWP
    3102: pyrolysis plant in every county, calculate biochar break-even prices - max NPV and min GWP
    3103: plant in every county, calculate biochar break-even prices - max NPV and min GWP
    3200: 1 pyrolysis plant for the whole state - min GWP
    3201: 2 pyrolysis plants for the whole state - min GWP
    3202: pyrolysis plant in every county, calculate biochar break-even prices - min GWP
    3203: plant in every county, calculate biochar break-even prices - min GWP
    10003: first calculated optimal plants at NPV max in every county, then implemented constraints on those plants for sensitivity analysis
    10103: first calculated optimal plants at a tradeoff in every county, then implemented constraints on those plants for sensitivity analysis
    10203: first calculated optimal plants at GWP min in every county, then implemented constraints on those plants for sensitivity analysis
    '''

    S = [11, 12]

    for scenario in S:
        lca_type = "CLCA"
        print("scenario", scenario)
        # scenario defines LCA type

        # clear existing/old TEA data
        files = glob.glob("data/" + str(scenario) + "/TEA/*")
        for file in files:
            os.remove(file)

        # clear existing/old LCA data
        files = glob.glob("data/" + str(scenario) + "/LCA/*")
        for file in files:
            os.remove(file)

        print("old files removed")

        if 1000 < scenario < 1999:
            lca_type = "CLCA"
        elif 2000 < scenario < 2999:
            lca_type = "ALCA"

        # scenario defines midpoint type
        if (int(scenario / 10) % 10) == 1:
            midpoint = 'eutrophication: freshwater'
        else:
            midpoint = "climate change"

        # scenario defines the number of counties
        if scenario in [51, 52, 3000, 3100, 3200] or 1000 < scenario < 2999:
            number_counties = 1
        elif scenario in [3001, 3101, 3201] or 420 < scenario < 430:
            number_counties = 2
        elif 430 < scenario < 440:
            number_counties = 3
        elif 440 < scenario < 450:
            number_counties = 4
        else:
            number_counties = 62

        # initialize the model
        for j in range(number_counties):
            initialize_model(scenario, j, midpoint, lca_type)
