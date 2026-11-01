import numpy as np
import random
import math
import matplotlib.pyplot as plt
from numpy import arange
from numpy import exp
from numpy import sqrt
from numpy import cos
from numpy import pi
from numpy import sin
from numpy import absolute
from numpy import e
from numpy.linalg import inv
from sklearn.preprocessing import MinMaxScaler
from matplotlib.animation import FuncAnimation
import time

# two dimensional functions
"""This function returns the solution of a test function (Dependent on the variable parameter) using the x and y parameter values """


def fitness_func2d(x: float, y: float, function_name: str):
    a = (
        -20.0 * exp(-0.2 * sqrt(0.5 * (x**2 + y**2)))
        - exp(0.5 * (cos(2 * pi * x) + cos(2 * pi * y)))
        + e
        + 20
    )  # [-5,5] # ackley function
    if function_name == "Ackley":
        return a

    h = -absolute(
        sin(x) * cos(y) * exp(absolute(1 - (sqrt(x**2 + y**2) / pi)))
    )  # holder [-10,10]
    if function_name == "Holders":
        return h

    r = (
        100 * (y - x**2) ** 2 + (x - 1) ** 2
    )  # rosenbrock [-5,10] #https://www.indusmic.com/post/rosenbrock-function
    if function_name == "Rosenbrock":
        return r

    b = (x + (2 * y) - 7) ** 2 + ((2 * x) + y - 5) ** 2  # booth function [-10,10]
    # https://machinelearningmastery.com/2d-test-functions-for-function-optimization/
    if function_name == "Booth":
        return b
    easom = (
        -cos(x) * cos(y) * exp(-((x - pi) ** 2 + (y - pi) ** 2))
    )  # easoms function [-100,100] global optimium= -1 at (pi,pi).

    if function_name == "Easom":
        return easom

    himmel = (x**2 + y - 11) ** 2 + (
        x + y**2 - 7
    ) ** 2  # Himmelblau’s function [-5,5] gbest=0 at (3-0,2.0)
    if function_name == "Himmelblau":
        return himmel

    goldstein = (
        (
            1
            + (x + y + 1) ** 2
            * (19 - 14 * x + 3 * x**2 - 14 * y + 6 * x * y + 3 * y**2)
        )
        * (
            30
            + (2 * x - 3 * y) ** 2
            * (18 - 32 * x + 12 * x**2 + 48 * y - 36 * x * y + 27 * y**2)
        )
    )  # [-2,2]  global optimum = 3 at (0,-1)#https://gist.github.com/MiguelAngelHFlores/777062e58419e1458a1c1800d00b03d5
    if function_name == "Goldstein":
        return goldstein

    # for shaffer function
    part1 = np.square(x) - np.square(y)
    part2 = np.square(x) + np.square(y)
    s = 0.5 + (np.square(np.sin(part1)) - 0.5) / np.square(1 + 0.001 * part2)  # shaffer

    if function_name == "Schaffer":
        return s


# one dimensional functions
"""This function returns the solution of a one dimensional test function (Dependent on the variable parameter) using the x parameter values"""


def fitness_func1d(x: float, function_name: str) -> float:
    s = np.sum(np.square(x))
    r = np.sum([np.square(x) - 10 * np.cos(2 * np.pi * x) + 10])  # rastgrin

    if function_name == "Sphere":
        return s

    if function_name == "Rastrigin":
        return r


"""This function returns the euclidean distance between points x1 and x2 """


def eucldist(
    x1: float, x2: float
) -> float:  # calculates distance between bacteria in algorithm
    distance = np.sqrt(np.sum((x1 - x2) ** 2))
    return distance


"""This is the main ABCO algorithm"""


def BCO(
    pop_size: int,
    iteration: int,
    step_size: int,
    step_explore: int,
    step_exploit: int,
    step_tumble: int,
    dimension: int,
    lower_bound: int,
    upper_bound: int,
    threshold: float,
    pop_split: float,
    k_neighbours: int,
    weight_mode: str,
    function_name: str,
    function_x: float,
    function_y: float,
    unchangedthreshold: int,
) -> None:
    """ABCO algorithm

    Args:
        pop_size (int): size of population.
        iteration (int): number of iterations
        step_size (int): step size of bacteria
        step_explore (int): number of steps for exploration
        step_exploit (int): number of steps for exploitation
        step_tumble (int): number of bacteria tumbling
        dimension (int): number of dimensions. proportional to number of test function arguments 2 dimensions (x and y) or 1 dimension (x)
        lower_bound (int): lower bound of search space
        upper_bound (int): upper bound of search space
        threshold (float): fitness threshold. e in paper.
        pop_split (float): percentage of surviving bacteria in reproduction step..
        k_neighbours (int): number of neighbours
        weight_mode (str): min or max. Defines if it is a minimisation or maximisation problem
        function_name (str): name of fitness function
        function_x (float): x value of fitness function best solution
        function_y (float): y value of fitness function best solution
        unchangedthreshold (int): % of bacteria whose solutions didn’t change since previous check.

    Returns:
        _type_: None
    """

    start = time.monotonic()
    BCO.lb = lower_bound  # lower bound of search space
    BCO.ub = upper_bound  # upper bound of search space
    global bacteria_pos
    bacteria_pos = np.random.uniform(
        lower_bound, upper_bound, [pop_size, dimension]
    )  # #the initialisation of the bacterias in the search space
    cost = np.zeros(pop_size)  # initialisation of numpy array
    cost[:] = (
        [fitness_func2d(p[0], p[1], function_name) for p in bacteria_pos]
        if dimension == 2
        else [fitness_func1d(p[0], function_name) for p in bacteria_pos]
    )  # calculates the initial cost for the b bacteria in bacteria_pos(used in the loop)
    pbest = np.copy(bacteria_pos)  # holds the  personal best position of each bacteria
    pbest_cost = np.copy(
        cost
    )  # holds the personal pbest fitness value for all bacteria
    global count
    count = 0
    global unchanged
    global prevBacteria
    prevBacteria = np.copy(bacteria_pos)
    global step_change  # tracks how many times step_size changes
    step_change = 0
    generation_gap = iteration * 0.25  # 25 percent of iteration
    print(generation_gap)
    # threshold to check remaing unchanged bacteria. This is a percent, ie 80 percent of bacteria is unchanged
    print(unchangedthreshold)

    def explore() -> None:
        for b in range(pop_size):
            tumble = step_size * ((np.random.rand(1, dimension) - 0.5) * 2)
            bacteria_pos[b] = bacteria_pos[b] + tumble
            """for if bacterium goes outside of lb or ub"""
            if dimension == 2:
                if bacteria_pos[b][0] < lower_bound:
                    bacteria_pos[b][0] = np.random.randint(lower_bound, upper_bound)
                if bacteria_pos[b][0] > upper_bound:
                    bacteria_pos[b][0] = np.random.randint(lower_bound, upper_bound)
                if bacteria_pos[b][1] > upper_bound:
                    bacteria_pos[b][1] = np.random.randint(lower_bound, upper_bound)
                if bacteria_pos[b][1] < lower_bound:
                    bacteria_pos[b][1] = np.random.randint(lower_bound, upper_bound)
            elif dimension == 1:
                if bacteria_pos[b] < lower_bound:
                    bacteria_pos[b] = np.random.randint(lower_bound, upper_bound)
                if bacteria_pos[b] > upper_bound:
                    bacteria_pos[b] = np.random.randint(lower_bound, upper_bound)
            cost[b] = (
                fitness_func2d(bacteria_pos[b][0], bacteria_pos[b][1], function_name)
                if dimension == 2
                else fitness_func1d(bacteria_pos[b][0], function_name)
            )
            fitness = cost[b]  # get the fitness value of b bacteria in tumble step.
            difference = (
                pbest_cost[b] - fitness
            )  # difference between current cost and pbest cost
            dist = eucldist(bacteria_pos, pbest)  # distance between current b and pbest
            # check if current bacteria solution is better than global best found
            if weight_mode == "min" and fitness < pbest_cost[b]:
                pbest_cost[b] = (
                    fitness  # change the pbest to that of the fitness variable
                )
                pbest[b] = bacteria_pos[
                    b
                ]  # change the pbest position to that of the current bacteria
                if (
                    difference > threshold
                ):  # check if the difference between current bacteria solution and gbest is greater than e
                    newmv = step_size * (dist)
                    bacteria_pos[b] = (
                        bacteria_pos[b] + newmv
                    )  # move towords best solution
                else:
                    bacteria_pos[b] = (
                        bacteria_pos[b] + tumble
                    )  # move in random direction
            # check if current bacteria solution is better than global best found
            if weight_mode == "max" and fitness > pbest_cost[b]:
                pbest_cost[b] = (
                    fitness  # change the pbest to that of the fitness variable
                )
                pbest[b] = bacteria_pos[
                    b
                ]  # change the pbest position to that of the current bacteria
                if (
                    difference > threshold
                ):  # check if the difference between current bacteria solution and gbest is greater than e
                    newmv = step_size * (dist)
                    bacteria_pos[b] = (
                        bacteria_pos[b] + newmv
                    )  # move towords best solution
                else:
                    bacteria_pos[b] = (
                        bacteria_pos[b] + tumble
                    )  # move in random direction

    def exploit() -> None:
        for b in bacteria_pos:
            # calculate distance of bacteria b against other bacteria
            distance = [eucldist(b, a) for a in bacteria_pos]
            # get the closest k neighbour
            n_index = np.argsort(distance)[
                :k_neighbours
            ]  # get index of all the k neighbours of bacteria b
            cost_n = [
                pbest_cost[n] for n in n_index
            ]  # get the pbest of the k neighbours
            if weight_mode == "min":
                b_index = np.argmin(cost_n)  # get the index of the best solution
                movement = step_size * (eucldist(b, pbest[b_index]) / step_exploit)
                b = (
                    b + movement
                )  # update postition to go towards best solution amoung neighbours
            elif weight_mode == "max":
                b_index = np.argmax(cost_n)  # get the index of the best solution
                movement = step_size * (eucldist(b, pbest[b_index]) / step_exploit)
                b = (
                    b + movement
                )  # update postition to go towards best solution amoung neighbours

    def reproduction() -> None:
        # reproduction
        pcost_n = np.array([n for n in pbest_cost])  # get the pbest of the k neighbours
        if weight_mode == "min":
            sorted_index = np.argsort(pcost_n)  # get the index of the best solution
        elif weight_mode == "max":
            sorted_index = np.argsort(pcost_n)[::-1]

        bestofs = math.floor(
            len(pcost_n) * pop_split
        )  # get s best of the population. s is a float. use math.floor to get nearest integer
        # new population the same size all the time.(ie, s=60% or s=30% doesn't alter population size.)
        new_Sbacteria = np.array(
            [pbest[j] for j in sorted_index[:bestofs]]
        )  # holds the top S bacteria.
        newBacteria = []  # initialise lower s array variable
        i = pop_size - bestofs
        scaler = MinMaxScaler()  # for min-max scaling
        for b in sorted_index[:i]:  # go through the top pop_size - S  bacteria.
            distance = [
                eucldist(pbest[b], a) for a in pbest
            ]  # distance between low S and Top S bacteria.

            # get the closest k neighbour
            n_index = np.argsort(distance)[
                :k_neighbours
            ]  # get index of all the k neighbours of bacteria b
            nearest_n = np.array(
                [pbest[i] for i in n_index]
            )  # put neighbours in new array using their postions from n_index
            weight_prep = np.array([pbest_cost[c] for c in n_index]).reshape(
                (-1, 1)
            )  # reshape array to [[w1],[w2]]
            weights = scaler.fit_transform(weight_prep)

            if weight_mode == "min":
                weightn = weights + 1  # add 1 to scaled points
                w = np.array(
                    [1 / i for i in weightn]
                )  # np.invert(weightn).astype(np.int)
                n = np.array([w[i] * nearest_n[i] for i in range(k_neighbours)])
                newB = sum(n) / sum(w)
                newBacteria.append(newB)
            elif weight_mode == "max":
                weightn = weights + 1  # add 1 to scaled points
                n = np.array(
                    [weightn[i] * nearest_n[i] for i in range(k_neighbours)]
                )  # calculate the weigth[i]* bacteria[i]
                newB = sum(n) / sum(
                    weightn
                )  # sum of weighted calculation/sum(weight): weighted average.
                newBacteria.append(
                    newB
                )  # add new bacteria gotten from weighted average.

        newBacteria = np.array(
            newBacteria
        )  # convert list with new bacteria to numpy array.
        newPop = np.concatenate(
            (new_Sbacteria, newBacteria)
        )  # join S bacteria array and newBacteria(Generated from k nearest neighbours)

        bacteria_pos = np.copy(newPop)  # change the bacteria_pos to newPop
        BCO.actualBest = (
            fitness_func2d(function_x, function_y, function_name)
            if dimension == 2
            else fitness_func1d(function_x, function_name)
        )

        sol = sorted_index[:1]

        BCO.bestsolution = pcost_n[sol]
        BCO.bestbacteria = pbest[sol]
        BCO.error = BCO.bestsolution - BCO.actualBest

        "updated iteration"
        if (
            (o + 1) % generation_gap == 0
        ):  # since o starts from 0 to iteration count add plus 1 to properley check criterion
            print((o + 1) % generation_gap)
            global unchanged, prevBacteria
            unchanged = 0
            global oprevBacteria
            oprevBacteria = np.copy(prevBacteria)
            for i in range(
                len(bacteria_pos)
            ):  # loop through current position and compare to previous postion
                if (
                    (eucldist(bacteria_pos[i], prevBacteria[i])) < 0.001
                ):  # if bacteria has not moved over a threshold then add to unchanged
                    unchanged += 1

            if (unchanged / pop_size) * 100 > unchangedthreshold:
                end = time.monotonic()
                BCO.final_time = end - start

                return 0
            else:
                prevBacteria = np.copy(bacteria_pos)

    # Main
    for o in range(iteration):
        for i in range(step_explore):
            for x in range(step_tumble):
                explore()
        for i in range(step_exploit):
            exploit()
        reproduction()

    end = time.monotonic()
    BCO.final_time = end - start
    print(
        "final:",
        BCO.bestsolution,
        " at ",
        BCO.bestbacteria,
        "Actual best: ",
        BCO.actualBest,
        "error:",
        BCO.error,
    )

