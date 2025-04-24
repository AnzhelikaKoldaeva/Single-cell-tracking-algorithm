import itertools as it
import copy

def precompute_combs(possible_combs, n_max = 17):
    max_divs = 3
    for n in range(1, n_max + 1):
        res = [[] for _ in range( n + max_divs)]
        for x in it.product(range(3), repeat = n):
            m = sum(x)
            if m == 0 or m > n + max_divs:
                continue
            should_append = True
            first_nonzero = False
            last_nonzero = False
            for i in range(0, len(x)):
                if x[i] == 0 and first_nonzero:
                    last_nonzero = True
                if x[i] != 0:
                    first_nonzero = True
                if x[i] != 0 and first_nonzero and last_nonzero:
                    should_append = False
            if should_append:
                res[m - 1].append(x)
        for m in range(n + max_divs):
            print(n, m+1, len(res[m]))
            tup = n, m+1
            possible_combs[tup] = res[m]

def find_all_comb(n, m, possible_combs):
    tup = n, m
    #return if we already computed
    if tup in possible_combs:
        return possible_combs[tup]

    res = []
    for x in it.product(range(3), repeat = n):
        should_append =  sum(x) == m
        first_nonzero = False
        last_nonzero = False
        for i in range(0, len(x)):
            if x[i] == 0 and first_nonzero:
                last_nonzero = True
            if x[i] != 0:
                first_nonzero = True
            if x[i] != 0 and first_nonzero and last_nonzero:
                should_append = False
        if should_append:
            res.append(x)
    tup = n, m
    possible_combs[tup] = res
    return possible_combs[tup]



def feature_diff(f1, f2):
    return sum([
        abs(x - y) for x, y in zip(f1, f2)
                     ])

def average_feature(feature_list):
    #эту функцию ты, возможно захочешь модифицировать, так как
    #так как свойства усредняются по разному и ты можешь захотеть добавить свойство
    res = []
    #assuming you have two features
    lengths = 0.0
    orient = 0.0
    for cell in feature_list:
        lengths += cell[0]
        orient += cell[1]

    res.append(lengths)
    res.append(orient / len(feature_list))

    return res


def total_diff(prev_step, cur_step):
    #prev step is a bunch of features, i.e. if there are 2 cells in the prev step
    #prev_step = [[1, 0], [2, 0.5]],
    # where [1, 0] and [2, 0.5] are features for 1st and 2nd cells
    #in the previous row
    #cur_step is the feature list on the next step, i.e. there are 3 cells in the cur step
    #cur_step = [[[0.4, -1], [0.5, 1.1]], [[1.9, 0.6]]]
    # here 1st in prev_step corresponds to 1st and 2nd in cur_step and
    # 2nd in prev_step corresponds to 3rd in cur_step
    total = 0
    for i in range(len(prev_step)):
        if cur_step[i]:
            total += feature_diff(average_feature(cur_step[i]), prev_step[i])

    return total


def cur_step_from_comb(c, b, comb):
    cur_step = [[] for _ in range(len(b))]
    c_ = copy.deepcopy(c)
    #find cur_step corresponding to comb
    for idx, num_ancstrs in enumerate(comb):
        for num in range(num_ancstrs):
            feature = c_.pop(0)
            cur_step[idx].append(feature)
    return cur_step


def find_best_matching(b, c, possible_combs):
    n = len(b)
    m = len(c)
    combs = find_all_comb(n, m, possible_combs)
    best_score = 10**7
    best_comb = None

    for comb in combs:
        cur_step = cur_step_from_comb(c, b, comb)
        score = total_diff(b, cur_step)
        if score < best_score:
            best_score = score
            best_comb = comb

    idx = 1
    res = []
    for children in best_comb:
        for child in range(children):
            res.append(idx)
        idx += 1
    return res

