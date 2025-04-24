import copy
import time

def generate_combs_from_n(n, result, cur = [], first_nonzero = False, last_nonzero = False):
     if n == 0:
          result.append(cur)
          return
     if first_nonzero == False: #last_nonzero = False also
          for x in range(3):
               new_cur = copy.deepcopy(cur)
               new_cur.append(x)
               if x != 0:
                   generate_combs_from_n(n-1, result, new_cur, True, last_nonzero)
               else:
                   generate_combs_from_n(n-1, result, new_cur, False, last_nonzero)
     elif last_nonzero == False: #first_nonzero = True
          for x in range(3):
               new_cur = copy.deepcopy(cur)
               new_cur.append(x)
               if x == 0:
                   generate_combs_from_n(n-1, result, new_cur, first_nonzero, True)
               else:
                   generate_combs_from_n(n-1, result, new_cur, first_nonzero, False)
     else:
          x = 0
          new_cur = copy.deepcopy(cur)
          new_cur.append(x)
          generate_combs_from_n(n-1, result, new_cur, first_nonzero, last_nonzero)

def compute_combs_fast(possible_combs, n_max = 20):
    for n in range(3, n_max + 1):
         print(n)
         for idx in range(2 * n):
              tup = n, idx+1
              possible_combs[tup] = []
         res = []
         generate_combs_from_n(n, res)
         for comb in res:
              m = sum(comb)
              if m == 0 or m < n - 5 or m > n + 4:
                  continue
              tup = n, m
              possible_combs[tup].append(comb)

# possible_combs = dict()
# precomputing_time = time.time()
# compute_combs_fast(possible_combs)
# print(len(possible_combs[(10,13)]))
# print("--- precomputing time: %s seconds ---" % (time.time() - precomputing_time))
# print(len(possible_combs))
