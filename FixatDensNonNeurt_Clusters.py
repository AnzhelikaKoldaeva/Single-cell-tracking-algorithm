#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 09:53:12 2022

@author: anzhelika-koldaeva
"""

###non-neutral case; construct density functions for arbitrary initial distribution#########

import numpy as np
import math
import matplotlib.pyplot as plt
import scipy
import random
from scipy import special as sp

#file = open(f"Fixation_Selection_random_freq_{mut_freq}.txt","w")
#file.close()

def cell_ind_sampling(curr_list, mutants, b_s, b):
    modif_list = [0*(m in mutants)+1*(not m in mutants) for m in curr_list]
    norm_fact = sum(modif_list)*b+(len(curr_list)-sum(modif_list))*b_s
    prob = [((x==1)*b+(x==0)*b_s)/norm_fact for x in modif_list]
    
    ind = random.choices(list(range(len(curr_list))), weights=prob)
    return ind[0]


def model_non_neutr_mass_1D_simul(N, numb_runs, init_config, mutants, b_s, b): 

    win_count = np.zeros(N)
    transit_count_matr = []
    sum_iter = np.zeros(N)
    for j in range(numb_runs):
        curr_list = init_config[:]  #mutants reproduse at b_s!
        trans_count = 0
        while len(np.unique(curr_list)) > 1:
            cell_ind = cell_ind_sampling(curr_list, mutants, b_s, b)
            
            if random.randint(0, 1) == 1: 
                for i in range(cell_ind):
                    curr_list[i] = curr_list[i+1]
            else:                         #to the right
                for i in range(N-cell_ind-1):
                     curr_list[N-i-1] = curr_list[N-i-2]
            trans_count = trans_count+1
        print('iteration number %i' %j)
        
        transit_count_matr.append(trans_count)    
        win_count[curr_list[0]] = win_count[curr_list[0]]+1
        sum_iter[curr_list[0]] = sum_iter[curr_list[0]]+trans_count
    
    P_win = win_count/numb_runs
    
    return P_win

# 0s are the mutants!
N = int(input("Number of cells = "))
b_s = float(input("b+s = "))#0.5
mut_freq = float(input("Mut. fequency = "))#0.5
init_config_str = int(input("Initial configuration (rand(1)/end(2)/rand_restr(3)/end_restr(4))/complimentary(5)/one_by_one(6) = "))
if mut_freq > 1 or mut_freq < 0:
    raise ValueError('The frequency is not valid.')


init_list = list(range(N))
mut_number = int(mut_freq*N)
if init_config_str == 1:
    mutants = random.sample(range(N), mut_number)
elif init_config_str == 2:
    mutants = list(range(mut_number))    
elif init_config_str == 3:
    num_cells_inter = int(3*np.sqrt(N))
    num_mut_inter = int(num_cells_inter*mut_freq)
    start_int = int(N/2-num_cells_inter/2)
    end_int = int(N/2+num_cells_inter/2-1)
    mutants = random.sample(range(start_int, end_int+1), num_mut_inter)
elif init_config_str == 4:
    num_cells_inter = int(3*np.sqrt(N))
    num_mut_inter = int(num_cells_inter*mut_freq)
    start_int = int(N/2-num_cells_inter/2)
    mutants = list(range(start_int, mut_number))
elif init_config_str == 5:
    mutants1 = random.sample(range(int(N/2)), int(mut_number/2))
    mutants2 = []
    for mi in range(int(N/2)):
        if not mi in mutants1:
            mutants2.append(N-mi-1)
    mutants = mutants1+mutants2
elif init_config_str == 6:
    mutants = []
    for i in range(N):
        if i%2 == 0:
            mutants.append(i)

            
    
b = 1
init_config = list(range(N))
P_win = model_non_neutr_mass_1D_simul(N, 10000, init_config, mutants, b_s, b)

def fact(x):
    return math.factorial(x)

def Binom_coef(N,k):
    return fact(N)/(fact(N-k)*fact(k))

def P_N_k_neutr(k, N):
    return Binom_coef(N-1,k)/(2**(N-1))

P_win_neutr = np.array([P_N_k_neutr(k, N) for k in init_config])


f= open(f"Fixat_probab_positions_N{N}_p{mut_freq}_bs{b_s}_initDist{init_config_str}.txt","a")
f.write(f'mut_ind = {mutants} \n')
f.write(f'P_win = {P_win} \n')
f.close() 

plt.figure()
plt.scatter(init_config, P_win)
plt.scatter(mutants, P_win[mutants], color = 'r', label = f'non-neutral; b+s = {b_s}')
plt.plot(init_config, P_win_neutr, linewidth = 1, label = 'neutral')
plt.xlabel('k')
plt.ylabel('fix. probability')
plt.title(f"sum of mut. probab = {round(sum(P_win[mutants]),2)}; sum of corr. wild-type = = {round(sum(P_win_neutr[mutants]),2)}")
plt.legend()

plt.show()
   


    
    
    
    
    