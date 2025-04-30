import numpy as np
import math
import os
import argparse
from tqdm import tqdm
import matplotlib.pyplot as plt
import random
from scipy import special as sp

"""
How to run this script?
python3 FixatDensNonNeurt_Clusters.py --N 10 --b_s 1.01 --mut_freq 0.5 --init_config 1 --plot_flag False
"""

###non-neutral case; construct density functions for arbitrary initial distribution#########
# 0s are the mutants!

NUM_ITERATIONS = 10000
OUTPUT_DIR = 'fixation_probab_results'

# Create the output directory if it doesn't exist
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Simulate cell evolution.")
    parser.add_argument("--N", type=int, required=True,
                        help="Number of cells")
    parser.add_argument("--b_s", type=float, required=True,
                        help="b+s value (e.g., 1.01)")
    parser.add_argument("--mut_freq", type=float, required=True,
                        help="Mutation frequency (between 0 and 1)")
    parser.add_argument("--init_config", type=int, required=True,
                        choices=[1, 2, 3, 4, 5, 6],
                        help="Initial configuration: "
                             "rand(1)/end(2)/rand_restr(3)/end_restr(4)/complementary(5)/one_by_one(6)")
    parser.add_argument("--plot_flag", type=bool, required=False, default = True,
                    help="Plot flag (true = plot the results)")
    args = parser.parse_args()

    if not (0 <= args.mut_freq <= 1):
        raise ValueError("The mutation frequency must be between 0 and 1.")

    return args

def fact(x):
    """
    Returns factorial of x
    """
    return math.factorial(x)

def Binom_coef(N,k):
    """
    Returns the Binomial coefficient = N!/(k!(N-k)!)
    """
    return fact(N)/(fact(N-k)*fact(k))

def P_N_k_neutr(k, N):
    """
    Returns the fixation probability of k mutants in a neutral population of N cells 
    """
    return Binom_coef(N-1,k)/(2**(N-1))

def cell_ind_sampling(curr_list, mutants, b_s, b):
    """
    Randomly samples a cell to reproduce at one step of population dynamics. The cell is sampled with respect to it's
    reproduction rate.
    """
    # Generate a list of 0s and 1s: 0 = mutants, 1 = wild-type
    modif_list = [0*(m in mutants)+1*(not m in mutants) for m in curr_list]

    # Noramlization factor
    norm_fact = sum(modif_list)*b+(len(curr_list)-sum(modif_list))*b_s

    # Sampling weights = normalized reproduction rates
    prob = [((x==1)*b+(x==0)*b_s)/norm_fact for x in modif_list]
    
    # Randomly sample an index of a cell to reproduce
    ind = random.choices(list(range(len(curr_list))), weights=prob)
    return ind[0]

def model_non_neutr_mass_1D_simul(N, init_config, mutants, b_s, b):
    """
    Implements the main loop of the script and runs the evolution of the population for NUM_ITERATIONS.
    """
    win_count = np.zeros(N)
    transit_count_matr = []
    sum_iter = np.zeros(N)
    for j in tqdm(range(NUM_ITERATIONS), desc="Iterations..."):
        curr_list = init_config[:]  # mutants reproduce at b_s!
        trans_count = 0
        # Keep running until the number of unique cell types > 1, i.e. until the fixtion is not reached 
        while len(np.unique(curr_list)) > 1:
            # Sample a cell to reproduce
            cell_ind = cell_ind_sampling(curr_list, mutants, b_s, b)
            # Sample a direction of reproduction (left or right) with equal probabilities
            if random.randint(0, 1) == 1:  # to the left
                for i in range(cell_ind):
                    curr_list[i] = curr_list[i+1]
            else:  # to the right
                for i in range(N - cell_ind - 1):
                    curr_list[N - i - 1] = curr_list[N - i - 2]
            trans_count += 1
        transit_count_matr.append(trans_count)    
        win_count[curr_list[0]] += 1
        sum_iter[curr_list[0]] += trans_count
    
    P_win = win_count / NUM_ITERATIONS
    return P_win

def generate_initial_config(init_config_str, N, mut_freq):
    """
    Generates a list of indices of mutants in a population
    """
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
    return mutants             

def plot_fixat_probabil(P_win, P_win_neutr, mutants, init_config, b_s, filename_specs):

    plt.figure()
    plt.scatter(init_config, P_win)
    plt.scatter(mutants, P_win[mutants], color = 'r', label = f'non-neutral; b+s = {b_s}')
    plt.plot(init_config, P_win_neutr, linewidth = 1, label = 'neutral')
    plt.xlabel('k')
    plt.ylabel('fixation probabilities')
    plt.title(f"sum of mut. probab = {round(sum(P_win[mutants]),2)}; sum of neutral mut. probab ={round(sum(P_win_neutr[mutants]),2)}")
    plt.legend()

  #  plt.show()    

    plt.savefig(f"{OUTPUT_DIR}/{filename_specs}.pdf")
    plt.close()

if __name__ == '__main__':
    # 1.Parse the input arguments
    args = parse_arguments()
    N = args.N
    b_s = args.b_s
    mut_freq = args.mut_freq
    init_config_str = args.init_config
    plot_flag = args.plot_flag

    # 2. Generate the initial configuration of mutants
    mutants = generate_initial_config(init_config_str, N, mut_freq)

    # 3. Run the main loop an dfind the fixat. probabil. of cells
    b = 1
    init_config = [1]*N#list(range(N))
    for i in mutants:
        init_config[i] = 0
    
    print(init_config)
    P_win = model_non_neutr_mass_1D_simul(N, init_config, mutants, b_s, b)

    # 4. Generate the fixation probabilities of cells in the neutral population
    P_win_neutr = np.array([P_N_k_neutr(k, N) for k in init_config])

    # 5. Save the results to a txt file
    filename_specs = f'fixat_probab_positions_N{N}_p{mut_freq}_bs{b_s}_initDist{init_config_str}'
    f = open(f"{OUTPUT_DIR}/{filename_specs}.txt","a")
    f.write(f'mut_ind = {mutants} \n')
    f.write(f'P_win = {P_win} \n')
    f.close() 

    # 6. Plot the fixation probabilities if needed
    if plot_flag == True:
        plot_fixat_probabil(P_win, P_win_neutr, mutants, init_config, b_s, filename_specs)

    