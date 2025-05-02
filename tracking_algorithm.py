#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 11:18:50 2019
xx
@author: anzhelika-koldaeva

DESCRIPTION:
    
    Function with the tracking algorithm for bacteria in a microchannel

INPUT:
    - the .csv "contour" file with the following columns
        - NAME.id (unique name of each bacterium)
        - POSITION (index of a frame)
        - X (x-coordinate of a contour)
        - Y (y-coordinate of a contour)
    - the .csv "bacteria" file file with the following columns
        - NAME.id (unique name of each bacteria)
        - LOCATION.x (x-coordinate of the center)
        - LOCATION.y (y-coordinate of the center)
        - SHAPE.length (length of a bacterium)
        - SHAPE.width (width of a bacterium)
        - SHAPE.orientation (orientation of a bacterium)
OUTPUT:
    - width per frame
    - length per frame
    - orientation per frame
    - centers x per frame
    - centers y per frame
    - contours per frame
    - tracking names
    - traces
    - divisions ((x,y,z) coordinates of the mother cells)
    - max_y_frames
    - min_y_frames
    - numb_uniq_exper
    
"""

import numpy as np

from utils.read_files import read_contour_csv, read_bact_csv, find_contour_per_frame, rotate_contour_per_frame, min_max_contourper_frame, stats_per_frame, bact_poles
from utils.read_files import angle_three_points, inside_circle, mitosis_check, order_matched, euclid_dist_diff_frames, matched_edges, range_from_the_middle, DL_two, ME_edq, MC_two_str
from utils.read_files import best_of_two_candidates, end_cells, inside_rect, poles_per_frame

def tracking(nm_fl, dim):   #input files: "contour_nm_fl.csv" and "bacteria_nm_fl.csv"; input value dim = number of rows
    print("Reading the files...")

    file_name_con = 'input/contour_'+nm_fl+'.csv'
    file_name_bac = 'input/bacteria_'+nm_fl+'.csv'

    name_cont, posit_cont, x_cont, y_cont = read_contour_csv(file_name_con)

    name_bact, cent_bact_x, cent_bact_y, length_bact, width_bact, orient_bact = read_bact_csv(file_name_bac)
      
    print("Starting analysis...")
        

    contour_per_frame_centr, name_contour_per_frame = find_contour_per_frame(x_cont, y_cont, name_cont, posit_cont)
        
    contour_per_frame_centr = rotate_contour_per_frame(contour_per_frame_centr, nm_fl)

    min_max_contourper_frame(contour_per_frame_centr, nm_fl)
             
    name_centers_per_frame, centers_x_per_frame, centers_y_per_frame, orient_per_frame, width_per_frame, length_per_frame = stats_per_frame(name_contour_per_frame, contour_per_frame_centr, name_bact, length_bact, width_bact, orient_bact, nm_fl)

    left_pole_bact_per_frame_x, left_pole_bact_per_frame_y, right_pole_bact_per_frame_x, right_pole_bact_per_frame_y, left_pole_bact_per_frame_x_ori, left_pole_bact_per_frame_y_ori, right_pole_bact_per_frame_x_ori, right_pole_bact_per_frame_y_ori  =  bact_poles(contour_per_frame_centr, name_centers_per_frame, name_contour_per_frame)
    
    std_left_pole_x_glob, std_left_pole_y_glob, std_right_pole_x_glob, std_right_pole_y_glob = poles_per_frame(centers_x_per_frame, centers_y_per_frame, width_per_frame, orient_per_frame, left_pole_bact_per_frame_x, left_pole_bact_per_frame_y, right_pole_bact_per_frame_x, right_pole_bact_per_frame_y)
    
    print("Calculating local structures...")    

    print(f"number of cells per frames = {[len(name_centers_per_frame[j]) for j in range(len(name_centers_per_frame))]}")

    local_str = []
    for j in range(len(name_centers_per_frame)):
       # print("local structure for j = "+str(j))
        local_str.append([])
        for q in range(len(name_centers_per_frame[j])):
           # print(f"q = {q}")
            local_str[-1].append([])
           # print(f"len(name_centers_per_frame[j]) = {len(name_centers_per_frame[j])}")
            for q1 in range(len(name_centers_per_frame[j])):
                    for q2 in range(q1,len(name_centers_per_frame[j])):
                            if((q-q1)*(q-q2)*(q1-q2) !=0 and angle_three_points(j, q, q1, q2, centers_x_per_frame, centers_y_per_frame) < 170 and angle_three_points(j, q, q1, q2, centers_x_per_frame, centers_y_per_frame) > 2):
                                count = 0
                                for t in range(len(name_centers_per_frame[j])):
                                    if((t-q)*(t-q1)*(t-q2) != 0):
                                        if(inside_circle(j,q,q1,q2,t, centers_x_per_frame, centers_y_per_frame) == "not inside"):
                                            count = count+1
                                if(count == len(name_centers_per_frame[j])-3):
                                    local_str[j][-1].append((q,q1,q2))      

   # print(f"local_str = {local_str}")                                                         
        

    if dim != 1:
        #6. Main loop      
        print("Starting tracking...")
        
        alpha = 10  #pixels
        alpha_loop = 15 #pixels
        tracking_names = []     
        R = []
        R_global = []
        removed_global = []
        numb_uniq_exper = []
        tracking_names.append(np.ndarray.tolist(name_centers_per_frame[0]))
        mitos_global = []
        for j in range(0,len(name_centers_per_frame)-1): 
            R=[]
            tracking_names.append([0]*len(name_centers_per_frame[j+1]))
            alpha_loop = 10
            while(R == []):
                for q in range_from_the_middle(j, centers_x_per_frame):
                    MC = [0]*len(local_str[j+1])
                    for s in range_from_the_middle(j+1, centers_x_per_frame):
                        MC[s] = MC_two_str(j,q,s,alpha_loop, local_str, centers_x_per_frame, centers_y_per_frame, left_pole_bact_per_frame_y, right_pole_bact_per_frame_y)
                    if (MC.count(1) == 1):
                        ind = MC.index(1)
                        R.append((q,ind))   # 1st from j, 2nd from (j+1)
                    elif (MC.count(1) >=2):
                        print("Error! More than 2 candidates")
                        ind = best_of_two_candidates(j,q,MC, std_left_pole_x_glob, std_left_pole_y_glob, std_right_pole_x_glob, std_right_pole_y_glob, left_pole_bact_per_frame_x_ori, left_pole_bact_per_frame_y_ori, right_pole_bact_per_frame_x_ori, right_pole_bact_per_frame_y_ori)
                        R.append((q,ind)) 
                    else:
                        print("No robust candidates")
                alpha_loop = alpha_loop+0.1
        
            for r in range(len(R)):
                tracking_names[j+1][R[r][1]] = tracking_names[j][R[r][0]]
                
            removed = []
            if (min(left_pole_bact_per_frame_x[j+1]) - min(left_pole_bact_per_frame_x[j]) > 15):
                ind = left_pole_bact_per_frame_x[j].index(min(left_pole_bact_per_frame_x[j]))
                removed.append(ind)
                print("cell "+str(tracking_names[j][ind])+" is removed")
            if (max(right_pole_bact_per_frame_x[j]) - max(right_pole_bact_per_frame_x[j+1]) > 15):
                ind = right_pole_bact_per_frame_x[j].index(max(right_pole_bact_per_frame_x[j]))
                removed.append(ind)
                print("cell "+str(tracking_names[j][ind])+" is removed")
            error  =[]
            nodes_matched_j = []
            nodes_matched_j_pl_1 = []
            for r in range(len(R)):
                nodes_matched_j.append(R[r][0])
                nodes_matched_j_pl_1.append(R[r][1])
            nodes_not_matched_j = list(set(list(range(len(local_str[j])))) - set(nodes_matched_j) - set(removed))
            nodes_not_matched_j_pl_1 = list(set(list(range(len(local_str[j+1])))) - set(nodes_matched_j_pl_1))
            iterat = 0
            while(nodes_not_matched_j != [] and nodes_not_matched_j_pl_1 != [] and iterat<100):
                iterat = iterat+1
                q = order_matched(j,nodes_not_matched_j,nodes_matched_j,nodes_not_matched_j_pl_1, nodes_matched_j_pl_1, alpha, alpha_loop, centers_x_per_frame, centers_y_per_frame, left_pole_bact_per_frame_y, right_pole_bact_per_frame_y, local_str, R)
                mch_edg_j = matched_edges(j,nodes_not_matched_j[q],nodes_matched_j, local_str)
                if (mch_edg_j != []):
                    PMC = []
                    for s in range(len(nodes_not_matched_j_pl_1)):
                        mch_edg_j_pl_1 = matched_edges(j+1,nodes_not_matched_j_pl_1[s],nodes_matched_j_pl_1, local_str)
                        ME = ME_edq(mch_edg_j,mch_edg_j_pl_1, R)
                        MC = MC_two_str(j,nodes_not_matched_j[q],nodes_not_matched_j_pl_1[s],alpha_loop, local_str, centers_x_per_frame, centers_y_per_frame, left_pole_bact_per_frame_y, right_pole_bact_per_frame_y)
                        PMC.append(MC + ME/max(len(mch_edg_j), len(mch_edg_j_pl_1))+ 1 - euclid_dist_diff_frames(j,nodes_not_matched_j[q],nodes_not_matched_j_pl_1[s], centers_x_per_frame, centers_y_per_frame, left_pole_bact_per_frame_y, right_pole_bact_per_frame_y)/alpha_loop)   #(1 - euclid_dist_diff_frames(j,nodes_not_matched_j[q],nodes_not_matched_j_min_1[s])/alpha))
                    ind = PMC.index(max(PMC))
                    ind_ends = end_cells(j,3, centers_x_per_frame, centers_y_per_frame)
                    if(abs(max(PMC))>3 and nodes_not_matched_j[q] in ind_ends):
                        removed.append(nodes_not_matched_j[q])
                        print("cell "+str(tracking_names[j][q])+" is removed")
                    elif(abs(max(PMC))>3):
                        print('Error')
                        error.append(nodes_not_matched_j[q])
                    else:
                        R.append((nodes_not_matched_j[q],nodes_not_matched_j_pl_1[ind]))
                        print('for ' + str((nodes_not_matched_j[q],nodes_not_matched_j_pl_1[PMC.index(max(PMC))]))+' PMC '+str(max(PMC)))
                    
                nodes_matched_j = []
                nodes_matched_j_pl_1 = []
                for r in range(len(R)):
                    nodes_matched_j.append(R[r][0])
                    nodes_matched_j_pl_1.append(R[r][1])
                nodes_not_matched_j = list(set(list(range(len(local_str[j])))) - set(nodes_matched_j) - set(removed)- set(error))
                nodes_not_matched_j_pl_1 = list(set(list(range(len(local_str[j+1])))) - set(nodes_matched_j_pl_1))
            for r in range(len(R)):
                tracking_names[j+1][R[r][1]] = tracking_names[j][R[r][0]] 
            
            nodes_matched_j = []
            nodes_matched_j_pl_1 = []
            for r in range(len(R)):
                nodes_matched_j.append(R[r][0])
                nodes_matched_j_pl_1.append(R[r][1])        
            
            unassign_ind_j_pl_1 = []
            for u in range(len(tracking_names[j+1])):
                if (tracking_names[j+1][u] == 0):
                    unassign_ind_j_pl_1.append(u)
            
            order = []  
            rg_md = range_from_the_middle(j+1, centers_x_per_frame)     
            for i in range(len(unassign_ind_j_pl_1)):
                order.append(rg_md.index(unassign_ind_j_pl_1[i]))
            sr = sorted((e,i) for i,e in enumerate(order))
            
            unassign_ind_j_pl_1_sorted = []
            for s in sr:
                unassign_ind_j_pl_1_sorted.append(unassign_ind_j_pl_1[s[1]])
                
            for un_j_pl_1 in unassign_ind_j_pl_1_sorted:
                candidates_un_j = []
                thresh = 15
                while(candidates_un_j == [] and thresh < 40):
                    for q in range(len(local_str[j])):
                        x = centers_x_per_frame[j][q]-centers_x_per_frame[j+1][un_j_pl_1]
                        y = centers_y_per_frame[j][q]-centers_y_per_frame[j+1][un_j_pl_1]
                        orient = orient_per_frame[j][q]
                        if(inside_rect(x,y,orient,thresh) == "inside"):
                            candidates_un_j.append(q)
                            print("inside")
                        else:
                            print("not inside; threshold = "+str(thresh))
  
                    candidates_un_j = list(set(candidates_un_j)- set(removed))
                    thresh = thresh+1
                if(candidates_un_j != []):
                    DL = []
                    for cand_j in candidates_un_j:
                        DL.append(DL_two(j,un_j_pl_1,cand_j,centers_x_per_frame,centers_y_per_frame, local_str, R))
                    ind = DL.index(min(DL))
                    occ= 0
                    for i in range(len(R)):
                        if(R[i][0] == candidates_un_j[ind]):
                            occ = occ+1
                    if(occ<2):
                        R.append((candidates_un_j[ind],un_j_pl_1))
            for r in range(len(R)):
                tracking_names[j+1][R[r][1]] = tracking_names[j][R[r][0]]  
            R_global.append(R)
            removed_global.append(removed)
            numb_uniq_exper.append(len(set(tracking_names[-1])-set([0])))
            tr = 20 #pixels
            if(mitosis_check(j,R,tr, length_per_frame, centers_x_per_frame, centers_y_per_frame) != []):
                mitos_global.append(mitosis_check(j,R,tr, length_per_frame, centers_x_per_frame, centers_y_per_frame))
                
    else:     
        tracking_names = []     
        R = []
        R_global = []
        removed_global = []
        numb_uniq_exper = []
        tracking_names.append(np.ndarray.tolist(name_centers_per_frame[0]))
        mitos_global = []
        
        print('Starting precomputing...')
        
        from best_matching_v2 import find_best_matching
        from fast_combs import compute_combs_fast
        import time
        
        possible_combs = dict()
        
        precomputing_time = time.time()
        compute_combs_fast(possible_combs)
        print("--- precomputing time: %s seconds ---" % (time.time() - precomputing_time))

       # print(f"len(name_centers_per_frame) = {len(name_centers_per_frame)}")

       # print(f"possible_combs = {possible_combs}")

        
        # Iterating over all frames
        for j in range(len(name_centers_per_frame)-1):
            print(f"len(name_centers_per_frame[j]) = {len(name_centers_per_frame[j])}")
            R = []
            tracking_names.append([0]*len(name_centers_per_frame[j+1]))
            
            a_ind_orig = list(range(len(name_centers_per_frame[j])))
            a_ind_sorted = [x for _,x in sorted(zip(centers_x_per_frame[j],a_ind_orig))]
            
            b_ind_orig = list(range(len(name_centers_per_frame[j+1])))
            b_ind_sorted = [x for _,x in sorted(zip(centers_x_per_frame[j+1],b_ind_orig))]
            
            a = [[length_per_frame[j][l], orient_per_frame[j][l]] for l in a_ind_sorted]
            b = [[length_per_frame[j+1][l], orient_per_frame[j+1][l]] for l in b_ind_sorted]

           # 
            
            best_matching = find_best_matching(a, b, possible_combs)
            
            R = [(a_ind_sorted[best_matching[i]-1],b_ind_sorted[i]) for i in range(len(best_matching))]
            
            for r in range(len(R)):
                tracking_names[j+1][R[r][1]] = tracking_names[j][R[r][0]] 

            R_global.append(R)
            numb_uniq_exper.append(len(set(tracking_names[-1])-set([0])))
            tr = 15 #pixels
            if(mitosis_check(j,R,tr, length_per_frame, centers_x_per_frame, centers_y_per_frame) != []):
                mitos_global.append(mitosis_check(j,R,tr, length_per_frame, centers_x_per_frame, centers_y_per_frame))
                
    f = open(nm_fl+"numb_uniq_exper.txt","w")
    f.write("numb_uniq_exper = "+str(numb_uniq_exper)+ '\n')
    f.close()  
    
    f = open(nm_fl+"tracking_names.txt","w")
    f.write("tracking_names = "+str(tracking_names)+ '\n')
    f.close()  
    
    traces_names = tracking_names[0]
    traces = []
    for nm in range(len(traces_names)): #names
        traces.append([])
        for j in range(len(R_global)):
            for i in range(len(R_global[j])):
                if(tracking_names[j][R_global[j][i][0]] == traces_names[nm]):
                    p1_j = R_global[j][i][0]
                    p2_j_pl_1 = R_global[j][i][1]
                    traces[nm].append(((centers_x_per_frame[j][p1_j],centers_y_per_frame[j][p1_j]),(centers_x_per_frame[j+1][p2_j_pl_1],centers_y_per_frame[j+1][p2_j_pl_1]),(j,j+1)))

    f = open(nm_fl+"traces.txt","w")
    f.write("traces = "+str(traces)+ '\n')
    f.close()  
    
    divis = []
    for j in range(len(traces)):
        traces_rec = []
        for i in range(len(traces[j])):
            traces_rec.append((traces[j][i][0],traces[j][i][2][0]))
            traces_rec.append((traces[j][i][1],traces[j][i][2][1]))
        divis.append([])
        for i1 in range(len(mitos_global)):
            for i2 in range(len(mitos_global[i1])):
                if(mitos_global[i1][i2] in traces_rec):
                    divis[j].append(mitos_global[i1][i2])
    
    numb_divis = 0
    for j in range(len(mitos_global)):
        numb_divis = numb_divis + len(mitos_global[j])
        
    f = open( nm_fl+"divis.txt","w")
    f.write("divis = "+str(divis)+ '\n'+"numb_divis = "+str(numb_divis)+ '\n')
    f.close()      

if __name__ == '__main__':
    nm_fl = "250205_075um_03_combined"
    dim = 1#"1"
    tracking(nm_fl, dim)

    