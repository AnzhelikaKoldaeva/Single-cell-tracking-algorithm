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
import matplotlib.pyplot as plt 
import csv


def tracking(nm_fl, dim):   #input files: "contour_nm_fl.csv" and "bacteria_nm_fl.csv"; input value dim = number of rows
    print("Reading the files...")

    file_name_con = 'contour_'+nm_fl+'.csv'
    file_name_bac = 'bacteria_'+nm_fl+'.csv'
    
    with open(file_name_con,'rt') as f: 
        data = csv.reader(f)
        row1 = next(data)
        ind_name = row1.index('NAME.id')
        ind_posit = row1.index('POSITION')  
        ind_x_cont = row1.index('X')  
        ind_y_cont = row1.index('Y')  
        name_cont = []
        posit_cont = [] 
        x_cont = [] 
        y_cont = []
        for row in data:
              name_cont.append(row[ind_name])
              posit_cont.append(row[ind_posit])
              x_cont.append(row[ind_x_cont])
              y_cont.append(row[ind_y_cont]) 
        f.close()  
        
    with open(file_name_bac,'rt') as f:  
        data = csv.reader(f)
        row1 = next(data)
        ind_nm = row1.index('NAME.id')
        ind_cent_x = row1.index('LOCATION.x') 
        ind_cent_y = row1.index('LOCATION.y') 
        ind_len = row1.index('SHAPE.length')  
        ind_wid = row1.index('SHAPE.width') 
        ind_ori = row1.index('SHAPE.orientation')

        name_bact = []
        cent_bact_x = []
        cent_bact_y = []
        length_bact = []
        width_bact = []
        orient_bact = []
        
        for row in data:
            name_bact.append(row[ind_nm])
            cent_bact_x.append(row[ind_cent_x])
            cent_bact_y.append(row[ind_cent_y])
            length_bact.append(row[ind_len])
            width_bact.append(row[ind_wid])
            orient_bact.append(row[ind_ori])
        f.close()
        
    print("Starting analysis...")
        
    #change the bacteria names
    good_names_bact = []     
    for n in range(len(name_bact)):
        good_names_bact.append('b'+str(n))
    
    good_names_contour = []
    for n in range(len(name_cont)):
        indx = name_bact.index(name_cont[n])
        good_names_contour.append(good_names_bact[indx])
                
    name_bact = good_names_bact
    name_cont = good_names_contour
    del good_names_bact
    del good_names_contour
        
    contour_per_frame = []
    name_contour_per_frame = []
    
    for k in range(len(np.unique(posit_cont))):
        contour_per_frame.append([])
        name_contour_per_frame.append([])
        for i in range(len(posit_cont)):
            if (posit_cont[i] == str(k+1)):
                contour_per_frame[-1].append((float(x_cont[i]),float(y_cont[i])))
                name_contour_per_frame[-1].append(name_cont[i])
    contour_per_frame_centr = contour_per_frame
                
    for j in range(len(contour_per_frame)):
        sum_x = 0
        sum_y = 0 
        for i in range(len(contour_per_frame[j])):
            sum_x = sum_x + contour_per_frame[j][i][0]/len(contour_per_frame[j])
            sum_y = sum_y + contour_per_frame[j][i][1]/len(contour_per_frame[j])
        for i in range(len(contour_per_frame[j])):
            contour_per_frame_centr[j][i] = (contour_per_frame[j][i][0] - sum_x, contour_per_frame[j][i][1] - sum_y)
    del contour_per_frame
    
    #find max x and min x for each frame
    for j in range(len(contour_per_frame_centr)):
        max_x = -10000
        min_x = 10000
        for i in range(len(contour_per_frame_centr[j])):
            if (contour_per_frame_centr[j][i][0] > max_x):
                max_x = contour_per_frame_centr[j][i][0]
            elif(contour_per_frame_centr[j][i][0] < min_x):
                min_x = contour_per_frame_centr[j][i][0]
        if max_x*min_x > 0:
            leng_ch = abs(abs(max_x) - abs(min_x))
        elif max_x*min_x <= 0:
            leng_ch = abs(abs(max_x) + abs(min_x))
    
        arr_y_for_x_neg = []  
        arr_y_for_x_posit = []      
        count = 0
        for i in range(len(contour_per_frame_centr[j])):
            if contour_per_frame_centr[j][i][0] < min_x+50:
                arr_y_for_x_neg.append(contour_per_frame_centr[j][i][1])
            elif contour_per_frame_centr[j][i][0] >  max_x-50:
                arr_y_for_x_posit.append(contour_per_frame_centr[j][i][1])
    
        arr_y_for_x_neg_mn = np.mean(arr_y_for_x_neg)
        arr_y_for_x_posit_mn = np.mean(arr_y_for_x_posit)
        
        if arr_y_for_x_neg_mn*arr_y_for_x_posit_mn > 0:
            a = abs(abs(arr_y_for_x_neg_mn) - abs(arr_y_for_x_posit_mn))
        elif arr_y_for_x_neg_mn*arr_y_for_x_posit_mn <= 0:
            a = abs(abs(arr_y_for_x_neg_mn) + abs(arr_y_for_x_posit_mn))
        
        sin_alph = a/leng_ch
        cos_alph = np.sqrt(1-sin_alph**2)
    
        rotated = []
        for i in range(len(contour_per_frame_centr[j])):
            rotated.append((contour_per_frame_centr[j][i][0]*cos_alph-sin_alph*contour_per_frame_centr[j][i][1], sin_alph*contour_per_frame_centr[j][i][0]+cos_alph*contour_per_frame_centr[j][i][1]))
        contour_per_frame_centr[j] = rotated
         
    max_y = 0
    min_y = 10000
    for j in range(len(contour_per_frame_centr)):
        for i in range(len(contour_per_frame_centr[j])):
            if(contour_per_frame_centr[j][i][0]>-90 and contour_per_frame_centr[j][i][0]<90):
                if (contour_per_frame_centr[j][i][1] > max_y):
                    max_y = contour_per_frame_centr[j][i][1]
                elif(contour_per_frame_centr[j][i][1] < min_y):
                    min_y = contour_per_frame_centr[j][i][1]                   
    max_y_frames = []
    min_y_frames = []
    for j in range(len(contour_per_frame_centr)):
        max_y_loc = 0
        min_y_loc = 10000
        for i in range(len(contour_per_frame_centr[j])):
            if(contour_per_frame_centr[j][i][0]>-90 and contour_per_frame_centr[j][i][0]<90):
                if (contour_per_frame_centr[j][i][1] > max_y_loc):
                    max_y_loc = contour_per_frame_centr[j][i][1]
                elif(contour_per_frame_centr[j][i][1] < min_y_loc):
                    min_y_loc = contour_per_frame_centr[j][i][1]
        max_y_frames.append(max_y_loc)
        min_y_frames.append(min_y_loc)
            
    f = open(nm_fl+"max_y_frames.txt","w")
    f.write("max_y_frames = "+str(max_y_frames)+ '\n')
    f.close()   
    
    f = open("local_structure_modif_OUTPUT/"+nm_fl+"min_y.txt","w")
    f.write("min_y = "+str(min_y)+ '\n')
    f.close()   
    
    centers_x_per_frame = []
    centers_y_per_frame = []
    name_centers_per_frame = []
    length_per_frame = []
    width_per_frame = []
    area_per_frame = []
    orient_per_frame = []
             
    for j in range(len(name_contour_per_frame)):
        unq_nam_fr = np.unique(name_contour_per_frame[j])
        name_centers_per_frame.append(unq_nam_fr)
        centers_x_per_frame.append([])
        centers_y_per_frame.append([])
        length_per_frame.append([])
        width_per_frame.append([])
        area_per_frame.append([])
        orient_per_frame.append([])
        for nm in unq_nam_fr:
            sm_cent_x = 0
            sm_cent_y = 0
            count = 0
            for i in range(len(name_contour_per_frame[j])):
                if(name_contour_per_frame[j][i] == nm):
                    count = count+1
                    sm_cent_x = sm_cent_x + contour_per_frame_centr[j][i][0]
                    sm_cent_y = sm_cent_y + contour_per_frame_centr[j][i][1]
            centers_x_per_frame[-1].append(sm_cent_x/count)
            centers_y_per_frame[-1].append(sm_cent_y/count)
            for k in range(len(name_bact)):
                if(name_bact[k] == nm):
                    length_per_frame[-1].append(float(length_bact[k]))
                    width_per_frame[-1].append(float(width_bact[k]))
                    orient_per_frame[-1].append(float(orient_bact[k]))    
        
    f = open(nm_fl+"_centers_x_per_frame.txt","w")
    f.write("centers_x_per_frame = "+str(centers_x_per_frame)+ '\n')
    f.close()    
    
    f = open(nm_fl+"centers_y_per_frame.txt","w")
    f.write("centers_x_per_frame = "+str(centers_y_per_frame)+ '\n')
    f.close()  
    
    f = open(nm_fl+"width_per_frame.txt","w")
    f.write("width_per_frame = "+str(width_per_frame)+ '\n')
    f.close()  
    
    f = open(nm_fl+"length_per_frame.txt","w")
    f.write("length_per_frame = "+str(length_per_frame)+ '\n')
    f.close()  
    
    f = open(nm_fl+"orient_per_frame.txt","w")
    f.write("orient_per_frame = "+str(orient_per_frame)+ '\n')
    f.close()  
    
    left_pole_bact_per_frame_x = []
    left_pole_bact_per_frame_y = []
    right_pole_bact_per_frame_x = []
    right_pole_bact_per_frame_y = []
    
    for j in range(len(contour_per_frame_centr)):
        left_pole_bact_per_frame_x.append([])
        left_pole_bact_per_frame_y.append([])
        right_pole_bact_per_frame_x.append([])
        right_pole_bact_per_frame_y.append([])
        for i in range(len(name_centers_per_frame[j])):
            cont_bact_i_x = []
            cont_bact_i_y = []
            for k in range(len(name_contour_per_frame[j])):
                if(name_contour_per_frame[j][k] == name_centers_per_frame[j][i]):
                    cont_bact_i_x.append(contour_per_frame_centr[j][k][0])
                    cont_bact_i_y.append(contour_per_frame_centr[j][k][1])
            ind_max_x = cont_bact_i_x.index(max(cont_bact_i_x))
            ind_min_x = cont_bact_i_x.index(min(cont_bact_i_x))
            left_pole_bact_per_frame_x[-1].append(cont_bact_i_x[ind_min_x])
            left_pole_bact_per_frame_y[-1].append(cont_bact_i_y[ind_min_x])
            right_pole_bact_per_frame_x[-1].append(cont_bact_i_x[ind_max_x])
            right_pole_bact_per_frame_y[-1].append(cont_bact_i_y[ind_max_x])    
        
    
    f = open(nm_fl+"contour_per_frame_centr.txt","w")
    f.write("contour_per_frame_centr = "+str(contour_per_frame_centr)+ '\n')
    f.close()  
    
    left_pole_bact_per_frame_x_ori  = left_pole_bact_per_frame_x
    left_pole_bact_per_frame_y_ori = left_pole_bact_per_frame_y
    right_pole_bact_per_frame_x_ori = right_pole_bact_per_frame_x
    right_pole_bact_per_frame_y_ori = right_pole_bact_per_frame_y
    
    for j in range(len(left_pole_bact_per_frame_x)):
        for i in range(len(left_pole_bact_per_frame_x[j])):
            cent_x = (left_pole_bact_per_frame_x[j][i] + right_pole_bact_per_frame_x[j][i])/2
            cent_y = (left_pole_bact_per_frame_y[j][i] + right_pole_bact_per_frame_y[j][i])/2
            l = np.sqrt((right_pole_bact_per_frame_x[j][i] - left_pole_bact_per_frame_x[j][i])**2 + (right_pole_bact_per_frame_y[j][i] - left_pole_bact_per_frame_y[j][i])**2)/2
            left_pole_bact_per_frame_x_ori[j][i] = cent_x - l
            left_pole_bact_per_frame_y_ori[j][i] = cent_y
            right_pole_bact_per_frame_x_ori[j][i] = cent_x + l
            right_pole_bact_per_frame_y_ori[j][i] = cent_y
        
    centers_x_per_frame_global = []  
    centers_y_per_frame_global = [] 
    width_per_frame_global = []   
    orient_per_frame_global = []  
    left_pole_bact_per_frame_x_global = []  
    left_pole_bact_per_frame_y_global = []  
    right_pole_bact_per_frame_x_global = [] 
    right_pole_bact_per_frame_y_global = [] 
        
    for j in range(len(centers_x_per_frame)):
        centers_x_per_frame_global = centers_x_per_frame_global + centers_x_per_frame[j]
        centers_y_per_frame_global = centers_y_per_frame_global + centers_y_per_frame[j]
        width_per_frame_global = width_per_frame_global + width_per_frame[j]
        orient_per_frame_global = orient_per_frame_global + orient_per_frame[j]
        left_pole_bact_per_frame_x_global = left_pole_bact_per_frame_x_global + left_pole_bact_per_frame_x[j]
        left_pole_bact_per_frame_y_global = left_pole_bact_per_frame_y_global + left_pole_bact_per_frame_y[j]
        right_pole_bact_per_frame_x_global = right_pole_bact_per_frame_x_global + right_pole_bact_per_frame_x[j]
        right_pole_bact_per_frame_y_global = right_pole_bact_per_frame_y_global + right_pole_bact_per_frame_y[j]
    
    std_left_pole_x_glob = np.std(left_pole_bact_per_frame_x_global)
    std_left_pole_y_glob = np.std(left_pole_bact_per_frame_y_global)
    std_right_pole_x_glob = np.std(right_pole_bact_per_frame_x_global)
    std_right_pole_y_glob = np.std(right_pole_bact_per_frame_y_global)      
        
    
    def euclid_dist_diff_frames(j,q,p): #eucl distance between point q on (j)th frame and p on (j+1)th frame
        q1 = centers_x_per_frame[j][q]
        q2 = centers_y_per_frame[j][q]
        p1 = centers_x_per_frame[j+1][p]
        p2 = centers_y_per_frame[j+1][p]
        if (abs(q2 - p2) > 5 or abs(left_pole_bact_per_frame_y[j][q] - left_pole_bact_per_frame_y[j+1][p])>5  or abs(right_pole_bact_per_frame_y[j][q] -right_pole_bact_per_frame_y[j+1][p])>5):
            scale_y = 10
        else:
            scale_y = 1
        d = np.sqrt((q1-p1)**2+(scale_y**2)*(q2-p2)**2)
        
        return d

    
    def euclid_dist_diff_frames_resc(j,q,p,delta_x,delta_y): #eucl distance between point q on (j)th frame and p on (j+1)th frame
        q1 = centers_x_per_frame[j][q]
        q2 = centers_y_per_frame[j][q]
        p1 = centers_x_per_frame[j+1][p]-delta_x
        p2 = centers_y_per_frame[j+1][p]-delta_y
        if (abs(q2 - p2) > 5 or abs(left_pole_bact_per_frame_y[j][q] - left_pole_bact_per_frame_y[j+1][p])>5  or abs(right_pole_bact_per_frame_y[j][q] -right_pole_bact_per_frame_y[j+1][p])>5):
            scale_y = 1
        else:
            scale_y = 1
        d = np.sqrt((q1-p1)**2+scale_y*(q2-p2)**2)
        return d
    
    def Hausdorff_dist(j,D1, D2): #Hausdorff measure between triangl D1 from (j+1)th frame and D2 from jth frame
        arr = []
        min_str = []
        for p in range(3): #for D1
            arr.append([])
            for q in range(3): #for D2
                arr[-1].append(euclid_dist_diff_frames(j,D2[q],D1[p]))
            min_str.append(min(arr[p]))
        fi_D1_D2 = max(min_str)
        return fi_D1_D2
    
    def MC_two_str(j,q,s,alpha): # similarity between 2 local structures: around cell q from jth frame and cell s from (j+1)th frame
        MT = []
        for qi in range(len(local_str[j][q])):
            D2 = local_str[j][q][qi]
            for si in range(len(local_str[j+1][s])):
                D1 = local_str[j+1][s][si]
                H = Hausdorff_dist(j,D1, D2) #D1 from (j+1)th frame and D2 from jth frame
                if (H < alpha):
                    MT.append(1)
                else:
                    MT.append(0)
        MC = sum(MT)/max([len(local_str[j+1][s]), len(local_str[j][q])])
        return MC
    
    def matched_edges(j,q,nodes_matched): #find edges which belong to the local structure of a cell q and the second end = matched point from R
        match_edg = []
        for l in range(len(local_str[j][q])):
            for n in nodes_matched:
                if (local_str[j][q][l][1] == n):
                    match_edg.append((q,local_str[j][q][l][1]))
                elif (local_str[j][q][l][2] == n):
                    match_edg.append((q,local_str[j][q][l][2]))  
        return list(set(match_edg))
    
    def range_from_the_middle(j): 
        cent = np.mean(centers_x_per_frame[j])
        curr = cent
        ls = []
        for i in range(len(centers_x_per_frame[j])):
            ls.append(centers_x_per_frame[j][i])
        arr = []
        while(len(ls) > 0):
            dist = []
            for i in range(len(ls)):
                dist.append(abs(curr - ls[i]))
            ind = dist.index(min(dist))
            arr.append(centers_x_per_frame[j].index(ls[ind]))
            ls.remove(ls[ind])
        return arr
        
    def len_same_frame(j, e1, e2): #length of an edge woth nodes e1 and e2 on the same frame j
        x1 = centers_x_per_frame[j][e1]
        y1 = centers_y_per_frame[j][e1]
        x2 = centers_x_per_frame[j][e2]
        y2 = centers_y_per_frame[j][e2]
        l = np.sqrt((x1-x2)**2+(y1-y2)**2)
        return l
    
    def angle_horiz_axis(j, q1, q2): #angle between vector (q1,q2) and horiz axis on jth frame
        x1 = centers_x_per_frame[j][q2]
        y1 = centers_y_per_frame[j][q2]
        x = centers_x_per_frame[j][q1]
        y = centers_y_per_frame[j][q1]    
        x2 = x1
        y2 = y
        a = (x-x1, y-y1)
        b = (x-x2, y-y2)
        arc = np.arccos((a[0]*b[0]+a[1]*b[1])/(np.sqrt(a[0]**2+a[1]**2)*np.sqrt(b[0]**2+b[1]**2)))
        return arc*180/np.pi
    
    def ME_edq(mch_edg_j,mch_edg_j_pl_1):
        a = []
        for i in range(len(mch_edg_j)):
            a.append(mch_edg_j[i][1])
        b = []
        for j in range(len(mch_edg_j_pl_1)):
            b.append(mch_edg_j_pl_1[j][1])
        count = 0
        for ai in a:
            for bi in b:
                if((ai,bi) in R):
                    count = count+1        
        return count 
        
    def DL_two(j,u,v):  # similarity between unassigned u from (j+1) and candidate v from j
        l_u_a = []
        for m in range(len(local_str[j+1][u])):
            l_u_a.append((local_str[j+1][u][m][0],local_str[j+1][u][m][1]))
            l_u_a.append((local_str[j+1][u][m][0],local_str[j+1][u][m][2]))
        l_u_a = list(set(l_u_a))
        
        l_v_b =[]
        for m in range(len(local_str[j][v])):
            l_v_b.append((local_str[j][v][m][0],local_str[j][v][m][1]))
            l_v_b.append((local_str[j][v][m][0],local_str[j][v][m][2]))
        l_v_b = list(set(l_v_b))
        
        list_for_sum = []
        for a in range(len(l_v_b)):
            for b in range(len(l_u_a)):
                if((l_v_b[a][1],l_u_a[b][1]) in R):
                    list_for_sum.append((a,b))
        if(list_for_sum != []):            
            sm_len = 0
            ln = []
            for i in range(len(list_for_sum)):
                ln.append(len_same_frame(j,l_v_b[list_for_sum[i][0]][0],l_v_b[list_for_sum[i][0]][1]))
                ln.append(len_same_frame(j+1,l_u_a[list_for_sum[i][1]][0],l_u_a[list_for_sum[i][1]][1]))
                sm_len = sm_len + abs(len_same_frame(j,l_v_b[list_for_sum[i][0]][0],l_v_b[list_for_sum[i][0]][1])-len_same_frame(j+1,l_u_a[list_for_sum[i][1]][0],l_u_a[list_for_sum[i][1]][1]))
            sm_len_resc = sm_len/max(ln)
            
            sm_angle = 0
            angle = []
            for i in range(len(list_for_sum)):
                angle.append(angle_horiz_axis(j,l_v_b[list_for_sum[i][0]][0], l_v_b[list_for_sum[i][0]][1]))
                angle.append(angle_horiz_axis(j+1,l_u_a[list_for_sum[i][1]][0], l_u_a[list_for_sum[i][1]][1]))
                sm_angle = sm_angle + abs(angle_horiz_axis(j,l_v_b[list_for_sum[i][0]][0], l_v_b[list_for_sum[i][0]][1]) - angle_horiz_axis(j+1,l_u_a[list_for_sum[i][1]][0], l_u_a[list_for_sum[i][1]][1]))
            sm_angle_resc = sm_angle/max(angle)
            
            return sm_len_resc+sm_angle_resc
        else:
            return 1000
            
    def order_matched(j,nodes_not_matched_j,nodes_matched_j,nodes_not_matched_j_pl_1):
        cand = []
        for q in range(len(nodes_not_matched_j)):
            cand.append(len(matched_edges(j,nodes_not_matched_j[q],nodes_matched_j)))   
        
        if(cand.count(max(cand)) ==1 ):
            return cand.index(max(cand))
        
        elif(cand.count(max(cand)) > 1 ):
            ind= []
            for i in range(len(cand)):
                if (cand[i] == max(cand)):
                    ind.append(i)
            PMC_tot = []
            for q in range(len(ind)):
                mch_edg_j = matched_edges(j,nodes_not_matched_j[q],nodes_matched_j)
                if (mch_edg_j != []):
                    PMC = []
                    for s in range(len(nodes_not_matched_j_pl_1)):
                        mch_edg_j_pl_1 = matched_edges(j+1,nodes_not_matched_j_pl_1[s],nodes_matched_j_pl_1)
                        ME = ME_edq(mch_edg_j,mch_edg_j_pl_1)
                        MC = MC_two_str(j,nodes_not_matched_j[q],nodes_not_matched_j_pl_1[s],alpha)
                        PMC.append(MC + ME/max(len(mch_edg_j), len(mch_edg_j_pl_1)) - euclid_dist_diff_frames(j,nodes_not_matched_j[q],nodes_not_matched_j_pl_1[s])/alpha_loop)   #(1 - euclid_dist_diff_frames(j,nodes_not_matched_j[q],nodes_not_matched_j_min_1[s])/alpha))
                    PMC_tot.append(max(PMC))
                else:
                    PMC_tot.append(-10000)
            q_good = PMC_tot.index(max(PMC_tot))
            return ind[q_good]
    
    def angle_three_points(j, q, q1, q2): #angle between 2 vectors (q,q1) and (q,q2) on jth frame
        x1 = centers_x_per_frame[j][q1]
        y1 = centers_y_per_frame[j][q1]
        x2 = centers_x_per_frame[j][q2]
        y2 = centers_y_per_frame[j][q2]
        x = centers_x_per_frame[j][q]
        y = centers_y_per_frame[j][q]    
        a = (x-x1, y-y1)
        b = (x-x2, y-y2)
        arc = np.arccos((a[0]*b[0]+a[1]*b[1])/(np.sqrt(a[0]**2+a[1]**2)*np.sqrt(b[0]**2+b[1]**2)))
        return arc*180/np.pi
    
    def inside_circle(j,p1,p2,p3,t): # frame j: is point t lies inside circle with points p1,p2,p3
        x1 = centers_x_per_frame[j][p1]
        y1 = centers_y_per_frame[j][p1]
        x2 = centers_x_per_frame[j][p2]
        y2 = centers_y_per_frame[j][p2]
        x3 = centers_x_per_frame[j][p3]
        y3 = centers_y_per_frame[j][p3]
        
        A = x1*(y2-y3)-y1*(x2-x3)+x2*y3-x3*y2
        B = (x1**2+y1**2)*(y3-y2)+(x2**2+y2**2)*(y1-y3)+(x3**2+y3**2)*(y2-y1)
        C = (x1**2+y1**2)*(x2-x3)+(x2**2+y2**2)*(x3-x1)+(x3**2+y3**2)*(x1-x2)
        D = (x1**2+y1**2)*(x3*y2-x2*y3)+(x2**2+y2**2)*(x1*y3-x3*y1)+(x3**2+y3**2)*(x2*y1-x1*y2)
        r = np.sqrt((B**2+C**2-4*A*D)/(4*A**2))
        c1 = -B/(2*A)
        c2 = -C/(2*A)
        
        x = centers_x_per_frame[j][t]
        y = centers_y_per_frame[j][t]
        if((x-c1)**2+(y-c2)**2 < r**2):
            return "inside"
        else:
            return "not inside"
    
    def mitosis_check(j,R,tr):
        mitos_checked = []
        R_j = []
        for i in range(len(R)):
            R_j.append(R[i][0])
        
        mitos_cand_j = []
        for r in np.unique(R_j):
            if(R_j.count(r) == 2):
                mitos_cand_j.append(r)
        
        mit_event = []      
        for cd in mitos_cand_j:
            mother_len = length_per_frame[j][cd]
            daugh_len = []
            for i in range(len(R)):
                if(R[i][0] == cd):
                    daugh_len.append(length_per_frame[j+1][R[i][1]])
                    mit_event.append(R[i])
            if(abs(mother_len - daugh_len[0]-daugh_len[1]) < tr):
                mitos_checked.append(((centers_x_per_frame[j][cd],centers_y_per_frame[j][cd]),j))
        return mitos_checked
    
    if dim != 1:
        print("Calculating local structures...")    
    
        local_str = []
        for j in range(len(name_centers_per_frame)):
            print("local structure for j = "+str(j))
            local_str.append([])
            for q in range(len(name_centers_per_frame[j])):
                local_str[-1].append([])
                for q1 in range(len(name_centers_per_frame[j])):
                        for q2 in range(q1,len(name_centers_per_frame[j])):
                                if((q-q1)*(q-q2)*(q1-q2) !=0 and angle_three_points(j, q, q1, q2) < 170 and angle_three_points(j, q, q1, q2) > 2):
                                    count = 0
                                    for t in range(len(name_centers_per_frame[j])):
                                        if((t-q)*(t-q1)*(t-q2) != 0):
                                            if(inside_circle(j,q,q1,q2,t) == "not inside"):
                                                count = count+1
                                    if(count == len(name_centers_per_frame[j])-3):
                                        local_str[j][-1].append((q,q1,q2))
        
        def best_of_two_candidates(j,q,MC): #q from j
            ind = []
            std_features_j = [std_left_pole_x_glob, std_left_pole_y_glob, std_right_pole_x_glob, std_right_pole_y_glob]
            for i in range(len(MC)):
                if(MC[i] == 1):
                    ind.append(i)
            feature_q = [left_pole_bact_per_frame_x_ori[j][q], left_pole_bact_per_frame_y_ori[j][q], right_pole_bact_per_frame_x_ori[j][q], right_pole_bact_per_frame_y_ori[j][q]]        
            dist = []
            for i in range(len(ind)):
                s = ind[i]
                feature_s = [left_pole_bact_per_frame_x_ori[j+1][s], left_pole_bact_per_frame_y_ori[j+1][s], right_pole_bact_per_frame_x_ori[j+1][s], right_pole_bact_per_frame_y_ori[j+1][s]]    
                sm = 0
                for v in range(len(feature_s)):
                    sm = sm+((feature_s[v]-feature_q[v])**2)/std_features_j[v]**2
                dist.append(np.sqrt(sm))
                
            ind_best = dist.index(min(dist))
            return ind[ind_best]
                    
        def end_cells(j,n): #frame j; 2*n - number of cells
            ind_pos = []
            ind_neg = []
            mn = np.mean(centers_y_per_frame[j])
            for yi in range(len(centers_y_per_frame[j])):
                if(centers_y_per_frame[j][yi] > mn):
                    ind_pos.append(yi)
                else:
                    ind_neg.append(yi)
            above = []
            below = []
            for i in range(len(centers_x_per_frame[j])):
                above.append(centers_x_per_frame[j][i])
            for i in range(len(centers_x_per_frame[j])):
                below.append(centers_x_per_frame[j][i]) 
            val = []
            for i in range(n):
                val.append(max(above))
                val.append(min(above))
                above.remove(max(above))
                above.remove(min(above))
            ind_end = []
            for v in range(len(val)):
                ind_end.append(centers_x_per_frame[j].index(val[v]))
            return ind_end
        
        def inside_rect(x,y,orient,thr):  #
            theta = (orient-90)*np.pi/180
            R_m = np.array([[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]])
            t1 = np.array([thr,thr/4])
            t2 = np.array([thr,-thr/4])
            t3 = np.array([-thr,-thr/4])
            t4 = np.array([-thr,thr/4])
            R_t1 = R_m @ t1
            R_t2 = R_m @ t2
            R_t3 = R_m @ t3
            R_t4 = R_m @ t4
            D1 = -(R_t4[1] - R_t1[1])*x + (R_t4[0] - R_t1[0])*y -(R_t1[1]*R_t4[0]-R_t1[0]*R_t4[1]) #must be >0
            D2 = -(R_t3[1] - R_t2[1])*x + (R_t3[0] - R_t2[0])*y -(R_t2[1]*R_t3[0]-R_t2[0]*R_t3[1]) # must be <0
            D3 = -(R_t2[1] - R_t1[1])*x + (R_t2[0] - R_t1[0])*y -(R_t1[1]*R_t2[0]-R_t1[0]*R_t2[1]) #must be <0
            D4 = -(R_t3[1] - R_t4[1])*x + (R_t3[0] - R_t4[0])*y -(R_t4[1]*R_t3[0]-R_t4[0]*R_t3[1]) #must be >0
            if(D1>0 and D2<0 and D3<0 and D4>0):
                return "inside"
            else:
                return "not inside"     
             
        print("Starting tracking...")
        
        alpha = 5  #pixels
        alpha_loop = 5 #pixels
        tracking_names = []     
        R = []
        R_global = []
        removed_global = []
        numb_uniq_exper = []
        tracking_names.append(np.ndarray.tolist(name_centers_per_frame[0]))
        mitos_global = []
        for j in range(0,len(name_centers_per_frame)-1): 
            print(j)
            R=[]
            tracking_names.append([0]*len(name_centers_per_frame[j+1]))
            alpha_loop = 5
            while(R == []):
                for q in range_from_the_middle(j):
                    MC = [0]*len(local_str[j+1])
                    for s in range_from_the_middle(j+1):
                        MC[s] = MC_two_str(j,q,s,alpha_loop)
                    if (MC.count(1) == 1):
                        ind = MC.index(1)
                        R.append((q,ind))   # 1st from j, 2nd from (j+1)
                    elif (MC.count(1) >=2):
                        print("Error! More than 2 candidates")
                        ind = best_of_two_candidates(j,q,MC)
                        R.append((q,ind)) 
                    else:
                        print("No robust candidates")
                alpha_loop = alpha_loop+0.1
                print(alpha_loop)
        
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
                q = order_matched(j,nodes_not_matched_j,nodes_matched_j,nodes_not_matched_j_pl_1)
                mch_edg_j = matched_edges(j,nodes_not_matched_j[q],nodes_matched_j)
                if (mch_edg_j != []):
                    PMC = []
                    for s in range(len(nodes_not_matched_j_pl_1)):
                        mch_edg_j_pl_1 = matched_edges(j+1,nodes_not_matched_j_pl_1[s],nodes_matched_j_pl_1)
                        ME = ME_edq(mch_edg_j,mch_edg_j_pl_1)
                        MC = MC_two_str(j,nodes_not_matched_j[q],nodes_not_matched_j_pl_1[s],alpha_loop)
                        PMC.append(MC + ME/max(len(mch_edg_j), len(mch_edg_j_pl_1))+ 1 - euclid_dist_diff_frames(j,nodes_not_matched_j[q],nodes_not_matched_j_pl_1[s])/alpha_loop)   #(1 - euclid_dist_diff_frames(j,nodes_not_matched_j[q],nodes_not_matched_j_min_1[s])/alpha))
                    ind = PMC.index(max(PMC))
                    ind_ends = end_cells(j,3)
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
            rg_md = range_from_the_middle(j+1)     
            for i in range(len(unassign_ind_j_pl_1)):
                order.append(rg_md.index(unassign_ind_j_pl_1[i]))
            sr = sorted((e,i) for i,e in enumerate(order))
            
            unassign_ind_j_pl_1_sorted = []
            for s in sr:
                unassign_ind_j_pl_1_sorted.append(unassign_ind_j_pl_1[s[1]])
                
            for un_j_pl_1 in unassign_ind_j_pl_1_sorted:
                candidates_un_j = []
                thresh = 5
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
                        DL.append(DL_two(j,un_j_pl_1,cand_j))
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
            if(mitosis_check(j,R,tr) != []):
                mitos_global.append(mitosis_check(j,R,tr))
                
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
        
        for j in range(0,len(name_centers_per_frame)-1):
            print(j)
            R = []
            tracking_names.append([0]*len(name_centers_per_frame[j+1]))
            
            a_ind_orig = list(range(len(name_centers_per_frame[j])))
            a_ind_sorted = [x for _,x in sorted(zip(centers_x_per_frame[j],a_ind_orig))]
            
            b_ind_orig = list(range(len(name_centers_per_frame[j+1])))
            b_ind_sorted = [x for _,x in sorted(zip(centers_x_per_frame[j+1],b_ind_orig))]
            
            a = [[length_per_frame[j][l], orient_per_frame[j][l]] for l in a_ind_sorted]
            b = [[length_per_frame[j+1][l], orient_per_frame[j+1][l]] for l in b_ind_sorted]
            
            best_matching = find_best_matching(a, b, possible_combs)
            
            R = [(a_ind_sorted[best_matching[i]-1],b_ind_sorted[i]) for i in range(len(best_matching))]
            
            for r in range(len(R)):
                tracking_names[j+1][R[r][1]] = tracking_names[j][R[r][0]] 

            R_global.append(R)
            numb_uniq_exper.append(len(set(tracking_names[-1])-set([0])))
            tr = 15 #pixels
            if(mitosis_check(j,R,tr) != []):
                mitos_global.append(mitosis_check(j,R,tr))
                
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
    