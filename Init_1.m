# SLAM 2D BEARING ONLY SLAM WITH LS OPTIMIZATION

close all
clear
clc

#load dependencies
#addpath "../../"
addpath "./tools/g2o_wrapper"
addpath "./tools/visualization"
source "./tools/utilities/geometry_helpers_2d.m"

addpath "./scripts"


#load your own dataset, without landmarks (first entry remains empty)
[_, poses, transitions, observations] = loadG2o("./dataset/slam2D_bearing_only_initial_guess.g2o");


# initial guess
pose_structure = poses
observation_structure = observations
#a = transitions.v

#poses(1).id
rob_pose = [poses(1).x; poses(1).y; poses(1).theta];
#disp(rob_pose);
printf("initial pose guess: [%f, %f, %f]\n", rob_pose(1), rob_pose(2), rob_pose(3));

#bookkeeping: to and from mapping between robot pose (x,y, theta) and landmark indices (i)
#all mappings are initialized with invalid value -1 (meaning that the index is not mapped)
#since we do not know how many landmarks we will observe, we allocate a large enough buffer
id_to_state_map = ones(300, 1)*-1;
state_to_id_map = ones(300, 1)*-1;

%# this is used to store previous bearing measurement or -1 if there was not
%last_maes = ones(20, 1)*-1;

curr_state = rob_pose;
#curr_state = [0;0;0];
#curr_state(end+1:length(state_to_id_map)+3) = state_to_id_map;



# 1st batch of landmarks ids
%observations(1).observation.bearing;
%obs = observations(1).observation
%num_obs = length(obs)
%for i = obs
%   i.id
% endfor


#observations(1).observation(1).id
%if any(observations(1).observation.id == observations(2).observation(3).id)
%   curr_state
%endif



#fetch the position in the state vector corresponding to the actual measurement
#n = id_to_state_map(measurement.id);

# vecchia funzione per il calcolo della posizione dei landmarks
function [x_l , y_l] = landmark_position(bearing1, bearing2, v_x, v_y, rob_pose)
   
   z1 = bearing1;
   z2 = bearing2;
   
   # landmark position wrt previous robot pose
   y_l = (sin(z1)*sin(z2)) * (1/sin(z2-z1)) * (v_x - v_y * cot(z2)) - v_y;
   x_l = y_l/tan(z1) - v_x;
   
   y_l = y_l + rob_pose(2);
   x_l = x_l + rob_pose(1);
   
%   R = [ cos(z1), -sin(z1);
%          sin(z1), cos(z1)];

end

%function [x_l , y_l] = land_pos(land_info1, land_info2)
   

# ci metto l'ultima e la penultima misurazione dei landmark
prev_land_meas = ones(300, 1)*-1;
last_land_meas = ones(300, 1)*-1;

# ci metto l'ultima e la penultima posa del robot per ogni landmark osservato
prev_rob_pose = nan(300, 3);
last_rob_pose = nan(300, 3);

# -------- 2 MISURAZIONI ------ #
%for t = 1:length(observations)
%  # struct array con id e posizione del landmark
%  obs = observations(t).observation;
%  # posa del robot
%  rob_pose = [poses(t+1).x, poses(t+1).y, poses(t+1).theta];
%  
%  for i = 1:length(obs)
%     # sostituisco la misura che c'è in prev con quella che c'è in last 
%     # e salvo l'ultima misura in lastlandmeas
%     if obs(i).id ==0
%        obs(i).id =1;
%     endif
%     prev_land_meas(obs(i).id) = last_land_meas(obs(i).id);
%     last_land_meas(obs(i).id) = obs(i).bearing;
%     
%     prev_rob_pose(obs(i).id, :) = last_rob_pose(obs(i).id, :);
%     last_rob_pose(obs(i).id, :) = rob_pose(1:3);
%
%   endfor
%    
%endfor 
# ---------------------------------------------------------#
n = 0;
for t = 1:length(observations)
  # struct array con id e bearing del landmark
  obs = observations(t).observation;
  # posa del robot
  rob_pose = [poses(t+1).x, poses(t+1).y, poses(t+1).theta];
  
  # info contiene struct array con: id del landmark, ogni osservazione del landmark, ogni posa del robot da cui ha fatto l'osservazione
  for i = 1:length(obs)
    if obs(i).id ==0
       strunz.id = 0;
       strunz.bearing(end+1,1) = obs(i).bearing;
       strunz.rob_pose(end+1, 1:3) = rob_pose;
    else
       info(obs(i).id).id = obs(i).id;
       info(obs(i).id).bearing(end+1,1) = obs(i).bearing;
       info(obs(i).id).rob_pose(end+1, 1:3) = rob_pose;
       
       if id_to_state_map(obs(i).id) == -1
          n++;
          id_to_state_map(obs(i).id) = n;
          state_to_id_map(n) = obs(i).id;
       endif
    endif
   endfor


endfor

max(info.id)

info(1)
length(info)
# angolo di bearing rispetto al RF del mondo
# theta_w = rob_pose.theta + land.bearing
# voglio y=mx+q in world -> y_w=tg(theta_w)(x_w - x_r) + y_r
# con x_r, y_r posizione del robot quando osserva l'angolo


state = [0;0;0];
# costruisco la matrice di coefficienti e ricavo la posizione dei landmarks (Ax=b -> x = pinv(A)*b)
for l=1:length(info)
   if length(info(l).id) !=0
      # angolo di bearing rispetto a RF world
      x_r = info(l).rob_pose(:, 1)
      y_r = info(l).rob_pose(:, 2)
      theta = info(l).rob_pose(:, 3)
      bearing_angle = info(l).bearing;
      
      theta_w = theta + bearing_angle;
      
%      rows(theta_w), rows(bearing_angle)
%      columns(theta_w), columns(bearing_angle)
      A = [-tan(theta_w), ones(length(theta_w),1)];
      b = [-tan(theta_w).*x_r + y_r];
      
      info(l).land_pos = pinv(A)*b;
      
   endif
endfor   
count = 0
for l =1:length(state_to_id_map)
   if state_to_id_map(l) != -1
      id = state_to_id_map(l);
      count ++;
      state(end+1: end +2, 1) = info(id).land_pos;
   endif
   
endfor
   
 count  
   
   
   
   
   
   
   
   
   
   
   
   

# numero totale di landmarks
total_land = 0;

# numero di landmarks con almeno 2 misurazioni
valid_total_land = 0;


for i=1:length(prev_land_meas)
   if last_land_meas(i) != -1
      total_land ++;
   endif
   #if last_land_meas(i) != -1
   
endfor

if isnan(last_rob_pose(1,1))
  rob_pose
endif  
  
  
  
  
  
  
  
  
  
  
  