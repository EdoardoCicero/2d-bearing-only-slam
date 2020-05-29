# 2D BEARING ONLY SLAM #

# 1) Landmarks initialization with 2 or more observations from all poses
# 2) Least Squares
# 3) Dynamic Plot
# 3.1) Static Plot (if uncommented)

# P.S. if the biggest id of a landmark is greater than 300 change dimension of id_to_state_map, state_to_id_map, bearing_count

close all
clear
clc

#load dependencies
addpath "./tools/g2o_wrapper"
addpath "./tools/visualization"
source "./scripts/functions.m"
source "./tools/utilities/geometry_helpers_2d.m"

#load dataset
[_, poses, transitions, observations] = loadG2o("./dataset/slam2D_bearing_only_initial_guess.g2o");


#bookkeeping: to and from mapping between robot pose (x,y, theta) and landmark indices (i)
#all mappings are initialized with invalid value -1 (meaning that the index is not mapped)
#since we do not know how many landmarks we will observe, we allocate a large enough buffer
id_to_state_map = ones(300, 1)*-1;
state_to_id_map = ones(300, 1)*-1;


# ---- creating a structure (info) for all landmarks infos needed ----- #

n = 0;
rob_poses = [];
for t = 1:length(observations)
   
  # struct array with id and bearing of robot
  obs = observations(t).observation;
  
  # robot pose at every observation
  rob_pose = [poses(t+1).x, poses(t+1).y, poses(t+1).theta];

  # all robot poses in one variable   
  state_poses(end+1:end+3, 1) = rob_pose;
  
  # info is a struct array containing: landmark id, every landmark observation, every robot pose from which
  # the observation was performed
  for i = 1:length(obs)
    if obs(i).id ==0 # handling id = 0
       zero.id = 0;
       zero.bearing(end+1,1) = obs(i).bearing;
       zero.rob_pose(end+1, 1:3) = rob_pose;
            
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

# adding landmark with id=0, now it's id=length(info)
info(end +1).id = zero.id;
info(end).bearing = zero.bearing;
info(end).rob_pose = zero.rob_pose;
n++;
id_to_state_map(length(info)) = n;
state_to_id_map(n) = length(info);

printf('\n');
printf('Number of landmarks: %i \n\n', n)


# --------- Get landmarks position ----------- #


# theta_w = bearing angle wrt to origin ( theta_w = rob_pose.theta + land.bearing )
# (x_r,y_r) = landmark position in origin RF when it measures the angle from landmark
# y_w=tg(theta_w)(x_w - x_r) + y_r


# basically i solve a system (Ax=b -> x = pinv(A)*b)
# i use the pseudoinverse because i have many observations, and so many equations, for the same landmark
# and they do not intersect in one point (that would be the landmark position)
discarded_landmarks = [];
for l=1:length(info)
   if length(info(l).id) !=0
      
      x_r = info(l).rob_pose(:, 1); 
      y_r = info(l).rob_pose(:, 2);
      theta = info(l).rob_pose(:, 3);
      bearing_angle = info(l).bearing;
      
      theta_w = theta + bearing_angle; # bearing angle in origin RF

      A = [-tan(theta_w), ones(length(theta_w),1)];
      b = [-tan(theta_w).*x_r + y_r];

      # add landmark position to info, only for VALID landmarks (at least 2 measurements)
      if rows(A) >= 2 
         pinv(A)*b;
         info(l).land_pos = pinv(A)*b;
      else
         info(l).id;
         discarded_landmarks(end+1) = info(l).id;

      endif

   endif

endfor
# print discarded landmarks
printf('The discarded landmarks with just 1 observation are:\n');
printf('%i\n', discarded_landmarks);
printf('\n');

# --------------------------- # 


# updating mapping, considering only valid landmarks
[id_to_state_map, state_to_id_map, state_landmark] = update_map(info, id_to_state_map, state_to_id_map);
#printf('Dimension of state landmarks: %i \n', length(state_landmark));

num_poses = size(state_poses,1)/3; # total number of poses
num_landmarks = size(state_landmark,1)/2; # total number of landmarks
      
        
# XR contains all the trasformations of the origin wrt robot to simplify the calculations
[XR, XL] = w2rob(state_poses, state_landmark);
system_size = num_poses*3 + num_landmarks*2;




# -------------------------------- Least Squares  +  Levenberg-Marquardt Algorithm ------------------------------------- #

# based on algorithm from "Notes on Least-Squares and SLAM" by Giorgio Grisetti, November 1, 2015

XR_backup = XR;
XL_backup = XL;


F_hat = F(XR, XL, observations, state_to_id_map, id_to_state_map);
F_new = F_hat;

lambda = 0.0001;
noise = 0; # additive noise in observation
damping_iteration = 0; # this var controls how many times the inner loop iterates
max_damping_iter = 30;
start = 1; #control variable

printf('Calculating');


while (F_hat - F_new > 0.000000001 && damping_iteration < max_damping_iter) || start
   bearing_count = zeros(1, 300);
   F_hat = F_new;
   
   H=zeros(system_size, system_size);
   b=zeros(system_size,1);


   for c =1:length(observations)
     p_index = c; # pose index in XR
     R = XR(1:2,1:2,p_index); # rotation
     t = XR(1:2,3,p_index); # traslation
     
     for ob = observations(c).observation 

        if (ob.id) == 0
           ob.id = length(id_to_state_map);
        endif
                
        l_index = id_to_state_map(ob.id) ; # index of landmark in XL
        
        if l_index == 0 
           continue
        endif
        
        bearing_count(1,ob.id) +=1;
           
        if bearing_count(1,ob.id) >=2 # optimize only if the landmark has been seen more than once

           landmark = XL(:, l_index);

           p_i = R*landmark + t; # landmark position
           x_i = p_i(1);
           y_i = p_i(2);
         
           h = atan2(p_i(2), p_i(1)); # prediction
           z = ob.bearing + noise*rand(1,1)  ; # measurement
           
           e = h - z; # error

           Jr = 1/(x_i**2 + y_i**2) * [-y_i, x_i] * [eye(2) , [-y_i; x_i]];
           Jl = 1/(x_i**2 + y_i**2) * [-y_i, x_i] * R;
      
           Hrr=Jr'*Jr;
           Hrl=Jr'*Jl;
           Hll=Jl'*Jl;
           br=Jr'*e;
           bl=Jl'*e;
           
           pose_dim = 3;
           pose_index = 3*(p_index -1) + 1; # index of robot pose in H
           landmark_dim = 2;
           landmark_index = num_poses*3 + 2*(l_index -1) + 1; # index of landmark in H
           
           H(pose_index:pose_index+pose_dim-1,
           pose_index:pose_index+pose_dim-1)+=Hrr;
           
           H(pose_index:pose_index+pose_dim-1,
           landmark_index:landmark_index+landmark_dim-1)+=Hrl;
           
           H(landmark_index:landmark_index+landmark_dim-1,
           landmark_index:landmark_index+landmark_dim-1)+=Hll;
            
           H(landmark_index:landmark_index+landmark_dim-1,
           pose_index:pose_index+pose_dim-1)+=Hrl';
            
           b(pose_index:pose_index+pose_dim-1)+=br;
            
           b(landmark_index:landmark_index+landmark_dim-1)+=bl;

        endif
     endfor


   endfor
   
   damping_iteration = 0;
   
   
   # REPEAT #
    start_2 = 1; # control variable
    
    while (damping_iteration>0 && damping_iteration < max_damping_iter) || start_2
       
      dx = -((H + eye(system_size)*lambda) \b);

      [XR, XL] = boxPlus(XR, XL, num_poses, num_landmarks, dx);
      
      F_new = F(XR, XL, observations, state_to_id_map, id_to_state_map);
      F_hat - F_new;
      
      if F_new < F_hat
         lambda = lambda/2;
         XR_backup = XR;
         XL_backup = XL;
         damping_iteration = -1;

      else
         lambda = lambda *4;
         XR = XR_backup;
         XL = XL_backup;
         damping_iteration +=1;
         
      endif
      #damping_iteration
      
      start_2 = 0;
      endwhile
start = 0;

#printf("Error: %f\n", F(XR, XL, observations, state_to_id_map, id_to_state_map));
printf('.');
endwhile



# -------------------------- VISUALIZATION PART----------------------- #
 
dynamic_plot(XR, XL, state_to_id_map, id_to_state_map, observations);

#static_plot(XR, XL,state_to_id_map);

