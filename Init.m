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
#disp(poses);
#disp(observations);
#a = transitions.v

#poses(1).id
rob_pose = [poses(1).x; poses(1).y; poses(1).theta];
#disp(rob_pose);
printf("initial pose guess: [%f, %f, %f]\n", rob_pose(1), rob_pose(2), rob_pose(3));

#bookkeeping: to and from mapping between robot pose (x,y, theta) and landmark indices (i)
#all mappings are initialized with invalid value -1 (meaning that the index is not mapped)
#since we do not know how many landmarks we will observe, we allocate a large enough buffer
id_to_state_map = ones(1000, 1)*-1;
state_to_id_map = ones(1000, 1)*-1;

%# this is used to store previous bearing measurement or -1 if there was not
%last_maes = ones(20, 1)*-1;

curr_state = rob_pose;
#curr_state = [0;0;0];
#curr_state(end+1:length(state_to_id_map)+3) = state_to_id_map;


# 1st batch of landmarks ids
observations(1).observation.bearing;
obs = observations(1).observation
num_obs = length(obs)
for i = obs
   i.id
 endfor
#m = cell2mat (ids)
%find()

#observations(1).observation(1).id
%if any(observations(1).observation.id == observations(2).observation(3).id)
%   curr_state
%endif



#fetch the position in the state vector corresponding to the actual measurement
#n = id_to_state_map(measurement.id);

function [x_l , y_l] = landmark_position(bearing1, bearing2, v_x, v_y, rob_pose)
   
   z1 = bearing1;
   z2 = bearing2;
   
   y_l = (sin(z1)*sin(z2)) * (1/sin(z2-z1)) * (v_x - v_y * cot(z2)) - rob_pose()
   
   x_l = y_l/tg(z1)
   
   
   
end

for t = 1:length(transitions)
   
  #obtain current transition
  transition = transitions(t);
 
  if t!=1
     
     previous_pose = poses(t-1);
     
     rob_pose = poses(t);
     
     delta_pose = rob_pose - previous_pose;
     
     #retrieve info about the observed landmark
     measurements = observations(t).observation;
     
     #get bearings  
     bearings = measurements.bearing;
     
     #determine how many landmarks we have seen in this step
     M = length(observations(t));
     
     if M > 0
        # loop in previous bearings 
        for prev_obs = prev_observation.observation
           prev_id = prev_obs.id
           
           #loop in current bearings
           for curr_obs = observations(t).observation
              curr_id = curr_obs.id
              
              # if there's a match calculate landmark position
              if (prev_id == curr_id)
                 
               # bearing angle in world frame  
               world_bearing = rob_pose(3) + curr_obs.bearing
               
               
             
             endif         
           endfor
        endfor
     endif
     #poses(t)
     #obtain current observations
     #observations_t = observations(t);
     
     #input vector
     #u = transition.v

  

     
     
     if M > 0
        
        #get landmark position
        
     endif
     break
     
     
     
  
     #get previous observations (used to calculate landmark position)
     prev_observations = observations_t;  

  else
  
    #get previous observations (used to calculate landmark position)
    prev_observations = observations(t);
    #prev_pose = 
    
    u = transition.v
    
    curr_state 
  endif  
endfor




%  #if I've seen no landmarks, i do nothing
%  if (M == 0)
%    continue;
%   
%  if t != 1
%     prev_meas = observations(t-1)
%     
%     

     
  
   
%   
%   
%   if t!=1 #fai differenza h(x) - z 
      
    





%
%#set initial pose at the origin - we don't know the map and neither our location
%mu = [0;  #x coordinate
%      0;  #y coordinate
%      0]; #orientation theta (yaw angle)
%printf("initial pose: [%f, %f, %f]\n", mu(1), mu(2), mu(3));
%
%#initialize covariance: high value means high uncertainty
%sigma = eye(3);
%
%#------------------------------------------ VISUALIZATION ONLY ------------------------------------------
%#initialize GUI with initial situation
%figure("name", "2d_bearing_slam",    #figure title
%       "numbertitle", "off"); #remove figure number
%trajectory = [mu(1), mu(2)];
%#------------------------------------------ VISUALIZATION ONLY ------------------------------------------
%
%#bookkeeping: to and from mapping between robot pose (x,y, theta) and landmark indices (i)
%#all mappings are initialized with invalid value -1 (meaning that the index is not mapped)
%#since we do not know how many landmarks we will observe, we allocate a large enough buffer
%id_to_state_map = ones(100, 1)*-1;
%state_to_id_map = ones(100, 1)*-1;
%
%
%#simulation cycle: for the number of transitions recorded in the dataset
%for t = 1:length(transitions)
%
%  #obtain current transition
%  transition = transitions(t);
%  
%  #obtain current observations, bearing angles
%  observations_t = observations(t);
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
