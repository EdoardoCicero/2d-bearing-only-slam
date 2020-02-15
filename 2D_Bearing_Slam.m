# SLAM 2D BEARING ONLY SLAM WITH LS OPTIMIZATION

close all
clear
clc

#load dependencies
addpath "../../"
addpath "./tools/g2o_wrapper"
addpath "./tools/visualization"
source "./tools/utilities/geometry_helpers_2d.m"

addpath "./scripts"

#load your own dataset, without landmarks (first entry remains empty)
[_, poses, transitions, observations] = loadG2o("../datasets/dataset_point.g2o");

#set initial pose at the origin - we don't know the map and neither our location
mu = [0;  #x coordinate
      0;  #y coordinate
      0]; #orientation theta (yaw angle)
printf("initial pose: [%f, %f, %f]\n", mu(1), mu(2), mu(3));

#initialize covariance: high value means high uncertainty
sigma = eye(3);

#------------------------------------------ VISUALIZATION ONLY ------------------------------------------
#initialize GUI with initial situation
figure("name", "2d_bearing_slam",    #figure title
       "numbertitle", "off"); #remove figure number
trajectory = [mu(1), mu(2)];
#------------------------------------------ VISUALIZATION ONLY ------------------------------------------

#bookkeeping: to and from mapping between robot pose (x,y, theta) and landmark indices (i)
#all mappings are initialized with invalid value -1 (meaning that the index is not mapped)
#since we do not know how many landmarks we will observe, we allocate a large enough buffer
id_to_state_map = ones(1000, 1)*-1;
state_to_id_map = ones(1000, 1)*-1;


#simulation cycle: for the number of transitions recorded in the dataset
for t = 1:length(transitions)

  #obtain current transition
  transition = transitions(t);
  
  #obtain current observations
  observations_t = observations(t);















