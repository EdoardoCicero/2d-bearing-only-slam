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
   

# l'ultima e la penultima misurazione dei landmark
prev_land_meas = ones(300, 1)*-1;
last_land_meas = ones(300, 1)*-1;

# l'ultima e la penultima posa del robot per ogni landmark osservato
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

function A = v2t(v)
   c = cos(v(3));
   s = sin(v(3));
   
   A = [c, -s, v(1) ;
        s, c,  v(2) ;
        0, 0,  1];
endfunction;

# ---- creazione di una struttura che contenga tutte le informazioni necessarie ----- #

n = 0;
rob_poses = [];
observations(1)
for t = 1:length(observations)
  global XR
  # struct array con id e bearing del landmark
  obs = observations(t).observation;
  # posa del robot
  rob_pose = [poses(t+1).x, poses(t+1).y, poses(t+1).theta];


     
  state_poses(end+1:end+3, 1) = rob_pose;
  # info contiene struct array con: id del landmark, ogni osservazione del landmark, ogni posa del robot da cui ha fatto l'osservazione
  for i = 1:length(obs)
    if obs(i).id ==0 # gestisco id del landmark 0
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

# aggiungo landmark con id=0
info(end +1).id = zero.id;
info(end).bearing = zero.bearing;
info(end).rob_pose = zero.rob_pose;
n++;
id_to_state_map(length(info)) = n;
state_to_id_map(n) = length(info);

printf('Num landmarks: %i \n', n)
# landmark 0 cambiato con length(info)
# ------------------------------------------------# 
info;
info(1);
info(2);






# ---- Ottengo posizione dei landmarks --- #


# angolo di bearing rispetto al RF del mondo  
# theta_w = rob_pose.theta + land.bearing
# voglio y=mx+q in world -> y_w=tg(theta_w)(x_w - x_r) + y_r
# con x_r, y_r posizione del robot quando osserva l'angolo




# costruisco la matrice di coefficienti e ricavo la posizione dei landmarks (Ax=b -> x = pinv(A)*b)
# uso la pinv perchè probabilmente non esiste una posizione del landmark che soddisfi tutte le equazioni
discarded_landmarks = [];
for l=1:length(info)
   if length(info(l).id) !=0
      
      x_r = info(l).rob_pose(:, 1); 
      y_r = info(l).rob_pose(:, 2);
      theta = info(l).rob_pose(:, 3);
      bearing_angle = info(l).bearing;
      
      theta_w = theta + bearing_angle; # angolo di bearing nel RF world

      A = [-tan(theta_w), ones(length(theta_w),1)];
      b = [-tan(theta_w).*x_r + y_r];

      
      if rows(A) >= 2  # aggiungo posizione del landmark per quelli che hanno almeno 2 misurazioni di bearing angle
         pinv(A)*b;
         info(l).land_pos = pinv(A)*b;
      else
         info(l).id;
         discarded_landmarks(end+1) = info(l).id;

      endif
      
   endif

endfor
discarded_landmarks
%printf('Discarded Landmarks: %i', discarded_landmarks)

%for t =1:length(info)
%   info(t)
%endfor
length(info(46).land_pos) != 0;
# update id_to_state e state_to id considerando solo i landmarks che sono stati inizializzati (minimo 2 measurements)
info(46);
info(46).land_pos;
state_to_id_map;
function [id_to_state, state_to_id, state_landmark] =  update_map( info, id_to_state_map, state_to_id_map)
   #global info
   #global XL
   global land_count
   global state_landmark
%   state_to_id_map_new = ones(300)*-1
%   id_to_state_map_new = ones(300)*-1
   
%   state = curr_state; 
   land_count = 0;
   #id_to_state = ones(300, 1)*-1;
   #state_to_id = ones(300, 1)*-1;
   for l =1:length(state_to_id_map)
      id = state_to_id_map(l);


      if id != -1 && size(info(id).land_pos,1) !=0        
            land_count ++;
            state_to_id_map_new(land_count) = id;
            id_to_state_map_new(id) = land_count;
            info(id).land_pos;
            #XL(:,land_count) = info(id).land_pos;
            
            state_landmark(end+1: end +2, 1) = info(id).land_pos; # aggiungo landmark allo stato
            
       endif

   endfor
   id_to_state = id_to_state_map_new;
   state_to_id = state_to_id_map_new;
   #printf('Landmarks in XL: %i \n', size(XL,2));
   printf('Valid landmarks: %i \n', land_count);  
endfunction;

[id_to_state_map, state_to_id, state_landmark] = update_map(info, id_to_state_map, state_to_id_map);
printf('Dimension of state landmarks: %i \n', length(state_landmark));

id_to_state_map(:)

for j=1:length(id_to_state_map)
   if id_to_state_map(j) == 0
      j, id_to_state_map(j)
   endif
   
   endfor



# ------ LEAST SQUARES ------- #







function [XR, XL]=boxPlus(XR, XL, num_poses, num_landmarks, dx)
   global pose_dim;
   global landmark_dim;
   
   for(pose_index=1:num_poses)
      pose_matrix_index=poseMatrixIndex(pose_index,
      num_poses,
      num_landmarks);
      dxr=dx(pose_matrix_index:pose_matrix_index+pose_dim-1);
      XR(:,:,pose_index)=v2t(dxr)*XR(:,:,pose_index);
   endfor;
   
   for(landmark_index=1:num_landmarks)
      landmark_matrix_index=landmarkMatrixIndex(landmark_index,
      num_poses,
      num_landmarks);
      dxl=dx(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,:);
      XL(:,landmark_index)+=dxl;
   endfor;
endfunction;







function [e,J]=errorAndJacobian(x,p,z)
   t=x(1:2);
   theta=x(3);
   R=rotation2D(theta);
   z_hat=R*p+t;
   e=z_hatz;
   J=zeros(2,3);
   J(1:2,1:2)=eye(2);
   J(1:2,3)=rotation2Dgradient(theta)*p;
endfunction;





%system_size = 3*size(XR, 3) + 2*size(XL,2)
%H=zeros(system_size, system_size);
%b=zeros(system_size,1);
%chi_stats(iteration)=0;
%for (measurement_num=1:size(Z,2))
%   
%   pose_index=associations(1,measurement_num);
%   landmark_index=associations(2,measurement_num);
%   
%   z=Z(:,measurement_num);
%   Xr=XR(:,:,pose_index);
%   Xl=XL(:,landmark_index);
%   [e,Jr,Jl] = errorAndJacobian(Xr, Xl, z);
%   Hrr=Jr'*Jr;
%   Hrl=Jr'*Jl;
%   Hll=Jl'*Jl;
%   br=Jr'*e;
%   bl=Jl'*e;
%   
%   pose_matrix_index=poseMatrixIndex(pose_index, num_poses, num_landmarks);
%   
%   landmark_matrix_index=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks);
%   
%   H(pose_matrix_index:pose_matrix_index+pose_dim-1,
%   pose_matrix_index:pose_matrix_index+pose_dim-1)+=Hrr;
%   
%   H(pose_matrix_index:pose_matrix_index+pose_dim-1,
%   landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=Hrl;
%   
%   H(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,
%   landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=Hll;
%   
%   H(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,
%   pose_matrix_index:pose_matrix_index+pose_dim-1)+=Hrl';
%   
%   b(pose_matrix_index:pose_matrix_index+pose_dim-1)+=br;
%   
%   b(landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=bl;
%endfor







function rot = rotation(angle)
   
   rot = [ cos(angle), -sin(angle);
           sin(angle), cos(angle)];

endfunction







 # DECIDI COSA FARE 
state = state_poses;
num_poses = size(state_poses,1)/3
#state = curr_state;
printf('Dimension of state poses: %i\n', size(state_poses,1));
num_landmarks = size(state_landmark,1)/2
state(end+1:end+size(state_landmark,1),1) = state_landmark;
printf("State total size %i \n", size(state,1));
   


function transform = v2t(pose)
   transform(1:2,1:2) = rot(pose(3));
   transform(1:2,3) = [pose(1); pose(2)];
   transform(3,:) = [0,0,1];
endfunction   
   
   
   
   
# trasformazione dalla posa del robot wrt world ad XR e posa del world wrt robot
# + XL  
function [XR, XL] = w2rob(state_poses, state_landmark)
   for r = 1:length(state_poses)/3
      
      pose_index = 3*(r-1) +1;
      theta = state_poses(pose_index+2);
      R = [ cos(theta), -sin(theta);
              sin(theta), cos(theta)];
      R_t = R';
      pose = state_poses(pose_index:pose_index +2);
      XR(1:3,1:3,r) = [R_t, -R_t*pose(1:2);
       0, 0, 1];
      pose_w = XR(1:2,3,r);
      pose_w(3) = atan(XR(2,1,r)/XR(1,1,r));
   endfor
   XR;
   for l =1:length(state_landmark)/2
      
      land_index = 2*(l-1) +1;
      XL(:, l) = state_landmark(land_index:land_index+1);
   endfor   
      
   
endfunction
   
   
   
%   
%for g =1:length(state_poses)/3
%  pose_index = 3*(g-1) +1;
%  XR(:,:,g) = v2t(state_poses(pose_index:pose_index+2));
%
%for l =1:length(state_landmark)/2
%   land_index = 2*(g-1) +1;
%   XL(:, l) = state_landmark(land_index:land_index+1);
%
%




























   
   
system_size = size(state_poses, 1) + size(state_landmark,1);
H=zeros(system_size, system_size);
b=zeros(system_size,1);
[XR, XL] = w2rob(state_poses, state_landmark);

for c =1:length(observations)
  c;
  #pose_index = 3*(t-1) +1;
  #rob_pose = state(pose_index: pose_index + 2);
  pose_index = c;
  R = XR(1:2,1:2,pose_index);
  t = XR(1:2,3,pose_index);
%  H = 0;
%  b = 0;
%  H=zeros(system_size, system_size);
%  b=zeros(system_size,1);
%  R = [ cos(rob_pose(3)), -sin(rob_pose(3));
%           sin(rob_pose(3)), cos(rob_pose(3))];
%  R_t = R';
%  R_t_theta = [-sin(rob_pose(3)), cos(rob_pose(3));
%               -cos(rob_pose(3)), -sin(rob_pose(3))];
  
  for ob = observations(c).observation 
  # struct array con id e bearing del landmark
%  obs = observations(t).observation;
  # posa del robot
  
  
     
     if (ob.id) == 0
        ob.id = 200; # TODO: fix
     endif
     controlla se l'id del landmark è tra quelli accettati
%     pose_index = 3*(t-1) +1;
     
     landmark_index = id_to_state_map(ob.id) ; #cerco landmark nello stato
     if landmark_index == 0
        ob.id
        endif
     landmark = XL(:, landmark_index);
%     delta = landmark - state(pose_index:pose_index+1,1);
     
     #---#
%     p_i_w = landmark;
     p_i = R*landmark + t;
     x_i = p_i(1);
     y_i = p_i(2);
     #p_i = R_t*landmark + state_poses(pose_index:pose_index+1,1);
     
     
     h = atan2(p_i(2,1), p_i(1,1)); # prediction
     z = ob.bearing; # measurement
     
     e = h - z; # error
     rot_e = rotation(h)*rotation(z)';
     e = atan2(rot_e(2,1),rot_e(1,1));
%     e = acos((trace(rot_e) - 1)/2)
     
%     if e <0
%        e = -e
%     endif
     #---#
%     x_l = landmark(1,1);
%     y_l = landmark(2,1);
     
     
%     R_t_prime = ones(2,2);  
%     R_t_prime(1:2, 3) = R_t_theta*p_i;
%     f(1:2,1:2) = R_t;
%     f(1:2, 3) = -R_t*p_i;
%     f(3,1:3) = (0,0,1);
%     R_t_prime(1:2, 3) = f
     omega = 1;
     
     
%     [-R_t , R_t_theta*delta];

     Jr = 1/(x_i**2 + y_i**2) * [-y_i, x_i] * [ones(2) , [-y_i; x_i]];
     size(Jr);
     
     Jl = 1/(x_i**2 + y_i**2) * [-y_i, x_i] * R;
     size(Jl);
     
     
     Hrr=Jr'*Jr;
%     chol(Hrr);
     Hrl=Jr'*Jl;
     Hll=Jl'*Jl;
     br=Jr'*e;
     bl=Jl'*e;
     
     pose_dim = 3;
     landmark_dim = 2;
     
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
     
     
     
%     if c==3
%     break
%      endif
     
     
%     
%     H += J' * omega * J;
%  
%     b += J' * omega * e;
%     
     # costruire la matrice H come composizione di matrici, una per 
   
  endfor
  
%  size(H)
%  size(b)
%  chol(H)
%  dx = -H\b;
%%
  state(1:3);  
%  state += dx;
%  state(1:3);
%break 
endfor



%issymmetric(H)
%dx = -H\b;
%state += dx;
%for i=1:size(H,1)
%   if H(i,i) == 0
%      i, H(i,i);  
%   endif
%endfor
state; 
%dx = -pinv(H)*b;
dx = -(H\b);

























# ----------------------------------------------- #

%# OLD FUNCTION
%
%for t=1:length(observations)
%  t;
%  pose_index = 3*(t-1) +1;
%%  state(1:3,1) = [poses(t+1).x; poses(t+1).y; poses(t+1).theta];
%  rob_pose = state(pose_index: pose_index + 2);
%%  H = 0;
%%  b = 0;
%%  H=zeros(system_size, system_size);
%%  b=zeros(system_size,1);
%  R = [ cos(rob_pose(3)), -sin(rob_pose(3));
%           sin(rob_pose(3)), cos(rob_pose(3))];
%  R_t = R';
%  R_t_theta = [-sin(rob_pose(3)), cos(rob_pose(3));
%               -cos(rob_pose(3)), -sin(rob_pose(3))];
%  
%  for ob = observations(t).observation 
%  # struct array con id e bearing del landmark
%%  obs = observations(t).observation;
%  # posa del robot
%  
%  
%     
%     if (ob.id) == 0
%        ob.id = 200; # TODO: fix
%     endif
%     
%%     pose_index = 3*(t-1) +1;
%     
%     landmark_index = size(state_poses, 1) +1 + 2*(id_to_state_map(ob.id) -1); #cerco landmark nello stato
%     
%     landmark = state(landmark_index: landmark_index+1, 1);
%     delta = landmark - state(pose_index:pose_index+1,1);
%     
%     #---#
%     p_i_w = landmark;
%     p_i = R_t*delta;
%     x_i = p_i(1);
%     y_i = p_i(2);
%     #p_i = R_t*landmark + state_poses(pose_index:pose_index+1,1);
%     
%     
%     h = atan2(p_i(2,1), p_i(1,1)) # prediction
%     z = ob.bearing # measurement
%     
%     e = h - z; # error
%%     rot_e = rotation(z)'*rotation(h);
%%     
%%     e = acos((trace(rot_e) - 1)/2)
%     
%%     if e <0
%%        e = -e
%%     endif
%     #---#
%%     x_l = landmark(1,1);
%%     y_l = landmark(2,1);
%     
%     
%     R_t_prime = ones(2,2);  
%     R_t_prime(1:2, 3) = R_t_theta*p_i;
%%     f(1:2,1:2) = R_t;
%%     f(1:2, 3) = -R_t*p_i;
%%     f(3,1:3) = (0,0,1);
%%     R_t_prime(1:2, 3) = f
%     omega = 1;
%     
%     
%     [-R_t , R_t_theta*delta];
%
%     Jr = 1/(x_i**2 + y_i**2) * [-y_i, x_i] * [-R_t , R_t_theta*delta];
%     size(Jr);
%     
%     Jl = 1/(x_i**2 + y_i**2) * [-y_i, x_i] * R_t;
%     size(Jl);
%     
%     
%     Hrr=Jr'*Jr;
%%     chol(Hrr);
%     Hrl=Jr'*Jl;
%     Hll=Jl'*Jl;
%     br=Jr'*e;
%     bl=Jl'*e;
%     
%     pose_dim = 3;
%     landmark_dim = 2;
%     
%     H(pose_index:pose_index+pose_dim-1,
%     pose_index:pose_index+pose_dim-1)+=Hrr;
%     
%     H(pose_index:pose_index+pose_dim-1,
%     landmark_index:landmark_index+landmark_dim-1)+=Hrl;
%     
%     H(landmark_index:landmark_index+landmark_dim-1,
%     landmark_index:landmark_index+landmark_dim-1)+=Hll;
%      
%     H(landmark_index:landmark_index+landmark_dim-1,
%     pose_index:pose_index+pose_dim-1)+=Hrl';
%      
%     b(pose_index:pose_index+pose_dim-1)+=br;
%      
%     b(landmark_index:landmark_index+landmark_dim-1)+=bl;
%     
%     
%     
%     if t==3
%      break
%      endif
%     
%     
%%     
%%     H += J' * omega * J;
%%  
%%     b += J' * omega * e;
%%     
%     # costruire la matrice H come composizione di matrici, una per 
%   
%  endfor
%  
%%  size(H)
%%  size(b)
%%  chol(H)
%%  dx = -H\b;
%%%
%  state(1:3);  
%%  state += dx;
%%  state(1:3);
%break 
%endfor
%
%
%
%%issymmetric(H)
%%dx = -H\b;
%%state += dx;
%%for i=1:size(H,1)
%%   if H(i,i) == 0
%%      i, H(i,i);  
%%   endif
%%endfor
%state; 
%%dx = -pinv(H)*b;
%dx = -(H\b);

# ----------------------------------------------- #
 
   
%chol(H)
%state += dx;
%state(1:3); 
   
H(1170:end, 1170:end)  ; 
   


  
  
  
  
  
  
  
  
  
  
  