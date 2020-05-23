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
pose_structure = poses;
# robot observations
observation_structure = observations;

%rob_pose = [poses(1).x; poses(1).y; poses(1).theta];
%#disp(rob_pose);
%printf("initial pose guess: [%f, %f, %f]\n", rob_pose(1), rob_pose(2), rob_pose(3));

#bookkeeping: to and from mapping between robot pose (x,y, theta) and landmark indices (i)
#all mappings are initialized with invalid value -1 (meaning that the index is not mapped)
#since we do not know how many landmarks we will observe, we allocate a large enough buffer
id_to_state_map = ones(300, 1)*-1;
state_to_id_map = ones(300, 1)*-1;
rob_pose = [poses(1).x; poses(1).y; poses(1).theta];
curr_state = rob_pose;

# quest due variabili contengono
# l'ultima e la penultima misurazione dei landmark
prev_land_meas = ones(300, 1)*-1;
last_land_meas = ones(300, 1)*-1;

# quest due variabili contengono
# l'ultima e la penultima posa del robot per ogni landmark osservato
prev_rob_pose = nan(300, 3);
last_rob_pose = nan(300, 3);

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
for t = 1:length(observations)
  # struct array con id e bearing del landmark
  obs = observations(t).observation;
  # posa del robot ad ogni osservazione
  rob_pose = [poses(t+1).x, poses(t+1).y, poses(t+1).theta];

  # tutte le pose del robot   
  state_poses(end+1:end+3, 1) = rob_pose;
  
  # la variabile info contiene struct array con: id del landmark, ogni osservazione del landmark, ogni posa del robot
  # da cui è stata fatta ogni osservazione l'osservazione
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

# aggiungo landmark con id=0, ora è id=length(info)
info(end +1).id = zero.id;
info(end).bearing = zero.bearing;
info(end).rob_pose = zero.rob_pose;
n++;
id_to_state_map(length(info)) = n;
state_to_id_map(n) = length(info);

printf('Num landmarks: %i \n', n)


# ---- Get landmarks position --- #


# theta_w = angolo di bearing rispetto al RF del mondo ( theta_w = rob_pose.theta + land.bearing )
# (x_r,y_r) = posizione del robot nel sistema del mondo quando osserva il bearing angle del landmark
# voglio y=mx+q in world -> y_w=tg(theta_w)(x_w - x_r) + y_r


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

      # aggiungo posizione del landmark ad info per quelli che hanno almeno 2 misurazioni di bearing angle
      if rows(A) >= 2 
         pinv(A)*b;
         info(l).land_pos = pinv(A)*b;
      else
         info(l).id;
         discarded_landmarks(end+1) = info(l).id;

      endif

   endif

endfor
# print discarded landmarks due to few observations
discarded_landmarks



# ---------------- Update mapping ------------- #

# function that updates id_to_state and state_to id mapping considering only initialized landmarks (with at least 2 measurements)
# input:
# - info structure
# - id to state mapping
# - state to id mapping

# output:
# - new id to state mapping
# - new state to id mapping
# - state_landmark: vector containing all landmarks that are in the state
function [id_to_state, state_to_id, state_landmark] =  update_map( info, id_to_state_map, state_to_id_map)

   land_count = 0;
   for l =1:length(state_to_id_map)
      id = state_to_id_map(l);

      # accepting a landmark only if the corrisponding info has landmark position field
      # so the landmark is initialized
      if id != -1 && size(info(id).land_pos,1) !=0        
            land_count ++;
            state_to_id_map_new(land_count) = id;
            id_to_state_map_new(id) = land_count;
           
            state_landmark(end+1: end +2, 1) = info(id).land_pos; # adding landmark to landmarks in the state
            
       endif

   endfor
   id_to_state = id_to_state_map_new;
   state_to_id = state_to_id_map_new;
   printf('Valid landmarks: %i \n', land_count);  
   
endfunction;

[id_to_state_map, state_to_id_map, state_landmark] = update_map(info, id_to_state_map, state_to_id_map);

printf('Dimension of state landmarks: %i \n', length(state_landmark));



# ----------------- LEAST SQUARES ----------------- #


%function [e,J]=errorAndJacobian(x,p,z)
%   t=x(1:2);
%   theta=x(3);
%   R=rotation2D(theta);
%   z_hat=R*p+t;
%   e=z_hatz;
%   J=zeros(2,3);
%   J(1:2,1:2)=eye(2);
%   J(1:2,3)=rotation2Dgradient(theta)*p;
%endfunction;





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



# rotation function, from angle to rotation
function rot = rotation(angle)
   
   rot = [ cos(angle), -sin(angle);
           sin(angle), cos(angle)];
endfunction



state = state_poses;
num_poses = size(state_poses,1)/3;
printf('Robot poses: %i\n', num_poses);
printf('Dimension of robot poses: %i\n', size(state_poses,1));

num_landmarks = size(state_landmark,1)/2;

state(end+1:end+size(state_landmark,1),1) = state_landmark;
printf("State total size %i \n", size(state,1));
   

# function from vector to transformation
function T = v2t(pose)
   T(1:2,1:2) = rotation(pose(3));
   T(1:2,3) = [pose(1); pose(2)];
   T(3,:) = [0,0,1];
endfunction   


# function from transformation to vector
function v = t2v(T)
   v(1:2,1) = T(1:2, 3);
   v(3,1) = atan(T(2,1)/T(1,1));
   
endfunction
   
   
   
# boxplus function    
function [XR, XL]=boxPlus(XR, XL, num_poses, num_landmarks, dx)
   pose_dim = 3;
   landmark_dim = 2;
   
   for(pose_index=1:num_poses)
      dxr=dx(pose_index:pose_index+pose_dim-1);
      XR(:,:,pose_index)=v2t(dxr)*XR(:,:,pose_index);
   endfor;
   
   for(landmark_index=1:num_landmarks)
      l = num_poses*3+ 2*(landmark_index -1) +1 ;
      dxl=dx(l:l+landmark_dim-1,:);
      XL(:,landmark_index)+=dxl;
   endfor;
   
endfunction;   
   
   
   
   

   
   
   
   
   
   
   
   
   
   
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




























[XR, XL] = w2rob(state_poses, state_landmark);
system_size = num_poses*3 + num_landmarks*2;
 
num_iteration = 3;
for it=1:num_iteration
   
   
   H=zeros(system_size, system_size);
   b=zeros(system_size,1);


   for c =1:length(observations)
     c;
     #pose_index = 3*(t-1) +1;
     #rob_pose = state(pose_index: pose_index + 2);
     p_index = c;
     R = XR(1:2,1:2,p_index);
     t = XR(1:2,3,p_index);
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
           ob.id = 200;
                         # TODO: fix
        endif
        
        
   %     pose_index = 3*(t-1) +1;
        
        l_index = id_to_state_map(ob.id) ; # index of landmark in XL
        
        if l_index == 0 # cioè l'id non è mappato a nessun landmark nello stato
           ob.id;
           continue
        endif
        ob.id;
        l_index;
        state_to_id_map(l_index);
        landmark = XL(:, l_index);
   %     delta = landmark - state(pose_index:pose_index+1,1);
        
        #---#
   %     p_i_w = landmark;
        p_i = R*landmark + t;
        x_i = p_i(1);
        y_i = p_i(2);
        #p_i = R_t*landmark + state_poses(pose_index:pose_index+1,1);
        
        
        h = atan2(p_i(2), p_i(1)); # prediction
        z = ob.bearing; # measurement
        
        e = h - z; # error
%        rot_e = rotation(h)*rotation(z)';
%        e = atan2(rot_e(2,1),rot_e(1,1));

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
   %     omega = 1;
        
        
   %     [-R_t , R_t_theta*delta];

        Jr = 1/(x_i**2 + y_i**2) * [-y_i, x_i] * [eye(2) , [-y_i; x_i]];
        size(Jr);
        
        Jl = 1/(x_i**2 + y_i**2) * [-y_i, x_i] * R;
        size(Jl);
        
        
        Hrr=Jr'*Jr;
        Hrl=Jr'*Jl;
        Hll=Jl'*Jl;
        br=Jr'*e;
        bl=Jl'*e;
        
        pose_dim = 3;
        landmark_dim = 2;
%        p_index
        pose_index = 3*(p_index -1) + 1;
%        l_index;
        landmark_index = num_poses*3 + 2*(l_index -1) + 1;
        
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


     endfor

%   break
   endfor

%#   dx = -(H\b);
   #dx=zeros(system_size,1);
   size(H(pose_dim+1:end,pose_dim+1:end));
   size(b(pose_dim+1:end,1));
   large_value = 10;
   lambda = 0.0001;
%   H(1:pose_dim,1:pose_dim)+=eye(3)*large_value;
   dx = -((H + ones(system_size)*lambda) \b);
   dx(901:end);
   #dx(pose_dim+1:end) = -(H(pose_dim+1:end,pose_dim+1:end)\b(pose_dim+1:end,1));

%   state += dx;

   [XR, XL]=boxPlus(XR, XL, num_poses, num_landmarks, dx);
   XR(:,:,1:3); 
endfor
   
   
# ottengo le pose del robot wrt world in state_rob
for h = 1:size(XR,3)

   XR(1:2,1:2,h) = XR(1:2,1:2,h)';
   XR(1:2,3,h) = -XR(1:2,1:2,h)*XR(1:2,3,h);
   
   state_rob(1:3,h) = t2v(XR(:,:,h));
   
endfor

state_rob
%for s=1:length(state_to_id_map)
%   state_to_id_map(s),XL(:,s)
%endfor














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
   


  
  
  
  
  
  
  
  
  
  
  