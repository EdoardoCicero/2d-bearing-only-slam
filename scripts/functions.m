# ---------------- Update mapping ------------- #

# function that updates id_to_state and state_to id mapping 
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
   printf('Valid landmarks: %i \n\n', land_count);  
   
endfunction;



# -------------------------------------------------------------------------- #

# rotation function, from angle to rotation matrix
function rot = rotation(angle)
   
   rot = [ cos(angle), -sin(angle);
           sin(angle), cos(angle)];
endfunction

# -------------------------------------------------------------------------- #

# -- Boxplus function -- #
# pose_index, landmark_index = indices of pose and landmark in XR and XL
# p, l = indices of pose and landmarks in dx
# output:
# - new [XR, XL]    
function [XR, XL]=boxPlus(XR, XL, num_poses, num_landmarks, dx)
   pose_dim = 3;
   landmark_dim = 2;
   
   for(pose_index=1:num_poses)
      p = (pose_index -1)*3 +1;
      dxr=dx(p:p+pose_dim-1);
      XR(:,:,pose_index)=v2t(dxr)*XR(:,:,pose_index);
   endfor;
   
   for(landmark_index=1:num_landmarks)
      l = num_poses*3+ 2*(landmark_index -1) +1 ;
      dxl=dx(l:l+landmark_dim-1);
      XL(:,landmark_index)+=dxl;
   endfor;
   
endfunction;   



# -------------------------------------------------------------------------- #



# From pose of origin wrt robot to pose of robot wrt origin + XL  
function [XR, XL] = w2rob(state_poses, state_landmark)
   for r = 1:length(state_poses)/3
      
      pose_index = 3*(r-1) +1; # pose index in state_poses
      theta = state_poses(pose_index+2);
      
      R = [ cos(theta), -sin(theta);
              sin(theta), cos(theta)];
      R_t = R';
      t = state_poses(pose_index:pose_index+1);

      XR(1:3,1:3,r) = [R_t, -R_t*t ; 0, 0, 1]; # calculate inverse of XR (it's the same as inv(XR))

   endfor
   XR;
   for l =1:length(state_landmark)/2
      
      land_index = 2*(l-1) +1;
      XL(1:2, l) = state_landmark(land_index:land_index+1);
   endfor   
      
   
endfunction
   

# ------------ F: Error function -------------------- #

function error_sum = F(XR, XL, observations, state_to_id_map, id_to_state_map)
   error_sum = [];
   for c=1:length(observations)
      
       p_index = c;
       R = XR(1:2,1:2,p_index);
       t = XR(1:2,3,p_index);
       
      for ob = observations(c).observation 
           if (ob.id) == 0
              ob.id = length(id_to_state_map);
           endif
    
           l_index = id_to_state_map(ob.id) ; # index of landmark in XL
           
           if l_index == 0 # id not mapped
              ob.id;
              continue
           endif

           landmark = XL(:, l_index); # get landmark from XL

           p_i = R*landmark + t;
           x_i = p_i(1);
           y_i = p_i(2);
          
           h = atan2(y_i, x_i); # prediction
           z = ob.bearing; # measurement
           
           e = h - z; # error
           error_sum(end+1,1) = e; #adding e to vector of errors
      endfor
endfor

error_sum = sqrt(sumsq(error_sum)); # squared error
#printf("Error: %f\n", error_sum);

endfunction        

# -------------------------------------------------------------------------- #


# gettin the robot poses in origin reference frame   
function XR = rob2w(XR)
   for h = 1:size(XR,3)

      XR(1:2,1:2,h) = XR(1:2,1:2,h)';
      XR(1:2,3,h) = -XR(1:2,1:2,h)*XR(1:2,3,h);
      
      state_rob(1:3,h) = t2v(XR(:,:,h));
      XR(:,:,h) = v2t(state_rob(1:3,h));
   endfor
   state_rob;
   
endfunction


# ---------------------------------------------------------#

# PLOTTING

# static plot of all the map and the trajectory
function end_plot = static_plot(XR, XL,state_to_id_map)
   XR_ = rob2w(XR);
   
for u=1:size(XR_,3)
   pose = XR_(1:3,3,u);
   
   if u ==1
    trajectory = [pose(1), pose(2)];
   else
    trajectory = [trajectory; pose(1), pose(2)];
   endif
endfor  
#trajectory
for u=1:size(XL,2)
   landmarks(end+1).id = state_to_id_map(u);
   landmarks(end).x_pose = XL(1,u);
   landmarks(end).y_pose = XL(2,u);
endfor
plotState(landmarks, XR_(1:3,3,end), eye(2)*0.1, eye(2), trajectory);

pause(60);
fflush(stdout);

endfunction



 # plotting the trajectory of the robot, pose by pose, spawning each landmark from the second time it has been seen

function init_plot = dynamic_plot(XR, XL, state_to_id_map, id_to_state_map, observations)
XR = rob2w(XR);
bearing_count = zeros(1, 300);
plot_time=0;
for o=1:length(observations)
   
    p_index = o;
    R = XR(1:2,1:2,p_index);
    t = XR(1:2,3,p_index);
    pose = XR(1:3,3,o);
   
   if o ==1
    trajectory = [pose(1), pose(2)];
   else
    trajectory = [trajectory; pose(1), pose(2)];
   endif
   
   for ob = observations(o).observation 
        if (ob.id) == 0
           ob.id = length(id_to_state_map);
        endif

        l_index = id_to_state_map(ob.id) ; # index of landmark in XL
        
        if l_index == 0 
           continue
        endif
        
        bearing_count(1,ob.id) +=1;
           
        if bearing_count(1,ob.id) ==2
           plot_time = 1;

           landmark = XL(:, l_index);
           landmarks(end+1).id = state_to_id_map(l_index);
           landmarks(end).x_pose = XL(1,l_index);
           landmarks(end).y_pose = XL(2,l_index);
        endif
    endfor
    
if plot_time      
  plotState(landmarks, pose, eye(2)*0.1, eye(2), trajectory);
endif  
#pause(0.1);
fflush(stdout);
        
endfor  

endfunction
# ---------------------------------------------------------