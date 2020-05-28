# 2d-bearing-only-slam

From console: octave 2D_Slam.m

- From the dataset all the poses and the observations were extracted through "loadG2o", then it was built a structure array named info containing for each landmark its id, all the poses from which it was observed, all the bearing angles measured by the robot in those poses.
- Along with that the mapping from id to state and state to id were created.
- From info the initial guess of the landmarks position was calculated, considering all observations for the single landmarks.
- Landmarks with less than 2 observations were deleted and new mappings id-to-state and state-to-id were done to account the loss of landmarks.
- The poses of robot wrt origin were transformed into transformation matrices (XR) and inverted, to have the poses of the origin wrt robot, to simplify calculations.
- With XL and XR the Least-Squares algorithm is computed and the trajectory is plotted.
