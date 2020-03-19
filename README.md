# 2D-gait
Generate torque driven 2-D gait simulations in below 5s.

Model is based on the one presented in "An Introduction to Trajectory Optimization: How to Do Your Own Direct Collocation" by Matthew Kelly. (SIAM Review, Vol. 59, No. 4, pp. 849-904)
Model exists of 5 segments (2 tibia, 2 femur and trunk) and driven by ankle, knee and hip torques.
The stance 'foot' is fixed to the ground and the gait pattern is imposed to have no double stance and no flight phase avoiding the implementation of different phases and/or contact models.

We use CasADi to generate and solve an NLP from our collocation simulation. We don't need to provide analytic gradients of our constraints or cost function as CasADi does this for us using automatic differentiation.

If you run the 'main' function you will see how the minimial torque walking pattern looks like. There are some options to restrict the gait pattern and look at the effects of for example ('walking without knee muscles'). You can uncomment these and run to observe the effects.

The functions to generate the equations of motion and heel-strike map constraint (to have a consistent periodic gait) are included ('generateEOM.m'). 
