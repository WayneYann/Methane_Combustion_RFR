Files ending in No_Diffusion are for the simplified model that excludes diffusion.
ode15s_JPJ is a modefied version of ode15s to calculate the eigenvalues and output the condition numbers for the Jacobian.
ode15s_JPJ will need to be placed in the following MATLAB directory.

C:\Program Files\MATLAB\R2015b\toolbox\matlab\funfun.

1)You will need to temporarly change the properties of the folder to give you write permssion.
a) Right click in the funfun folder space and click on properties.
b) Go to security
c) select your user (such as Users (Josh-C\Users))
d) click edit
e) click allow for write permissions
f) unlick after saving the file to the space so you don't accidentally mess up your MATLAB home directory.


If you wish, you can replace ode15s_JPJ in the code with just ode15s.
If you do not replace it you will not get a graph of the condition numbers.


Files work like this
1) Run JPJ_RFR_No_Diffusion
2) This calls JPJ_ODEs_No_Diffusion and sends the function to ode15s_JPJ
3) JPJ_variables_no_Diffusion is run in both JPJ_RFR_No_Diffusion and JPJ_ODEs_No_Diffusion to save space in the main code.
