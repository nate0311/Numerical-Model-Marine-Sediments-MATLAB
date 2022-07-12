# Numerical-MATLAB-Model
A numerical approach to a partial differential equation that describes the transport of a radioisotope in marine sediments
This is a MATLAB script for Cochran and Krishnaswami (1980) reaction-transport equation:

∂C/∂t=(Ds+KDb)(∂^2 C)/(∂z)^2 - λ(1+K)C + P        for 0 ≤ z ≤ L					

where:
C = concentration of dissolved radioisotope in porewater (atoms cm-3)
Ds = molecular diffusion corrected for tortuosity (cm2 yr-1)
Db = bioturbation (cm2 yr-1); assumed constant in the upper layer 
λ = decay constant of daughter (yr-1)
P = production rate of daughter atoms per unit of porewater from parent decay (atoms cm-3 yr-1)
K = partition coefficient (dimensionless) 
z = depth in sediments (cm)
L = thickness of bioturbation zone (cm)	

The MATLAB code uses a numerical approach to solve the above reaction-transport equation and only uses the upper few cm of sediments to model 227Ac. 
Furthermore, the numerical model assumes S=0 since S is very small (<0.25 cm/kyr). Known variables are first defined. These include 227Ac and 231Pa 
decay constants, bioturbation rate (Db), molecular diffusion (Dm), F, 231Pa activity in mix layer, and distribution coefficient (kd). Next, an array 
of porosity for every depth modeled was created based on fitting a logarithmic function to the measured porosity vs. depth profile. Next K and Ds 
values were assigned to every depth since these values depend on porosity. Lastly, initial conditions were applied to the model and the model ran 
for 100 years, producing 227Ac values for every 0.1 cm in the top 10 cm.  

The model is created to change variables and observe how the 227Ac profiles behaves relative to the change. For this model, it is optimized to 
change F, Db, and kd. The model takes in the measured 227Ac values in the upper few cm of sediments (input 227Ac measured values into variable B1) 
and compares it to the model results for each cm. 


![image](https://user-images.githubusercontent.com/109116048/178394346-958ad301-7cfc-41d1-936f-64053e5775d0.png)
