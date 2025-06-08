Master's Thesis: Modeling and Control of Steam Bottoming Cycles with Moving Point of Vaporization

In this thesis, a dynamic model of a steam bottoming cycle for offshore implementations was developed in CasADi/MATLAB and solved using the IDAS solver from SUNDIALS. 
The model was based on simplifying assumptions such as ideal gas and constant heat capacities. 
the model has the ability to change equations used in the solver based on the vapor quality, allowing the point of vaporization to move along the once-through steam generator (OTSG).
Several control structures and operational strategies for power were tested and compared to relevant studies. 

The reader should consult the thesis before diving into the code. 
The given directory contains code, simulation result files, and plots used in the thesis. 
Each of the cases has its own folder. 

The model comprises of three files: an initModel.m, Model.m and U_calc.m, where Model denotes the name of the case.
initModel.m  loads constants, initial state of the system variables and calls the Model.m file to solve system over the given step time.
Model.m is a function that define, build and solves the DAE over the given step in time. 
U_calc.m contains a function for calculation of the overall heat transfer coefficent. This function is called by Model.m.
The SteamCycleCaseconstant.m and initSteamCycleCaseconstant.m contains a more in-depth explaination of each part of the code and one should start with this case to understand the code structure. 
