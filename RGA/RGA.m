% RGA Calculation Script

% Define the process gain matrix G 
%Openloop
%u1 = zm, u2 = zv
%y1 = Ts, y2 = ps
% g11 = zm-Ts, g12 = zv-Ts
% g21 = zm-ps, g22 = zv-ps

%FC closed
%u1 = mpSP, u2 = zv
%y1 = Ts,   y2 = ps
% g11 = mpSP-Ts, g12 = zv-Ts
% g21 = mpSP-ps, g22 = zv-ps

G = [-52.4, 2.9;  
     -1.25,  -1.05];

% Compute the Moore-Penrose pseudoinverse of G
G_inv = pinv(G);

% Calculate the RGA
RGA = G .* G_inv.';

% Display the results
disp('Process Gain Matrix G:');
disp(G);

disp('Pseudoinverse of G:');
disp(G_inv);

disp('Relative Gain Array (RGA):');
disp(RGA);
