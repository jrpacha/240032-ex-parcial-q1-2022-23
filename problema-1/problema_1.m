%% Problema 1 test options
% 

clearvars
close all
format short e
%
% Data
%
L1 = 3.0; L2 = 1.0; 
w = 0.01; d = 1.0; ks = 1.0; 

% #1
% h1 = 2.0; h2 = 4.0; E = 2.0;
% #2
% h1 = 2.0; h2 = 4.0; E = 2.0; 
% #3
% h1 = 1.0; h2 = 2.0; E = 1.0; 
% #4
% h1 = 2.0; h2 = 4.0; E = 2.0; 
% #5
% h1 = 2.0; h2 = 3.0; E = 2.0; 
% #6
% h1 = 2.0; h2 = 3.0; E = 2.0; 
% #7
% h1 = 2.0; h2 = 3.0; E = 2.0;
% #8
% h1 = 1.0; h2 = 2.0; E = 1.0;
% #9
% h1 = 2.0; h2 = 4.0; E = 2.0;
% #10
h1 = 1.0; h2 = 4.0; E = 1.0; 

%
% Geometry
%
A1 = L1*L2;
A2 = L2*L2;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part A
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Element 1
alpha = -(A1-A2)/h1;
beta = A1;                     %A(x) = alpha * x + beta
x1_elem_1 = 0; x2_elem_1 = h1; %Positions element one's nodes
K1 = 0.5*alpha*E*(x1_elem_1+x2_elem_1)*[1,-1;-1,1]/h1 + ...
    beta*E*[1,-1;-1,1]/h1;

%Element 2
K2 = E*A2*[1,-1;-1,1]/h2;

%Assemble matrices
K = zeros(3);
K(1:2,1:2) = K1;
K(2:3,2:3) = K(2:3,2:3)+K2
fprintf(['Part A\n',...
    'The entry K(2,2) of the stiff matrix of the global system is ',...
    'K(2,2) = %.4e\n',...
    'Hint 1. The value of K(2,3) is K(2,3) = %.4e\n'],K(2,2),K(2,3))

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part B 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Element 1
F1 = -alpha*w*h1*[2*x1_elem_1+x2_elem_1; x1_elem_1+2*x2_elem_1]/6 + ...
    -beta*w*h1*[1;1]/2;
% Element 2
F2 = -w*A2*h2*[1;1]/2.0;

% Assemble forces
F = zeros(3,1);
F(1:2) = F1;
F(2:3) = F(2:3) + F2
fprintf(['Part B\n',...
    'The value of F(2) is, F(2) = %.4e\n',...
    'Hint 2. The value of F(3) is, F(3) = %.4e\n'],F(2),F(3))

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part C
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Km = K+[ks,0,0;0,0,0;0,0,ks];
Qm = F+d*ks*[-1;0;1];
U = Km\Qm
fprintf(['Part C\n',...
    'Displacement U(2) of the global node 2, U(2) = %.4e\n',...
    'Hint 3. The value of U(1) is U(1) = %.4e\n'],U(2),U(1))