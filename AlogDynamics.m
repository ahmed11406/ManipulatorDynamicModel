%Alogrithm for computoing the dynamics model. 
%This code is created by Ahmed Ali
clearvars; 
%---------------------------------START------------------------------------
%--------------------------------INPUTS------------------------------------
    n=3;                       %n= #links.
    k=[0;0;0]; %rev=0, pris=1;
    syms dq1 dq2 dq3 real; % real dq3 dq4 dq5 real;
    dq = [dq1; dq2; dq3];% dq3; dq4];
    syms q1  q2  q3 real ;%q3  q4 real;
    q = [q1; q2; q3];% q3; q4];
    syms m1  m2 m3  real ;%m3  m4 real;
    m=[m1; m2; m3];% m3; m4];
    syms L1  L2  L3 real ;% L3  L4 real;
    syms Ic1x Ic1y Ic1z  Ic2x Ic2y Ic2z Ic3x Ic3y Ic3z real; % Ic3  Ic4 real;
    Ic=cell(3,1);
    Ic{1,1}=diag([Ic1x,Ic1y,Ic1z]);
    Ic{2,1} =diag([Ic2x,Ic2y,Ic2z]);
    Ic{3,1}=diag([Ic3x,Ic3y,Ic3z]);% Ic3; Ic4];
    syms gx gy gz real;
    g0_dh= [ 0; -gy; 0;];
    % Step(1) Assigning the DH-frames AND irci refereed to the same frame i
    syms rc1x rc1y rc1z rc2x rc2y rc2z rc3x rc3y rc3z rc4x rc4y rc4z real;
    rc1=[rc1x; 0; 0];
    rc2=[0; rc2y; 0];
    rc3=[rc3x; 0; 0];
    %rc4=[rc4x; rc4y; rc4z];
    rci_Matrix=[rc1, rc2, rc3];%, rc3, rc4];
    % Step(2) Optaining the DH-table:
    % alpha a d theta; alpha2 a2 .. ;
    DH_Table= [0 L1 0 q1; pi/2 0 L2 q2;0 L3 0 q3];
%--------------------------------------------------------------------------
%-----------------------------initializations------------------------------
    W = cell(n+1,1);
    W{1,1}= [0;0;0];
    V = cell(n+1,1);
    V{1,1}= [0;0;0];
    Vc = cell(n,1);
    U_link= cell(n,1);
    T_link= cell(n,1);
    U_tot=0; 
    T_tot=0;  
    M= cell(1,1);
    C= cell(n,1);
    c= cell(1,1);
    zero_rci_Matrix= cell(n,1);
%---------------------------------------------------------
% Step(3) Gompute the DH-homo-transformation matrix for each frame Ai:
    [Te,A]= DHMatrix(DH_Table);
for i=1:1:n
    %-----------------------------------------------------
    % Step(4) Compute 0_rci from i_rci for potential energy latter on:
    zero_rci_Matrix{i,1}= [rci_Matrix(:,i);1];
    for j=i:-1:1
    zero_rci_Matrix{i,1}= A{j}*zero_rci_Matrix{i,1};  
    end
    %-----------------------------------------------------
    % Step(5) Compute angular velocities and linear velocities: 
    RFi=transpose(A{i}(1:3, 1:3));
    W{i+1,1}= RFi*(W{i,1} + (1-k(i))*dq(i)*[0; 0;1]);
    V{i+1,1}= RFi*(V{i,1} + k(i)*dq(i)*[0;0;1]) + cross(W{i+1,1},RFi*A{i}(1:3,4));
    Vc{i,1}= V{i+1,1} + cross(W{i+1,1},rci_Matrix(:,i));
    %-----------------------------------------------------
    % Step(6) Compute Link energies and incrementing the total energy: 
    T_link{i,1} = simplify(0.5 *transpose(Vc{i,1})*m(i)*Vc{i,1}+  0.5*transpose(W{i+1,1})*Ic{i,1}*W{i+1,1});
    T_tot= simplify(T_tot+T_link{i,1});
    U_link{i,1}= simplify(-m(i)*transpose(g0_dh)*zero_rci_Matrix{i,1}(1:3,1));
    U_tot= simplify(U_tot+U_link{i,1});
    %-----------------------------------------------------
end
% Step(7) Compute the total energy/Hamiltonian and Lagrangian function: 
    T_tot = collect(T_tot, 1/2);
    T_tot = collect(T_tot, m);
    T_tot = collect(T_tot, dq);
    T_tot = simplify(T_tot); 
    H = simplify(T_tot+U_tot);
    L = simplify(T_tot-U_tot);
%---------------------------------------------------------
% Step(8) Compute the Inertia Matrix M(q): 
for i=1:1:n
    Ti= diff(T_tot,dq(i));
    M{1,1}(i,i)= diff(T_tot,dq(i),2);
    for j=i+1:1:n
        M{1,1}(i,j)=diff(Ti,dq(j));
        M{1,1}(j,i)=M{1,1}(i,j); %symetric.
    end
end
M{1,1}=simplify(M{1,1});
 %-------------------------------------------------------
 %Step(9) compute the Christoffel Matrix Cell aray with axis C{i,1}=Ci: 
for i=1:1:n
   M1=M{1,1}(:,i);
   JMq=jacobian(M1,q);
   C{i,1}=(1/2)*(JMq+transpose(JMq)-diff(M,q(i)));
   C{i,1}=simplify(C{i,1});
   c{1,1}(i,1)= transpose(dq)*C{i,1}*dq; % Coriolis Vector c(q, dq);
   c{1,1}(i,1)= simplify(c{1,1}(i,1));
end
%--------------------------------------------------------
 %Step(10) compute the gravity vector g(q): 
     g=transpose(jacobian(U_tot,q));
     g=simplify(g);
 %-------------------------------------------------------
 %Step(11) Displaying information:
 disp('----------------------------------------The Inertia Matrix M(q):-------------------- ');
     disp(M{1,1});
     disp('');
 disp('----------------------------------------The Coriolis Vector C(q,dq):---------------- ');
     disp(c{1,1});
     disp('');
 disp('----------------------------------------The Gravity  Vector g(q):------------------- ');
     disp(g);
 disp('------------------------------------------------------------------------------------ ');

 %--------------------END--------------------------------------------------
 %T_tot = collect(T_tot, dq(i)^2);
% % Subs value 
% F (1st col=f1 and so on) = simplify(subs([f1 f2 ..],{q1,q2,q3},{pi/6,pi/6,pi/6})); double(J);
% simplify(norm(qdd))
