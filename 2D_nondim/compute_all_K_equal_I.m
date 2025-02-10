function [x,y, c_total,m_total,n_total,h_total] = compute_all_K_equal_I()
clear all;
clc;
close all;
% This is the function, which calls all the related function and
% simulates the MP model (4.4.39 of the thesis) in 2-D(non-dimensional form) 
% and saves the involved variables
% in this case Kij’s equal to K_cm: Experiment 4.4 of the thesis (model 4.4.39)
%% x,y and t domain in dimensional form
dx_dim = 5;  % spatial step length
x_dim = 0:dx_dim:1000; % x domain
y_dim = x_dim; % y,same as x for a square domain
dt_dim = 0.001; % time spacing
t_end_dim = 15; % end time %
t_dim = 0:dt_dim:t_end_dim; % array of time steps
%% Parameters (dimensional)
kcm_dim = 0.042;
kmn_dim = 0.042;
kcn_dim = 0.042;
alpha_c_dim = 1;
alpha_m_dim = 1;
chi_dim = 2e+7;
theta = 1;
scale = 60*60*24*30;
b1_dim = 1e-9*scale;
b2_dim = 1e-11*scale;
D_h_dim = 10000;
h_max = 10^(-6.4);
c1_dim = (0.2/(60*60*24))*scale;
c2_dim = 1.5*c1_dim;
c3_dim = 0.01*c2_dim;
%% Parameters (non-dimensional),
x_scale = sqrt(b2_dim/D_h_dim);
dt = dt_dim*b2_dim;
t = t_dim*b2_dim;
dx = dx_dim*x_scale;
x = x_dim*x_scale;
y = x;
c1 = c1_dim/b2_dim;
c2 = c2_dim/b2_dim;
c3 = c3_dim/b2_dim;
kmn = kmn_dim*(D_h_dim/alpha_c_dim);
kcn = kcn_dim*(D_h_dim/alpha_c_dim);
kcm = kcm_dim*(D_h_dim/alpha_c_dim);
chi = chi_dim*h_max/alpha_c_dim;
alpha_m = alpha_m_dim/alpha_c_dim;
b1 = b1_dim/(b2_dim*h_max);
D_h = 1;

%% Grid generation and memory allocations
[X,Y] = meshgrid(x,y);
Lx = length(x);
Lt = length(t);
N = (Lx-2)^2; % number of unknowns in the domain
adv_c_x = zeros(Lx,Lx); % Glioma advection term x-component
adv_c_y = zeros(Lx,Lx); % Glioma advection term y-component
adv_m_x = zeros(Lx,Lx); % Normal cells advection term x-component
adv_m_y = zeros(Lx,Lx); % Normal cells advection term y-component
c_new = zeros(Lx,Lx); %
m_new = zeros(Lx,Lx);
% these  lines to save 10 values of solutions for 1 time interval(month/week/day)
temp = 1/(10*dt_dim);
c_total = zeros(Lx,Lx,t_end_dim*10+1);
m_total = zeros(Lx,Lx,t_end_dim*10+1);
n_total = zeros(Lx,Lx,t_end_dim*10+1);
h_total = zeros(Lx,Lx,t_end_dim*10+1);

%% initial conditions
centerx = 500*x_scale; centerx2 = 600*x_scale; centerx3 =300*x_scale; %centerx3 = 500;
centery = 500*x_scale; centery2 = 500*x_scale; centery3 =400*x_scale; %centery3 = 600;

% for tumor
sigma1  = 25*x_scale;   sigma2  = 20*x_scale;   sigma3  = 10*x_scale;
exponent = ((X-centerx).^2 + (Y-centery).^2)./(2*sigma1^2);
exponent2 = ((X-centerx2).^2 + (Y-centery2).^2)./(2*sigma2^2);
exponent3 = ((X-centerx3).^2 + (Y-centery3).^2)./(2*sigma3^2);
c_old       =  0.05*(exp(-exponent)+exp(-exponent2)+exp(-exponent3));
% c_old       =  0.05*(exp(-exponent)+exp(-exponent2));
c_old = c_old';

% for necritic cells
sigma1  = 5*x_scale;   sigma2  = 2*x_scale;   sigma3  = 1*x_scale;
exponent = ((X-centerx).^2 + (Y-centery).^2)./(2*sigma1^2);
exponent2 = ((X-centerx2).^2 + (Y-centery2).^2)./(2*sigma2^2);
exponent3 = ((X-centerx3).^2 + (Y-centery3).^2)./(2*sigma3^2);
n_old       =  0.9*(exp(-exponent)+exp(-exponent2)+exp(-exponent3));
% n_old       =  0.95*(exp(-exponent)+exp(-exponent2));
n_old = n_old';

% for normal cells
m_old = 1 - c_old - n_old;

radius = 2.5*x_scale;   radius2 = 1*x_scale;   radius3 = 1*x_scale;

for i=1:length(x)
    for j=1:length(y)
        dx1 =   (x(i)-centerx).^2;
        dy1 =   (y(j)-centery).^2;
        dx2 =  (x(i)-centerx2).^2;
        dy2=   (y(j)-centery2).^2;
        dx3 =  (x(i)-centerx3).^2;
        dy3 =  (y(j)-centery3).^2;
        distance = sqrt(dx1+dy1);
        distance2 = sqrt(dx2+dy2);
        distance3 = sqrt(dx3+dy3);
        
        
        if (distance<radius)
            m_old(i,j) = 0;
        elseif (distance2<radius2)
            m_old(i,j) = 0;
        elseif (distance3<radius3)
            m_old(i,j) = 0;
        end
        
    end
end

% for acidity
sigma1  = 15*x_scale;   sigma2  = 10*x_scale;   sigma3  = 7.5*x_scale;
exponent = ((X-centerx).^2 + (Y-centery).^2)./(2*sigma1^2);
exponent2 = ((X-centerx2).^2 + (Y-centery2).^2)./(2*sigma2^2);
exponent3 = ((X-centerx3).^2 + (Y-centery3).^2)./(2*sigma3^2);
h_old       =  (1/h_max)*(10^(-7)*exp(-exponent)+10^(-7)*exp(-exponent2)+10^(-6.4)*exp(-exponent3));
% h_old       =  exp(-exponent)+ exp(-exponent2);
h_old = h_old';

%% save the copy of initial cond
c_total (:,:,1) = c_old;
m_total (:,:,1) = m_old;
n_total (:,:,1) = n_old;
h_total (:,:,1) = h_old;

%% calling the diffusion matrix functions
A_H = set_const_diff(x,D_h,dt); % diff matrix for Acidity
%%

count = 2;
count2 = 1;
for j = 2: length(t)
    % the coefficients
    D11 = (kcm*c_old + kmn*(1 - c_old)).*(1 - c_old) - (kcm - kcn)*c_old.*m_old;
    D12 = -(kcm*c_old + kmn*(1 - c_old)).*c_old + (kcm - kcn)*c_old.*(1-m_old);
    D21 = -(kcm*m_old + kcn*(1 - m_old)).*m_old + (kcm - kmn)*m_old.*(1-c_old);
    D22 = (kcm*m_old + kcn*(1 - m_old)).*(1 - m_old) - (kcm - kmn)*c_old.*m_old;
    
    % S(c,m)
    S = kcm*kcn*c_old + kcm*kmn*m_old + kcn*kmn*(1-c_old-m_old);
    
    % G(h)
    G_h = ((chi*h_old)./(1 + (h_old)));
    % diff coeff for c
    D_C = (D11.*(2*c_old + G_h) + theta.*alpha_m.*D12.*m_old.*m_old)./S;
    
    % 5 stencils (explicit) for C
    D_C1 = (0.5/(dx*dx))*(D_C(1:end-2,2:end-1) + D_C(2:end-1,2:end-1));
    D_C2 = -(0.5/(dx*dx))*(4*D_C(2:end-1,2:end-1) + D_C(3:end,2:end-1)...
        + D_C(1:end-2,2:end-1) + D_C(2:end-1,1:end-2) + D_C(2:end-1,3:end));
    D_C3 = (0.5/(dx*dx))*(D_C(3:end,2:end-1) + D_C(2:end-1,2:end-1));
    D_C4 = (0.5/(dx*dx))*(D_C(2:end-1,1:end-2) + D_C(2:end-1,2:end-1));
    D_C5 = (0.5/(dx*dx))*(D_C(2:end-1,2:end-1) + D_C(2:end-1,3:end));
    
    % discretization of diff part
    diff_temp_C = dt*(D_C1.*c_old(1:end-2,2:end-1))...
        + (1 + dt*D_C2).*c_old(2:end-1,2:end-1)...
        + dt*D_C3.*c_old(3:end,2:end-1)...
        + dt*D_C4.*c_old(2:end-1,1:end-2)...
        + dt*D_C5.*c_old(2:end-1,3:end);
    
    %diff coeff for m
    D_M = (D22*alpha_m.*(1 + theta*c_old).*(2*m_old))./S;
    
    % 5 stencils (explicit) for m
    D_M1 = (0.5/(dx*dx))*(D_M(1:end-2,2:end-1) + D_M(2:end-1,2:end-1));
    D_M2 = -(0.5/(dx*dx))*(4*D_M(2:end-1,2:end-1) + D_M(3:end,2:end-1) ...
        + D_M(1:end-2,2:end-1) + D_M(2:end-1,1:end-2) + D_M(2:end-1,3:end));
    D_M3 = (0.5/(dx*dx))*(D_M(3:end,2:end-1) + D_M(2:end-1,2:end-1));
    D_M4 = (0.5/(dx*dx))*(D_M(2:end-1,1:end-2) + D_M(2:end-1,2:end-1));
    D_M5 = (0.5/(dx*dx))*(D_M(2:end-1,2:end-1) + D_M(2:end-1,3:end));
    
    % discretization of diff part
    diff_temp_M = dt*(D_M1.*m_old(1:end-2,2:end-1))...
        + (1+dt*D_M2).*m_old(2:end-1,2:end-1)...
        + dt*(D_M3.*m_old(3:end,2:end-1))...
        + dt*(D_M4.*m_old(2:end-1,1:end-2))...
        + dt*(D_M5.*m_old(2:end-1,3:end));
    
    % advective flux for c in x-direction
    adv_c_x(2:end-1,2:end-1)= (D11(2:end-1,2:end-1).*c_old(2:end-1,2:end-1).*((G_h(3:end,2:end-1)-G_h(1:end-2,2:end-1))./(2*dx))...
        + D12(2:end-1,2:end-1)*alpha_m.*(1+theta*c_old(2:end-1,2:end-1))...
        .*(2*m_old(2:end-1,2:end-1)).*((m_old(3:end,2:end-1)-m_old(1:end-2,2:end-1))./(2*dx)))./S(2:end-1,2:end-1);
    
    % advective flux for c in y-direction
    adv_c_y(2:end-1,2:end-1) = (D11(2:end-1,2:end-1).*c_old(2:end-1,2:end-1).*((G_h(2:end-1,3:end)-G_h(2:end-1,1:end-2))./(2*dx))...
        + D12(2:end-1,2:end-1)*alpha_m.*(1+theta*c_old(2:end-1,2:end-1))...
        .*(2*m_old(2:end-1,2:end-1)).*((m_old(2:end-1,3:end)-m_old(2:end-1,1:end-2))./(2*dx)))./S(2:end-1,2:end-1);
    
    % \tau_c * c
    temp_m1 = (c_old + G_h).*c_old;
    
    % advective flux for m in x-direction
    adv_m_x(2:end-1,2:end-1)= (D21(2:end-1,2:end-1).*((temp_m1(3:end,2:end-1) - temp_m1(1:end-2,2:end-1))./(2*dx)) ...
        + D22(2:end-1,2:end-1)*theta*alpha_m.*m_old(2:end-1,2:end-1).*m_old(2:end-1,2:end-1)...
        .*((c_old(3:end,2:end-1)-c_old(1:end-2,2:end-1))./(2*dx)))./S(2:end-1,2:end-1);
    
    % advective flux for m in y-direction
    adv_m_y(2:end-1,2:end-1) = (D21(2:end-1,2:end-1).*((temp_m1(2:end-1,3:end) - temp_m1(2:end-1,1:end-2))./(2*dx)) ...
        + D22(2:end-1,2:end-1)*theta*alpha_m.*m_old(2:end-1,2:end-1).*m_old(2:end-1,2:end-1)...
        .*((c_old(2:end-1,3:end)-c_old(2:end-1,1:end-2))./(2*dx)))./S(2:end-1,2:end-1);
    
    % sign test for upwind
    test_up_c_x = (((kcn*kmn - kcm*kcn)./(S(2:end-1,2:end-1).^2)).*D11(2:end-1,2:end-1).*c_old(2:end-1,2:end-1)...
        + (1./S(2:end-1,2:end-1)).*c_old(2:end-1,2:end-1).*(kcm*(1-2*c_old(2:end-1,2:end-1)-m_old(2:end-1,2:end-1))...
        - kmn*2*(1-c_old(2:end-1,2:end-1)) + kcn*m_old(2:end-1,2:end-1)) + ...
        D11(2:end-1,2:end-1)./S(2:end-1,2:end-1)).*((G_h(3:end,2:end-1)-G_h(1:end-2,2:end-1))./(2*dx)) + ...
        (2*alpha_m*m_old(2:end-1,2:end-1)*((m_old(3:end,2:end-1)-m_old(1:end-2,2:end-1))./(2*dx))).*...
        ((1./S(2:end-1,2:end-1)).*(1+theta*c_old(2:end-1,2:end-1)).*(kcm*(1-m_old(2:end-1,2:end-1)-2*c_old(2:end-1,2:end-1)) ...
        + kmn*(2*c_old(2:end-1,2:end-1)-1) + kcn*(m_old(2:end-1,2:end-1)-1) ) + (D12(2:end-1,2:end-1)*theta./S(2:end-1,2:end-1)) ...
        + D12(2:end-1,2:end-1).*(1+theta*c_old(2:end-1,2:end-1))*((kcn*kmn - kcm*kcn)./(S(2:end-1,2:end-1).^2)));
    
    test_up_c_y = (((kcn*kmn - kcm*kcn)./(S(2:end-1,2:end-1).^2)).*D11(2:end-1,2:end-1).*c_old(2:end-1,2:end-1)...
        + (1./S(2:end-1,2:end-1)).*c_old(2:end-1,2:end-1).*(kcm*(1-2*c_old(2:end-1,2:end-1)-m_old(2:end-1,2:end-1))...
        - kmn*2*(1-c_old(2:end-1,2:end-1)) + kcn*m_old(2:end-1,2:end-1)) + ...
        D11(2:end-1,2:end-1)./S(2:end-1,2:end-1)).*((G_h(2:end-1,3:end)-G_h(2:end-1,1:end-2))./(2*dx)) + ...
        (2*alpha_m*m_old(2:end-1,2:end-1)*((m_old(2:end-1,3:end)-m_old(2:end-1,1:end-2))./(2*dx))).*...
        ((1./S(2:end-1,2:end-1)).*(1+theta*c_old(2:end-1,2:end-1)).*(kcm*(1-m_old(2:end-1,2:end-1)-2*c_old(2:end-1,2:end-1))...
        + kmn*(2*c_old(2:end-1,2:end-1)-1) + kcn*(m_old(2:end-1,2:end-1)-1) ) + (D12(2:end-1,2:end-1)*theta./S(2:end-1,2:end-1)) ...
        + D12(2:end-1,2:end-1).*(1+theta*c_old(2:end-1,2:end-1))*((kcn*kmn - kcm*kcn)./(S(2:end-1,2:end-1).^2)));
    
    test_up_m_x = ((temp_m1(3:end,2:end-1) - temp_m1(1:end-2,2:end-1))./(2*dx)).*...
        ( (D21(2:end-1,2:end-1).*((kcn*kmn - kcm*kmn)./(S(2:end-1,2:end-1).^2)))...
        + (1./S(2:end-1,2:end-1)).*( kcm*(1-c_old(2:end-1,2:end-1) -2*m_old(2:end-1,2:end-1))...
        + kcn*(2*m_old(2:end-1,2:end-1) -1) + kmn*(c_old(2:end-1,2:end-1)-1))) + ...
        (alpha_m*theta*((c_old(3:end,2:end-1) - c_old(1:end-2,2:end-1))./(2*dx)).*...
        ((D22(2:end-1,2:end-1).*(m_old(2:end-1,2:end-1).^2).*((kcn*kmn - kcm*kmn)./(S(2:end-1,2:end-1).^2))) ...
        + (1./S(2:end-1,2:end-1)).*(2*D22(2:end-1,2:end-1).*m_old(2:end-1,2:end-1)...
        + (m_old(2:end-1,2:end-1).^2).*(kcm*(1-2*m_old(2:end-1,2:end-1)-c_old(2:end-1,2:end-1))...
        - kcn*(2-2*m_old(2:end-1,2:end-1)) + kmn*c_old(2:end-1,2:end-1) ))  )  );
    
    test_up_m_y = ((temp_m1(2:end-1,3:end) - temp_m1(2:end-1,1:end-2))./(2*dx)).*...
        ( (D21(2:end-1,2:end-1).*((kcn*kmn - kcm*kmn)./(S(2:end-1,2:end-1).^2)))...
        + (1./S(2:end-1,2:end-1)).*( kcm*(1-c_old(2:end-1,2:end-1) -2*m_old(2:end-1,2:end-1))...
        + kcn*(2*m_old(2:end-1,2:end-1) -1) + kmn*(c_old(2:end-1,2:end-1)-1))) + ...
        (alpha_m*theta*((c_old(2:end-1,3:end) - c_old(2:end-1,1:end-2))./(2*dx)).*...
        ((D22(2:end-1,2:end-1).*(m_old(2:end-1,2:end-1).^2).*((kcn*kmn - kcm*kmn)./(S(2:end-1,2:end-1).^2)))+ ...
        (1./S(2:end-1,2:end-1)).*(2*D22(2:end-1,2:end-1).*m_old(2:end-1,2:end-1)...
        + (m_old(2:end-1,2:end-1).^2).*(kcm*(1-2*m_old(2:end-1,2:end-1)-c_old(2:end-1,2:end-1))...
        - kcn*(2-2*m_old(2:end-1,2:end-1)) + kmn*c_old(2:end-1,2:end-1) ))  )  );
    
    %     test_c_x = (adv_c_x(3:end,2:end-1) - adv_c_x(2:end-1,2:end-1))./(c_old(3:end,2:end-1) - c_old(2:end-1,2:end-1));
    %     test_c_y = (adv_c_y(2:end-1,3:end) - adv_c_y(2:end-1,2:end-1))./(c_old(2:end-1,3:end) - c_old(2:end-1,2:end-1));
    %     test_m_x = (adv_m_x(3:end,2:end-1) - adv_m_x(2:end-1,2:end-1))./(m_old(3:end,2:end-1) - m_old(2:end-1,2:end-1));
    %     test_m_y = (adv_m_y(2:end-1,3:end) - adv_m_x(2:end-1,2:end-1))./(m_old(2:end-1,3:end) - m_old(2:end-1,2:end-1));
    
    mask1 = test_up_c_x<=0;
    mask2 = test_up_c_y<=0;
    mask3 = test_up_m_x<=0;
    mask4 = test_up_m_y<=0;
    
    % upwinding in both directions
    advection_c_x = (max(mask1,0)).*((adv_c_x(2:end-1,2:end-1)- adv_c_x(1:end-2,2:end-1))/dx) + ...
        (-min(mask1-1,0)).*((adv_c_x(3:end,2:end-1)- adv_c_x(2:end-1,2:end-1))/dx);
    
    advection_c_y = (max(mask2,0)).*((adv_c_y(2:end-1,2:end-1)- adv_c_y(2:end-1,1:end-2))/dx) + ...
        (-min(mask2-1,0)).*((adv_c_y(2:end-1,3:end)- adv_c_y(2:end-1,2:end-1))/dx);
    
    advection_m_x = (max(mask3,0)).*((adv_m_x(2:end-1,2:end-1)- adv_m_x(1:end-2,2:end-1))/dx) + ...
        (-min(mask3-1,0)).* ((adv_m_x(3:end,2:end-1)- adv_m_x(2:end-1,2:end-1))/dx);
    
    advection_m_y = (max(mask4,0)).*((adv_m_y(2:end-1,2:end-1)- adv_m_y(2:end-1,1:end-2))/dx) + ...
        (-min(mask4-1,0)).*((adv_m_y(2:end-1,3:end)- adv_m_y(2:end-1,2:end-1))/dx);
    
    
    
    %source for c
    source_c = c1*c_old.*n_old.*((h_old)-1);%source
    
    %source for m
    %     source_m = (c2*(max(0,(h_old)-1)).*m_old)./(1+m_old+(h_old));
    %     source_m = c2*m_old.*(1 - m_old).*((h_old)-1) + c3*m_old;%source
    source_m = (c2*(max(0,(h_old)-1)).*m_old)./(1+m_old+(h_old))  + c3*m_old;%source
    
    % glioma(C) cell update (explicit)
    c_new(2:end-1,2:end-1) = diff_temp_C + dt*(advection_c_x + advection_c_y) - dt*source_c(2:end-1,2:end-1);
    
    % normal cell update (explicit)
    m_new(2:end-1,2:end-1) = diff_temp_M + dt*(advection_m_x + advection_m_y) - dt*source_m(2:end-1,2:end-1);
    
    B_C = reshape(c_old(2:end-1,2:end-1),N,1);
    B_H = reshape(h_old(2:end-1,2:end-1),N,1);
    
    RHS_H = B_H + dt*b1*B_C - dt*B_H;% RHS containing source & uptake for acidity(h) eq
    
    [sol_H,flag] = bicgstab(A_H,RHS_H); % solution for acidity(h) at interior nodes
    
    %% boundary conditions
    
    c_new(1,:) = c_old(1,:) + dt*(((D_C(1,:) + D_C(2,:)).*(c_old(2,:) - c_old(1,:))*(1/(dx*dx)))...
        - source_c(1,:));
    
    c_new(end,:) = c_old(end,:) + dt*((D_C(end,:) + D_C(end-1,:)).*(c_old(end-1,:) - c_old(end,:))*(1/(dx*dx))...
        - source_c(end,:));
    
    c_new(:,1) = c_old(:,1) + dt*(((D_C(:,1) + D_C(:,2)).*(c_old(:,2) - c_old(:,1))*(1/(dx*dx)))...
        - source_c(:,1));
    
    c_new(:,end) = c_old(:,end) + dt*((D_C(:,end) + D_C(:,end-1)).*(c_old(:,end-1) - c_old(:,end))*(1/(dx*dx))...
        - source_c(:,end));
    
    m_new(1,:) = m_old(1,:) + dt*(((D_M(1,:) + D_M(2,:)).*(m_old(2,:) - m_old(1,:))*(1/(dx*dx)))...
        - source_m(1,:));
    
    m_new(end,:) = m_old(end,:) + dt*((D_M(end,:) + D_M(end-1,:)).*(m_old(end-1,:) - m_old(end,:))*(1/(dx*dx))...
        - source_m(end,:));
    
    m_new(:,1) = m_old(:,1) + dt*(((D_M(:,1) + D_M(:,2)).*(m_old(:,2) - m_old(:,1))*(1/(dx*dx)))...
        - source_m(:,1));
    
    m_new(:,end) = m_old(:,end) + dt*((D_M(:,end) + D_M(:,end-1)).*(m_old(:,end-1) - m_old(:,end))*(1/(dx*dx))...
        - source_m(:,end));
    
    H_sol = reshape(sol_H,Lx-2,Lx-2);
    new_H = [H_sol(1,1:end);H_sol;H_sol(end,1:end)];
    h_new = [new_H(1:end,1),new_H, new_H(1:end,end)];
    
    % nectoric cells
    n_new = 1 - c_new - m_new;
    
    if (mod(count2,temp)==0)
        c_total (:,:,count) = c_new;
        m_total (:,:,count) = m_new;
        n_total (:,:,count) = n_new;
        h_total (:,:,count) = h_new;
        
        count = count+1;
        
    end
    % resetting new solution as old for time marching
    c_old = c_new;
    m_old = m_new;
    n_old = n_new;
    h_old = h_new;
    
    
    count2 = count2+1;
end

end