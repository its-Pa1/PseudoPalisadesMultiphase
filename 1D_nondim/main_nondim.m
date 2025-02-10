% This is the main script, which calls all the related function and
% computes the MP model (base case model 4.4.39, Experiment 1) in 1-D and
% saves the results in a video.
% The corresponding reults have been included in Figure 4.10 of the thesis.

clear all;
clc;
close all;
tic
%% Grid generation, Memory allocation
dx_dim = 5;
x_dim = 0:dx_dim:1000;
dt_dim = 0.001;
t_dim = 0:dt_dim:15;
%% Parameters (dimensional)
kcm_dim = 0.042;
kmn_dim = 0.043;
kcn_dim = 0.044;
alpha_c_dim = 1;
alpha_m_dim = 1;
chi_dim = 2e+7;
theta = 1;
scale = 60*60*24*30;
b1_dim = 1e-9*scale;
b2 = 1e-11*scale;
D_h_dim = 10000;
h_max = 10^(-6.4);
c1_dim = (0.2/(60*60*24))*scale;
c2_dim = 1.5*c1_dim;
c3_dim = (1/100)*c2_dim;
%% Parameters (non-dimensional)
x_scale = sqrt(b2/D_h_dim);
dt = dt_dim*b2;
t = t_dim*b2;
dx = dx_dim*x_scale;
x = x_dim*x_scale;
Lx = length(x);
Lt = length(t);
temp = 1/(10*dt_dim);
c_new = zeros(1,Lx);
m_new = zeros(1,Lx);
c_total = zeros(Lx,Lt);
m_total = zeros(Lx,Lt);
n_total = zeros(Lx,Lt);
h_total = zeros(Lx,Lt);
diff_c = zeros(size(x));
adv_c = zeros(size(x));
adv_m = zeros(size(x));
velocity_c = zeros(size(x));

c1 = c1_dim/b2;
c2 = c2_dim/b2;
c3 = c3_dim/b2;
kmn = kmn_dim*(D_h_dim/alpha_c_dim);
kcn = kcn_dim*(D_h_dim/alpha_c_dim);
kcm = kcm_dim*(D_h_dim/alpha_c_dim);
chi = chi_dim*h_max/alpha_c_dim;
alpha_m = alpha_m_dim/alpha_c_dim;
b1 = b1_dim/(b2*h_max);
D_h = 1;
%% initial conditions
sigma_c  = 25*x_scale;    sigma_n = 5*x_scale;  sigma_h  = 10*x_scale;
center_c = 500*x_scale;     center_n = 500*x_scale;  center_h = 500*x_scale;
exponent_c = ((x-center_c).^2)./(2*sigma_c^2);
exponent_n = ((x-center_n).^2)./(2*sigma_n^2);
exponent_h = ((x-center_h).^2)./(2*sigma_h^2);
c_old = 0.0025*exp(-exponent_c);
n_old = 0.9*exp(-exponent_n);
m_old = 1 - (c_old + n_old); % 1-(c+n)
m_old(101) = 0; % normal tissue 0 at x=0
h_old = exp(-exponent_h);
%% save the copy of initial cond
c_total(:,1) = c_old;
m_total(:,1) = m_old;
n_total(:,1) = n_old;
h_total(:,1) = h_old;
%%
% acidity diff matrix (for IMEX method)
A2  = diag((1) + (dt*2*D_h/(dx*dx)) * ones(Lx,1), 0) + diag(-dt*D_h/(dx*dx)* ones(Lx-1,1), -1) ...
    + diag(- (dt*D_h/(dx*dx))* ones(Lx-1,1), 1);
A2(1,2) = A2(1,2) - dt*D_h/(dx*dx); % BC in diff matrix
A2(end,end-1) = A2(end,end-1) - dt*((D_h/(dx*dx))); %BC in diff matrix
B_h = zeros(1,size(A2,2)); % RHS for acidity

%% time loop
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
    diff_c = (D11.*(2*c_old + G_h) + theta.*alpha_m.*D12.*m_old.*m_old)./S;
    
    % advective flux for c
    adv_c(2:end-1) = (D11(2:end-1).*c_old(2:end-1).*((G_h(3:end)-G_h(1:end-2))./(2*dx))...
        + D12(2:end-1)*alpha_m.*(1+theta*c_old(2:end-1))...
        .*(2*m_old(2:end-1)).*((m_old(3:end)-m_old(1:end-2))./(2*dx)))./S(2:end-1);
    
    
    %diff coeff for m
    diff_m = (D22*alpha_m.*(1 + theta*c_old).*(2*m_old))./S;
    
    % \tau_c * c
    temp_m1 = (c_old + G_h).*c_old;
    
    % advective flux for m
    adv_m(2:end-1) = (D21(2:end-1).*((temp_m1(3:end) - temp_m1(1:end-2))./(2*dx)) ...
        + D22(2:end-1)*theta*alpha_m.*m_old(2:end-1).*m_old(2:end-1).*((c_old(3:end)-c_old(1:end-2))./(2*dx)))./S(2:end-1);
    
    %source for c
    source_c = c1*c_old.*n_old.*(h_old-1);%source
    
    %source for m
    %     source_m = (c2*(max(0,h_old - 1)).*m_old)./(1 + m_old + h_old );
    %     source_m = c2*m_old.*(1 - m_old).*((h_old)-1) + c3*m_old;%source
    source_m = (c2*(max(0,h_old - 1)).*m_old)./(1 + m_old + h_old ) + c3*m_old;%source
    
    
    
    for i = 2:Lx-1
        
        %%%%%%%% condition for conservative first order upwind   %%%%%%%%%%%%%
        if(c_old(i)==c_old(i+1))
            temp_up_c = (((1/(2*dx))*(G_h(i+1)-G_h(i-1)))*(D11(i) +...
                c_old(i)*(2*(c_old(i)-1)*kmn + kcn*m_old(i) - kcm*(m_old(i)+2*c_old(i)-1))))...
                + (((2*c_old(i)-1)*kmn + kcn*(m_old(i)-1)-kcm*(-1+2*c_old(i) + m_old(i)))*alpha_m...
                *(1+theta*c_old(i)) + D12(i)*alpha_m*theta)*2*m_old(i)*((1/(2*dx))*(m_old(i+1)-m_old(i-1)))...
                + adv_c(i)*(kcn*kmn - kcm*kcn);
            
        else
            temp_up_c = (adv_c(i+1)- adv_c(i))/(c_old(i+1)-c_old(i));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % explicit scheme for tumor cells(c) on internal nodes
        if (temp_up_c<=0)
            c_new(i) = c_old(i) ...
                + dt*(((0.5/(dx*dx))*(((diff_c(i)+ diff_c(i+1))*(c_old(i+1)-c_old(i)))...% diffusion
                -((diff_c(i-1)+diff_c(i))*(c_old(i) -c_old(i-1)))))...% diffusion till here
                + (adv_c(i)-adv_c(i-1))*(1/dx)... % upwind for taxis
                - source_c(i));%source
            
        else
            c_new(i) = c_old(i) ...
                + dt*(((0.5/(dx*dx))*(((diff_c(i)+ diff_c(i+1))*(c_old(i+1)-c_old(i)))...% diffusion
                -((diff_c(i-1)+diff_c(i))*(c_old(i) -c_old(i-1)))))...% diffusion till here
                + (adv_c(i+1)-adv_c(i))*(1/dx)... upwind for taxis
                - source_c(i));%source
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % condition for conservative first order upwind for normal cells, m
        if(m_old(i)==m_old(i+1))
            temp_up_m = ((c_old(i)-1)*kmn + kcn*(2*m_old(i)-1) - kcm*(c_old(i)+2*m_old(i)-1))...
                *((1/(2*dx))*(temp_m1(i+1)-temp_m1(i-1))) ...
                + (alpha_m*((1/(2*dx))*(c_old(i+1)-c_old(i-1))))*(m_old(i)*m_old(i)*(c_old(i)*kmn...
                +2*kcn*(m_old(i)-1)-kcm*(c_old(i)+2*m_old(i)-1)) + D22(i)*2*m_old(i))...
                + adv_m(i)*(kmn*(kcm-kcn));
        else
            temp_up_m = (adv_m(i+1)- adv_m(i))/(m_old(i+1)-m_old(i));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % explicit scheme for normal cells(m)
        if(temp_up_m<=0)
            m_new(i) = m_old(i)...
                + dt*(((0.5/(dx*dx))*(((diff_m(i)+ diff_m(i+1))*(m_old(i+1)-m_old(i)))...% diffusion
                -((diff_m(i-1)+diff_m(i))*(m_old(i) -m_old(i-1)))))...% diffusion till here
                + (adv_m(i)-adv_m(i-1))*(1/dx)... % upwind for taxis
                - source_m(i)); % source
        else
            m_new(i) = m_old(i)...
                + dt*(((0.5/(dx*dx))*(((diff_m(i)+ diff_m(i+1))*(m_old(i+1)-m_old(i)))...% diffusion
                -((diff_m(i-1)+diff_m(i))*(m_old(i) -m_old(i-1)))))...% diffusion till here
                + (adv_m(i+1)-adv_m(i))*(1/dx)... % upwind for taxis
                - source_m(i)); % source
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %RHS for acidity (for IMEX scheme)
        B_h(i) = (h_old(i)) + dt*(b1*c_old(i) - h_old(i));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    
    % boundary conditions
    
    
    c_new(1) = c_old(1) + dt*(((diff_c(1) + diff_c(2))*(c_old(2) - c_old(1))*(1/(dx*dx)))...
        - source_c(1));
    
    c_new(end) = c_old(end) + dt*((diff_c(end) + diff_c(end-1))*(c_old(end-1) - c_old(end))*(1/(dx*dx))...
        - source_c(end));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    m_new(1) = m_old(1) + dt*(((diff_m(1) + diff_m(2))*(m_old(2) - m_old(1))*(1/(dx*dx)))...
        - source_m(1));
    
    m_new(end) = m_old(end) + dt*((diff_m(end) + diff_m(end-1))*(m_old(end-1) - m_old(end))*(1/(dx*dx))...
        - source_m(end));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % IMEX for h
    B_h(1) = (h_old(1)) + dt*(b1*c_old(1) - h_old(1));
    B_h(end) = (h_old(end)) + dt*(b1*c_old(end) - h_old(end));
    h_new = (A2\B_h')';
    
    
    % nectoric cells
    n_new = 1 - c_new - m_new;
    
    % saves only 10 solution per time step
    % if (mod(count2,temp)==0)
    c_total (:,count) = c_new;
    m_total (:,count) = m_new;
    n_total (:,count) = n_new;
    h_total (:,count) = h_new;
    count = count+1;
    
    % end
    
    count2 = count2+1;
    
    % update for time loop
    c_old = c_new;
    m_old = m_new;
    h_old = h_new;
    n_old = n_new;
    
end
%% Plots initial
% % subplot(2,2,1);
% % plot(x,c_new(:,1),'Linewidth', 1.5);
% % title('Tumor cells ')
% % subplot(2,2,2);
% % plot(x,m_new(:,1),'Linewidth', 1.5);
% % title('Normal cells ')
% % subplot(2,2,3);
% % plot(x,n_new(:,1),'Linewidth', 1.5);
% % title('Necrotic cells ')
% % subplot(2,2,4);
% % plot(x,h_new(:,1),'Linewidth', 1.5);
% % title('Acidity ')

%% Video% to check which video profile supports available in the machine
% if mp4 is not supported then avi format will be used
profiles = VideoWriter.getProfiles();
check_mp4_support = find(ismember({profiles.Name},'MPEG-4'));

if isempty(check_mp4_support)
    video_ext = '.avi';
    v_pro  = 'Motion JPEG AVI';
else
    video_ext = '.mp4';
    v_pro = 'MPEG-4';
end

videofile = VideoWriter(strcat('MPM_1D_nondim', video_ext),v_pro);
videofile.Quality = 100;
videofile.FrameRate = 5; % no of frame per sec
open(videofile);
%nn = 10; % no. of plots for 1 time interval required e.g 2 means 2 plot per month
%temp = 1/(nn*dt_dim);
for j = 1:temp:size(c_total,2)
    
    figure(1)
    subplot(2,2,1);
    plot(x,c_total(:,j),'Linewidth', 1.5);
    if (j==1)
        title('Initial glioma cells');
    else
        title(['Glioma cells at t = ', num2str((3/temp)*(j-1)), ' days '], 'Fontsize', 6);
    end
    axis tight
    drawnow
    subplot(2,2,2);
    plot(x,m_total(:,j),'Linewidth', 1.5);
    if (j==1)
        title('Initial normal tissue');
    else
        title(['Normal tissue at t = ', num2str((3/temp)*(j-1)), ' days '], 'Fontsize', 7);
    end
    axis tight
    drawnow
    
    subplot(2,2,3);
    plot(x,n_total(:,j),'Linewidth', 1.5);
    if (j==1)
        title('Initial necrotic matter');
    else
        title(['Necrotic matter at t = ', num2str((3/temp)*(j-1)), ' days '],'Fontsize', 5);
    end
    axis tight
    drawnow
    
    subplot(2,2,4);
    plot(x,h_total(:,j),'Linewidth', 1.5);
    if (j==1)
        title('Initial acidity');
    else
        title(['Acidity at t = ', num2str((3/temp)*(j-1)), ' days '],'Fontsize', 5);
    end
    axis tight
    drawnow
    
    F= getframe(gcf);
    writeVideo(videofile,F);
    
end
close(videofile)
%% Plots
% plots of glioma, Acidity, normal and necrotic cells with time and space

% checking if the directory exists, otherwise it creates one
if not(isfolder('Plots_pattern'))
    mkdir('Plots_pattern')
end

figure(2)
surf(x,t,(c_total)')
view(0,90)
colorbar
set(gca,'ColorScale','log') % in this case we use nonlinear color scale
%caxis([0 1e-5]) % as matlab colormap is only linear and has only few colors,
% a differnce of order 100 make eveything blue
% so this command colors red all the values above its second argument(e.g here 1e-5)
shading interp
colormap jet
xlabel('X' , 'Fontsize', 15);
ylabel('Time' , 'Fontsize', 15);
title('Glioma cells' , 'Fontsize', 15);
axis tight
% axis([280 440 0, t(end)]);
saveas(gcf,'Plots_pattern/Tumor_pattern_I.png');
%%
figure(3)
surf(x,t,m_total')
view(0,90)
colorbar
shading interp
colormap jet
xlabel('X' , 'Fontsize', 15);
ylabel('Time' , 'Fontsize', 15);
title('Normal tissue' , 'Fontsize', 15);
axis tight
% axis([280 440 0, t(end)]);
saveas(gcf,'Plots_pattern/Normal_cells_pattern_I.png');

figure(4)
surf(x,t,n_total')
view(0,90)
colorbar
shading interp
colormap jet
xlabel('X' , 'Fontsize', 15);
ylabel('Time' , 'Fontsize', 15);
title('Necrotic matter' , 'Fontsize', 15);
axis tight
% axis([280 440 0, t(end)]);
saveas(gcf,'Plots_pattern/Necrotic_pattern_I.png');
%%
figure(5)
surf(x,t,h_total')
view(0,90)
colorbar
% caxis([0 1.5]) % as matlab colormap has only few colors, a differnce of order 100 make eveything blue
% so this command colors red all the values above its second argument(e.g here 1.5)
shading interp
colormap jet
xlabel('X' , 'Fontsize', 15);
ylabel('Time' , 'Fontsize', 15);
title('Acidity' , 'Fontsize', 15);
axis tight
% axis([280 440 0, t(end)]);
saveas(gcf,'Plots_pattern/Acidity_pattern_I.png');

toc