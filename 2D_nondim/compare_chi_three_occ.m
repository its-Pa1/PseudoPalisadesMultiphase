% This is the main script, which calls all the related function and
% computes the MP model (4.4.39 of the thesis) in 2-D(non-dimensional form) 
% in three cases:
%      1. Base case (compute_three_occ()): Experiment 4.1 of the thesis (model 4.4.39)
%      2. with \chi = 0 (compute_three_occ_chi_zero()): Experiment 4.3 of the thesis (model 4.4.39)
%      3. with \chi is doubled than base case(compute_three_occ_chi_double()): Experiment 4.3 of the thesis (model 4.4.39)

% and saves the results in the folder Plots_diff_chi_zero, Plots_diff_chi_double,
% Videos_diff_chi_zero and Videos_diff_chi_double.
% The corresponding reults are included in Figures 4.5 and 4.6 of the thesis

clear all;
clc;
close all;
tic
[x,y,c,m,n,h] = compute_three_occ();
[x0,y0,c0,m0,n0,h0] = compute_three_occ_chi_zero();
[x2,y2,c2,m2,n2,h2] = compute_three_occ_chi_double();

c_diff0 = c - c0;
m_diff0 = m - m0;
n_diff0 = n - n0;
h_diff0 = h - h0;

c_diff2 = c - c2;
m_diff2 = m - m2;
n_diff2 = n - n2;
h_diff2 = h - h2;
%%
toc
fprintf('Computation Done! Saving the plots and videos...\n')
tic;
%% Videos
% In this section Videos are saved in Videos__diff_chi_zero folder
% To check if Videos_diff_chi_zero folder exists otherwise to create
if not(isfolder('Videos_diff_chi_zero'))
    mkdir('Videos_diff_chi_zero')
end

% to check which video profile supports available in the machine
% if mp4 is not supported then avi format will be used
profiles = VideoWriter.getProfiles();
check_mp4_support = find(ismember({profiles.Name},'MPEG-4'),1);

if isempty(check_mp4_support)
    video_ext = '.avi';
    v_pro  = 'Motion JPEG AVI';
else
    video_ext = '.mp4';
    v_pro = 'MPEG-4';
end

% These loops save the solution video  each 3 days
videofile = VideoWriter(strcat('Videos_diff_chi_zero/Tumor_diff_chi_zero', video_ext),v_pro);
videofile.Quality = 100;
videofile.FrameRate = 10;
open(videofile);
for i = 1: size(c_diff0,3)
    figure(4)
    surf(x,y,c_diff0(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    axis tight
    drawnow
    if (i==1)
        title('Initial glioma cells difference', 'Fontsize', 13);
    else
        title(['Glioma cells difference at t = ', num2str(3*(i-1)), ' days '], 'Fontsize', 13);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    
    F= getframe(gcf);
    writeVideo(videofile,F);
end
close(videofile);

videofile = VideoWriter(strcat('Videos_diff_chi_zero/Acidity_diff_chi_zero', video_ext),v_pro);
videofile.Quality = 100;
videofile.FrameRate = 10;
open(videofile);
for i = 1: size(h_diff0,3)
    figure(5)
    surf(x,y,h_diff0(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    axis tight
    drawnow
    if (i==1)
        title('Initial acidity difference', 'Fontsize', 13);
    else
        title(['Acidity difference at t = ', num2str(3*(i-1)),' days '], 'Fontsize', 13);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    axis tight
    F= getframe(gcf);
    writeVideo(videofile,F)
    
end
close(videofile);

videofile = VideoWriter(strcat('Videos_diff_chi_zero/Normal_diff_chi_zero', video_ext),v_pro);
videofile.Quality = 100;
videofile.FrameRate = 10;
open(videofile);
for i = 1: size(m_diff0,3)
    figure(6)
    surf(x,y,m_diff0(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    axis tight
    drawnow
    if (i==1)
        title('Initial normal tissue difference', 'Fontsize', 13);
    else
        title(['Normal tissue difference at t = ', num2str(3*(i-1)),' days '], 'Fontsize', 13);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    
    F= getframe(gcf);
    writeVideo(videofile,F)
    
end
close(videofile);

videofile = VideoWriter(strcat('Videos_diff_chi_zero/Necrotic_diff_chi_zero', video_ext),v_pro);
videofile.Quality = 100;
videofile.FrameRate = 10;
open(videofile);
for i = 1: size(n_diff0,3)
    figure(8)
    surf(x,y,n_diff0(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    axis tight
    drawnow
    if (i==1)
        title('Initial necrotic matter difference', 'Fontsize', 13);
    else
        title(['Necrotic matter difference at t = ', num2str(3*(i-1)),' days '], 'Fontsize', 13);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    
    F= getframe(gcf);
    writeVideo(videofile,F)
    
end
close(videofile);
%% % This loop saves the results for each 30 days in folder Plots__diff_chi_zero in png
% or eps format with days in file names

% to check wheather Plots__diff_chi_zero folder exists otherwise it makes a folder Plots__diff_chi_zero
if not(isfolder('Plots_diff_chi_zero'))
    mkdir('Plots_diff_chi_zero')
end

for i = 1:10:size(c_diff0,3)
    figure(9)
    surf(x,y,c_diff0(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    axis tight
    drawnow
    if (i==1)
        title('Initial glioma cells difference', 'Fontsize', 13);
    else
        title(['Glioma cells difference at t = ', num2str(3*(i-1)), ' days '], 'Fontsize', 13);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    
    saveas(gcf,sprintf('Plots_diff_chi_zero/Tumor_diff_chi_zero_%ddays.png',3*(i-1)));%png
    
    figure(10)
    surf(x,y,h_diff0(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    axis tight
    drawnow
    if (i==1)
        title('Initial acidity difference', 'Fontsize', 13);
    else
        title(['Acidity difference at t = ', num2str(3*(i-1)),' days '], 'Fontsize', 13);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    %     saveas(gcf,sprintf('Plots__diff_chi_zero/Acidity__diff_chi_zero_%ddays',3*(i-1)),'epsc');
    
    saveas(gcf,sprintf('Plots_diff_chi_zero/Acidity_diff_chi_zero_%ddays.png',3*(i-1)))
    
    figure(11)
    surf(x,y,m_diff0(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    axis tight
    drawnow
    if (i==1)
        title('Initial normal tissue difference', 'Fontsize', 13);
    else
        title(['Normal tissue difference at t = ', num2str(3*(i-1)),' days '], 'Fontsize', 13);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    %     saveas(gcf,sprintf('Plots__diff_chi_zero/Normal__diff_chi_zero_%ddays',3*(i-1)),'epsc');
    
    saveas(gcf,sprintf('Plots_diff_chi_zero/Normal_diff_chi_zero_%ddays.png',3*(i-1)))
    
    figure(12)
    surf(x,y,n_diff0(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    axis tight
    drawnow
    if (i==1)
        title('Initial necrotic matter difference', 'Fontsize', 13);
    else
        title(['Necrotic matter difference at t = ', num2str(3*(i-1)),' days '], 'Fontsize', 13);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    %     saveas(gcf,sprintf('Plots__diff_chi_zero/Necrotic__diff_chi_zero_%ddays',3*(i-1)),'epsc');
    
    saveas(gcf,sprintf('Plots_diff_chi_zero/Necrotic_diff_chi_zero_%ddays.png',3*(i-1)))
end


%% Videos
% In this section Videos are saved in Videos__diff_chi_double folder
% To check if Videos_diff_chi_double folder exists otherwise to create
if not(isfolder('Videos_diff_chi_double'))
    mkdir('Videos_diff_chi_double')
end

% to check which video profile supports available in the machine
% if mp4 is not supported then avi format will be used
profiles = VideoWriter.getProfiles();
check_mp4_support = find(ismember({profiles.Name},'MPEG-4'),1);

if isempty(check_mp4_support)
    video_ext = '.avi';
    v_pro  = 'Motion JPEG AVI';
else
    video_ext = '.mp4';
    v_pro = 'MPEG-4';
end

% These loops save the solution video  each 3 days
videofile = VideoWriter(strcat('Videos_diff_chi_double/Tumor_diff_chi_double', video_ext),v_pro);
videofile.Quality = 100;
videofile.FrameRate = 10;
open(videofile);
for i = 1: size(c_diff2,3)
    figure(4)
    surf(x,y,c_diff2(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    axis tight
    drawnow
    if (i==1)
        title('Initial glioma cells difference', 'Fontsize', 13);
    else
        title(['Glioma cells difference at t = ', num2str(3*(i-1)), ' days '], 'Fontsize', 13);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    
    F= getframe(gcf);
    writeVideo(videofile,F);
end
close(videofile);

videofile = VideoWriter(strcat('Videos_diff_chi_double/Acidity_diff_chi_double', video_ext),v_pro);
videofile.Quality = 100;
videofile.FrameRate = 10;
open(videofile);
for i = 1: size(h_diff2,3)
    figure(5)
    surf(x,y,h_diff2(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    axis tight
    drawnow
    if (i==1)
        title('Initial acidity difference', 'Fontsize', 13);
    else
        title(['Acidity difference at t = ', num2str(3*(i-1)),' days '], 'Fontsize', 13);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    axis tight
    F= getframe(gcf);
    writeVideo(videofile,F)
    
end
close(videofile);

videofile = VideoWriter(strcat('Videos_diff_chi_double/Normal_diff_chi_double', video_ext),v_pro);
videofile.Quality = 100;
videofile.FrameRate = 10;
open(videofile);
for i = 1: size(m_diff2,3)
    figure(6)
    surf(x,y,m_diff2(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    axis tight
    drawnow
    if (i==1)
        title('Initial normal tissue difference', 'Fontsize', 13);
    else
        title(['Normal tissue difference at t = ', num2str(3*(i-1)),' days '], 'Fontsize', 13);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    
    F= getframe(gcf);
    writeVideo(videofile,F)
    
end
close(videofile);

videofile = VideoWriter(strcat('Videos_diff_chi_double/Necrotic_diff_chi_double', video_ext),v_pro);
videofile.Quality = 100;
videofile.FrameRate = 10;
open(videofile);
for i = 1: size(n_diff2,3)
    figure(8)
    surf(x,y,n_diff2(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    axis tight
    drawnow
    if (i==1)
        title('Initial necrotic matter difference', 'Fontsize', 13);
    else
        title(['Necrotic matter difference at t = ', num2str(3*(i-1)),' days '], 'Fontsize', 13);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    
    F= getframe(gcf);
    writeVideo(videofile,F)
    
end
close(videofile);
%% % This loop saves the results for each 30 days in folder Plots__diff_chi_double in png
% or eps format with days in file names

% to check wheather Plots__diff_chi_double folder exists otherwise it makes a folder Plots__diff_chi_double
if not(isfolder('Plots_diff_chi_double'))
    mkdir('Plots_diff_chi_double')
end

for i = 1:10:size(c_diff2,3)
    figure(9)
    surf(x,y,c_diff2(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    axis tight
    drawnow
    if (i==1)
        title('Initial glioma cells difference', 'Fontsize', 13);
    else
        title(['Glioma cells difference at t = ', num2str(3*(i-1)), ' days '], 'Fontsize', 13);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    
    saveas(gcf,sprintf('Plots_diff_chi_double/Tumor_diff_chi_double_%ddays.png',3*(i-1)));%png
    
    figure(10)
    surf(x,y,h_diff2(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    axis tight
    drawnow
    if (i==1)
        title('Initial acidity difference', 'Fontsize', 13);
    else
        title(['Acidity difference at t = ', num2str(3*(i-1)),' days '], 'Fontsize', 13);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    %     saveas(gcf,sprintf('Plots__diff_chi_double/Acidity__diff_chi_double_%ddays',3*(i-1)),'epsc');
    
    saveas(gcf,sprintf('Plots_diff_chi_double/Acidity_diff_chi_double_%ddays.png',3*(i-1)))
    
    figure(11)
    surf(x,y,m_diff2(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    axis tight
    drawnow
    if (i==1)
        title('Initial normal tissue difference', 'Fontsize', 13);
    else
        title(['Normal tissue difference at t = ', num2str(3*(i-1)),' days '], 'Fontsize', 13);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    %     saveas(gcf,sprintf('Plots__diff_chi_double/Normal__diff_chi_double_%ddays',3*(i-1)),'epsc');
    
    saveas(gcf,sprintf('Plots_diff_chi_double/Normal_diff_chi_double_%ddays.png',3*(i-1)))
    
    figure(12)
    surf(x,y,n_diff2(:,:,i)')
    view(0,90)
    colorbar
    shading interp
    colormap jet
    axis tight
    drawnow
    if (i==1)
        title('Initial necrotic matter difference', 'Fontsize', 13);
    else
        title(['Necrotic matter difference at t = ', num2str(3*(i-1)),' days '], 'Fontsize', 13);
    end
    xlabel('X' , 'Fontsize', 15);
    ylabel('Y' , 'Fontsize', 15);
    %     saveas(gcf,sprintf('Plots__diff_chi_double/Necrotic__diff_chi_double_%ddays',3*(i-1)),'epsc');
    
    saveas(gcf,sprintf('Plots_diff_chi_double/Necrotic_diff_chi_double_%ddays.png',3*(i-1)))
end

%% To save the workspace
if not(isfolder('Mat_files'))
    mkdir('Mat_files')
end
save('Mat_files/main_2D_diff_chi.mat');
%%
toc