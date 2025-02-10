% This is the main script, which calls all the related function and
% computes the MP model (4.4.39 of the thesis) in 2-D(non-dimensional form) 
% in three cases:
%      1. Base case (compute_three_occ()): Experiment 4.1 of the thesis (model 4.4.39)
%      2. Kij’s equal to K_cm (compute_all_K_equal_I()): Experiment 4.4 of the thesis (model 4.4.39)
%      3. Kij’s equal to K_mn (compute_all_K_equal_II()): Experiment 4.4 of the thesis (model 4.4.39)
%      4. Kij’s equal to K_cn (compute_all_K_equal_III()): Experiment 4.4 of the thesis (model 4.4.39)

% and saves the results in the folder Plots_diff_K_I, Plots_diff_K_II, Plots_diff_K_III,
% Videos_diff_K_I, Videos_diff_K_II and Videos_diff_K_III.
% The corresponding reults are included in Figures 4.7, 4.8 and 4.9 of the thesis

clear all;
clc;
close all;
tic
[x,y,c,m,n,h] = compute_three_occ();
[x1,y1,c1,m1,n1,h1] = compute_all_K_equal_I();
[x2,y2,c2,m2,n2,h2] = compute_all_K_equal_II();
[x3,y3,c3,m3,n3,h3] = compute_all_K_equal_III();

c_diff1 = c - c1;
m_diff1 = m - m1;
n_diff1 = n - n1;
h_diff1 = h - h1;

c_diff2 = c - c2;
m_diff2 = m - m2;
n_diff2 = n - n2;
h_diff2 = h - h2;

c_diff3 = c - c3;
m_diff3 = m - m3;
n_diff3 = n - n3;
h_diff3 = h - h3;


%%
toc
fprintf('Computation Done! Saving the plots and videos...\n')
tic
%% Videos
% In this section Videos are saved in Videos__diff_K_low folder
% To check if Videos_diff_K_low folder exists otherwise to create
if not(isfolder('Videos_diff_K_equal_I'))
    mkdir('Videos_diff_K_equal_I')
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
videofile = VideoWriter(strcat('Videos_diff_K_equal_I/Tumor_diff_K_equal_I', video_ext),v_pro);
videofile.Quality = 100;
videofile.FrameRate = 10;
open(videofile);
for i = 1: size(c_diff1,3)
    figure(4)
    surf(x,y,c_diff1(:,:,i)')
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

videofile = VideoWriter(strcat('Videos_diff_K_equal_I/Acidity_diff_K_equal_I', video_ext),v_pro);
videofile.Quality = 100;
videofile.FrameRate = 10;
open(videofile);
for i = 1: size(h_diff1,3)
    figure(5)
    surf(x,y,h_diff1(:,:,i)')
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

videofile = VideoWriter(strcat('Videos_diff_K_equal_I/Normal_diff_K_equal_I', video_ext),v_pro);
videofile.Quality = 100;
videofile.FrameRate = 10;
open(videofile);
for i = 1: size(m_diff1,3)
    figure(6)
    surf(x,y,m_diff1(:,:,i)')
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

videofile = VideoWriter(strcat('Videos_diff_K_equal_I/Necrotic_diff_K_equal_I', video_ext),v_pro);
videofile.Quality = 100;
videofile.FrameRate = 10;
open(videofile);
for i = 1: size(n_diff1,3)
    figure(8)
    surf(x,y,n_diff1(:,:,i)')
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
%% % This loop saves the results for each 30 days in folder Plots__diff_K_low in png
% or eps format with days in file names

% to check wheather Plots__diff_K_low folder exists otherwise it makes a folder Plots__diff_K_low
if not(isfolder('Plots_diff_K_equal_I'))
    mkdir('Plots_diff_K_equal_I')
end

for i = 1:10:size(c_diff1,3)
    figure(9)
    surf(x,y,c_diff1(:,:,i)')
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
    
    saveas(gcf,sprintf('Plots_diff_K_equal_I/Tumor_diff_K_equal_I_%ddays.png',3*(i-1)));%png
    
    figure(10)
    surf(x,y,h_diff1(:,:,i)')
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
    %     saveas(gcf,sprintf('Plots__diff_K_low/Acidity__diff_K_low_%ddays',3*(i-1)),'epsc');
    
    saveas(gcf,sprintf('Plots_diff_K_equal_I/Acidity_diff_K_equal_I_%ddays.png',3*(i-1)))
    
    figure(11)
    surf(x,y,m_diff1(:,:,i)')
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
    %     saveas(gcf,sprintf('Plots__diff_K_low/Normal__diff_K_low_%ddays',3*(i-1)),'epsc');
    
    saveas(gcf,sprintf('Plots_diff_K_equal_I/Normal_diff_K_equal_I_%ddays.png',3*(i-1)))
    
    figure(12)
    surf(x,y,n_diff1(:,:,i)')
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
    %     saveas(gcf,sprintf('Plots__diff_K_low/Necrotic__diff_K_low_%ddays',3*(i-1)),'epsc');
    
    saveas(gcf,sprintf('Plots_diff_K_equal_I/Necrotic_diff_K_equal_I_%ddays.png',3*(i-1)))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Videos
% In this section Videos are saved in Videos__diff_K_low folder
% To check if Videos_diff_K_low folder exists otherwise to create
if not(isfolder('Videos_diff_K_equal_II'))
    mkdir('Videos_diff_K_equal_II')
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
videofile = VideoWriter(strcat('Videos_diff_K_equal_II/Tumor_diff_K_equal_II', video_ext),v_pro);
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

videofile = VideoWriter(strcat('Videos_diff_K_equal_II/Acidity_diff_K_equal_II', video_ext),v_pro);
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

videofile = VideoWriter(strcat('Videos_diff_K_equal_II/Normal_diff_K_equal_II', video_ext),v_pro);
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

videofile = VideoWriter(strcat('Videos_diff_K_equal_II/Necrotic_diff_K_equal_II', video_ext),v_pro);
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
%% % This loop saves the results for each 30 days in folder Plots__diff_K_low in png
% or eps format with days in file names

% to check wheather Plots__diff_K_low folder exists otherwise it makes a folder Plots__diff_K_low
if not(isfolder('Plots_diff_K_equal_II'))
    mkdir('Plots_diff_K_equal_II')
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
    
    saveas(gcf,sprintf('Plots_diff_K_equal_II/Tumor_diff_K_equal_II_%ddays.png',3*(i-1)));%png
    
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
    %     saveas(gcf,sprintf('Plots__diff_K_low/Acidity__diff_K_low_%ddays',3*(i-1)),'epsc');
    
    saveas(gcf,sprintf('Plots_diff_K_equal_II/Acidity_diff_K_equal_II_%ddays.png',3*(i-1)))
    
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
    %     saveas(gcf,sprintf('Plots__diff_K_low/Normal__diff_K_low_%ddays',3*(i-1)),'epsc');
    
    saveas(gcf,sprintf('Plots_diff_K_equal_II/Normal_diff_K_equal_II_%ddays.png',3*(i-1)))
    
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
    %     saveas(gcf,sprintf('Plots__diff_K_low/Necrotic__diff_K_low_%ddays',3*(i-1)),'epsc');
    
    saveas(gcf,sprintf('Plots_diff_K_equal_II/Necrotic_diff_K_equal_II_%ddays.png',3*(i-1)))
end
%----------------------------------------------------------------------------------------------------------------
%% Videos
% In this section Videos are saved in Videos__diff_K_low folder
% To check if Videos_diff_K_low folder exists otherwise to create
if not(isfolder('Videos_diff_K_equal_III'))
    mkdir('Videos_diff_K_equal_III')
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
videofile = VideoWriter(strcat('Videos_diff_K_equal_III/Tumor_diff_K_equal_III', video_ext),v_pro);
videofile.Quality = 100;
videofile.FrameRate = 10;
open(videofile);
for i = 1: size(c_diff3,3)
    figure(4)
    surf(x,y,c_diff3(:,:,i)')
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

videofile = VideoWriter(strcat('Videos_diff_K_equal_III/Acidity_diff_K_equal_III', video_ext),v_pro);
videofile.Quality = 100;
videofile.FrameRate = 10;
open(videofile);
for i = 1: size(h_diff3,3)
    figure(5)
    surf(x,y,h_diff3(:,:,i)')
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

videofile = VideoWriter(strcat('Videos_diff_K_equal_III/Normal_diff_K_equal_III', video_ext),v_pro);
videofile.Quality = 100;
videofile.FrameRate = 10;
open(videofile);
for i = 1: size(m_diff3,3)
    figure(6)
    surf(x,y,m_diff3(:,:,i)')
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

videofile = VideoWriter(strcat('Videos_diff_K_equal_III/Necrotic_diff_K_equal_III', video_ext),v_pro);
videofile.Quality = 100;
videofile.FrameRate = 10;
open(videofile);
for i = 1: size(n_diff3,3)
    figure(8)
    surf(x,y,n_diff3(:,:,i)')
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
%% % This loop saves the results for each 30 days in folder Plots__diff_K_low in png
% or eps format with days in file names

% to check wheather Plots__diff_K_low folder exists otherwise it makes a folder Plots__diff_K_low
if not(isfolder('Plots_diff_K_equal_III'))
    mkdir('Plots_diff_K_equal_III')
end

for i = 1:10:size(c_diff3,3)
    figure(9)
    surf(x,y,c_diff3(:,:,i)')
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
    
    saveas(gcf,sprintf('Plots_diff_K_equal_III/Tumor_diff_K_equal_III_%ddays.png',3*(i-1)));%png
    
    figure(10)
    surf(x,y,h_diff3(:,:,i)')
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
    %     saveas(gcf,sprintf('Plots__diff_K_low/Acidity__diff_K_low_%ddays',3*(i-1)),'epsc');
    
    saveas(gcf,sprintf('Plots_diff_K_equal_III/Acidity_diff_K_equal_III_%ddays.png',3*(i-1)))
    
    figure(11)
    surf(x,y,m_diff3(:,:,i)')
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
    %     saveas(gcf,sprintf('Plots__diff_K_low/Normal__diff_K_low_%ddays',3*(i-1)),'epsc');
    
    saveas(gcf,sprintf('Plots_diff_K_equal_III/Normal_diff_K_equal_III_%ddays.png',3*(i-1)))
    
    figure(12)
    surf(x,y,n_diff3(:,:,i)')
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
    %     saveas(gcf,sprintf('Plots__diff_K_low/Necrotic__diff_K_low_%ddays',3*(i-1)),'epsc');
    
    saveas(gcf,sprintf('Plots_diff_K_equal_III/Necrotic_diff_K_equal_III_%ddays.png',3*(i-1)))
end

%% To save the workspace
if not(isfolder('Mat_files'))
    mkdir('Mat_files')
end
save('Mat_files/main_2D_K_equal.mat');
%%
toc