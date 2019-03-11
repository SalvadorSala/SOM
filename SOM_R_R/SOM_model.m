% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
% %%%Description:
% A mathematical model of the development of correctly oriented topographic
% maps in visual cortex is presented that takes into account the 
% synchronization between the two retinas mediated by retino-retinal 
% projections and the effect of the molecular guidance cues.
% 
% We used a simplified model of self-organizing map (SOM) as described by
% Kohonen [1] and a modified version of the script from the book [2]
% For more detailed information about the model refer to the paper[3].

% %%%%%%%%%%%%      References
% [1]T. Kohonen, \Self-organized formation of topologically correct feature
% maps. biological cybernetics 43(1), 59-69," vol. 43, pp. 59{69, 01 1982.
% [2]Trappenberg, Thomas. (2002). Fundamentals of Computational Neuroscience. 
% [3]


% %%%%%%%%%%%%%    Parameters
% 
% lambda_aux.........................................Weight decay term 0.1
% Tau......................................Time constant of ? decay 0.00018
% sig_SOM...............................Lateral interaction influence decay 
% nn................................................Number of cells 11 × 11
% sigma noise..........................Standard deviation (SD) of position noise 
% sigma molecular.........................SD of the molecular gradient [0.1, 10]
% 
%
% To study the function of the retino-retinal protections as an activity
% dependent mechanism that contributes to the formation of correct topographic
% maps, we introduced a parameter that varied in 4 ways how the presynaptic
% layer stimulated the postsynaptic layer,

%%%%%%%%%%%%%%%    Stimulus type
%  stim={'ret_wave','sincro_ini','sincro_azar','azar'};
% Ret_wave....................................................Retinal Wave
% Sincro_ini..............................................Initial sinchrony
% sincro_azar...................................Equal Random in both sheets
% azar...................................................Random Simtulation

%%%%%%%%%%%%%     Initial Weights

% weights='ephrin_gradient'; %'azar', ephrin_gradient
% type 'azar' for initial random weights
% type 'ephrin_gradient' for initial molecular gradient

% By modelling the development of the right and left postsynaptic sheet
% simultaneously we are able to study how these parameters affect the correct
% outcome of bilateral symmetry. The model returns correct results when the
% orientation of both postsynaptic sheets is the same as the orientation in the
% presynaptic sheet. On the other hand,Incorrect results indicate different
% orientations between pre and postsynaptic sheets or incorrect unfolding, which
% produce disruptions on the topographic map

%%%%%%%%%%%      Steps to follow.

% 1) Run the script
% 2) Time 5.53h with Intel Core i7 3.70GHz
% 3) Run figure1.m and figure2.m to obtain the figures from[2].

%%%%%%%%%                                         16/10/18   %%%%%%%%%%%%%%
%%%%%%Authors:                              AJV and SSP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
tic

addpath([cd, '/functions_SOM']);
load('correct_ori.mat')

%%%%%%%%%%%%        Parameters                             
%%%%%%%%%%%%       

iter=100; %number of iterations
counter=0;
sponta=100; %number of similiar initial localiced random inputs
Tau=0.00018; %decaying rate of lambda
tot_trial=12000; %total number of trials 

%%%%%%%%%%%%%%%    Stimulus type
%%%%%%%%%%%%%%%
stim={'ret_wave','sincro_ini','sincro_azar','azar'};

%%%%%%%%%%%%%%%    Initial Weights
%%%%%%%%%%%%%%%
weights='ephrin_gradient'; %'azar', ephrin_gradient

model=cell(1,6);model{1}='lambda';model{1,2}='coor ret';model{1,3}='alfa'; model{1,4}='sig_ini';
model{1,5}='sig'; model{1,6}='Good ori';  model{1,7}='quality';
calidad=[];

for alfan=0:0.03:0.21
for sig_ephrin=[0.1,1:10]
    for input=1:length(stim)
for z=1:iter
     
clearvars -except model calidad alfan z  stim tot_trial estimulo Tau nstim sponta ret_wave1 ret_wave2 z cum_stim1...
    cum_stim2 R_c1_cor L_c1_cor input R_c2_cor L_c2_cor input ini_stim counter iter weights alfan sig_ephrin  ;...

estimulo=stim(input);
if strcmp('ret_wave',estimulo)
      
angle=0:0.001:2*pi; %angle of retinal wave
r_wave=0:0.01:0.7;  % radius retinal wave
[ret_wave1]=ret_wave_simple5(r_wave,50,angle,0.4,0.4);
[ret_wave2]=ret_wave_simple5(r_wave,50,angle,0.6,0.4);
nstim1=size(ret_wave1,2); %n stimuly that form the retinal wave
nstim2=size(ret_wave2,2); %n stimulytant form the retinal wave
elseif strcmp('azar',estimulo)
nstim1=0;nstim2=0;
elseif strcmp('sincro_azar',estimulo)
nstim1=0;nstim2=0;
elseif strcmp('sincro_ini',estimulo)
ax = 0.45;ay=0.45;
bx = 0.35;by=0.35;
L_ini_stim =[(bx-ax).*rand(1000,1) + ax,(by-ay).*rand(1000,1) + ay ]';
ax = 0.65; ay=0.45;
bx = 0.55; by=0.35;
R_ini_stim = [(bx-ax).*rand(1000,1) + ax,(by-ay).*rand(1000,1) + ay ]';
nstim1=0;nstim2=0;
end


nn=11; lambda_aux=0.1; ntrial=0; alfa=1;  sig1=1/(2*sig_ephrin^2); grid=[1:nn]./10-0.1;
sig_SOM=2; sig2=1/(2*sig_SOM^2);

[Rxmat_pos,Rymat_pos]=meshgrid(grid,grid); Rxmat_pos_nois=Rxmat_pos+(alfa*(normrnd(0,alfan,nn,nn))); Rymat_pos_nois=Rymat_pos+(alfa*(normrnd(0,alfan,nn,nn)));
[Rxmat_pre,Rymat_pre]=meshgrid(grid,grid); Rxmat_pre_nois=Rxmat_pre+(alfa*(normrnd(0,alfan,nn,nn))); Rymat_pre_nois=Rymat_pre+(alfa*(normrnd(0,alfan,nn,nn)));

[Lxmat_pos,Lymat_pos]=meshgrid(grid,grid); Lxmat_pos_nois=Lxmat_pos+(alfa*(normrnd(0,alfan,nn,nn))); Lymat_pos_nois=Lymat_pos+(alfa*(normrnd(0,alfan,nn,nn)));
[Lxmat_pre,Lymat_pre]=meshgrid(grid,grid); Lxmat_pre_nois=Lxmat_pre+(alfa*(normrnd(0,alfan,nn,nn))); Lymat_pre_nois=Lymat_pre+(alfa*(normrnd(0,alfan,nn,nn)));


if strcmp('ephrin_gradient',weights)
pesos=zeros(nn,nn,nn^2);

for i=1:nn^2
Rx=Rxmat_pos_nois(i);Lx=Lxmat_pos_nois(i); Ry=Rymat_pos_nois(i);Ly=Lymat_pos_nois(i);

gauss = @(x,y) exp(-((Rxmat_pre_nois-x).^2+(Rymat_pre_nois-y).^2)*sig1);
R_pesos(:,:,i)=gauss(Rx,Ry);

gauss = @(x,y) exp(-((Lxmat_pre_nois-x).^2+(Lymat_pre_nois-y).^2)*sig1);
L_pesos(:,:,i)=gauss(Lx,Ly);

end
else
R_pesos=rand(nn,nn,nn^2);
L_pesos=rand(nn,nn,nn^2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Calculate initial weight of each postsynaptic neuron, center of mass%%%%

[Rpesoi]=peso_ini(Rxmat_pre,Rymat_pre,R_pesos);
[Lpesoi]=peso_ini(Lxmat_pre,Lymat_pre,L_pesos);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Initial centres of prefered features %%%%%%%%%%%%%%%%%%%%%%%%

R_c1=reshape(Rpesoi(:,1),[nn,nn]);
R_c2=reshape(Rpesoi(:,2),[nn,nn]);

L_c1=reshape(Lpesoi(:,1),[nn,nn]);
L_c2=reshape(Lpesoi(:,2),[nn,nn]);


[X,Y]=meshgrid(1:nn,1:nn); ntrial=0; 

%%%%%%%%% training session  
%%%%%%%%% 
   while ntrial<tot_trial

lambda=exp(-Tau*ntrial)*lambda_aux;

if ntrial<nstim1
    r_in=ret_wave1(:,ntrial+1);
    
elseif strcmp('sincro_ini',estimulo) && ntrial<sponta;
    
    r_in=L_ini_stim(:,ntrial+1);
else
     [r_in]=[rand;rand];
end
       r=exp(-(L_c1-r_in(1)).^2-(L_c2-r_in(2)).^2); 
      [rmax,x_winner]=max(max(r)); [rmax,y_winner]=max(max(r'));
      r=exp(-((X-x_winner).^2+(Y-y_winner).^2)*sig2); 

      L_c1=L_c1+lambda*r.*(r_in(1)-L_c1); 
      L_c2=L_c2+lambda*r.*(r_in(2)-L_c2); 
      
    r_in_L=r_in;
if ntrial<nstim2;
          r_in=ret_wave2(:,ntrial+1);

elseif strcmp('sincro_azar',estimulo)
elseif strcmp('sincro_ini',estimulo) && ntrial<sponta;
    
    r_in=R_ini_stim(:,ntrial+1);
    
else
     [r_in]=[rand;rand];
end

      r=exp(-(R_c1-r_in(1)).^2-(R_c2-r_in(2)).^2); 
     [rmax,x_winner]=max(max(r)); [rmax,y_winner]=max(max(r'));
     r=exp(-((X-x_winner).^2+(Y-y_winner).^2)*sig2);
  
      R_c1=R_c1+lambda*r.*(r_in(1)-R_c1); 
      R_c2=R_c2+lambda*r.*(r_in(2)-R_c2); 
      
    r_in_R=r_in;

      
       ntrial=ntrial+1; 
   end
   
 
%%%%%%%%%%%%% Quantify map orientation
%%%%%%%%%%%%%

if    (corr2(R_c1,L_c1)>0.95 &&  corr2(R_c1_cor,L_c1)>0.95)   &&    (corr2(R_c2,L_c2)>0.95 && corr2(R_c2_cor,L_c2)>0.95)
    
counter=counter+1;

end

Lq = quality(nn, nn^2, Y,X,L_c1,L_c2);
Rq = quality(nn, nn^2, Y,X,R_c1,R_c2);
[L_plega,Lcruce,Lgirado]=correct_plega(L_c1,L_c2,nn,0.7);
[R_plega,Rcruce,Rgirado]=correct_plega(R_c1,R_c2,nn,0.7);


calidad=[calidad ;Lq,L_plega,Lcruce,Lgirado,Rq,R_plega,Rcruce,Rgirado,(corr2(R_c1,L_c1)>0.95 &&  corr2(R_c1_cor,L_c1)>0.95)   &&    (corr2(R_c2,L_c2)>0.95 && corr2(R_c2_cor,L_c2)>0.95)];
toc
end

toc
good_ori=(counter/iter)*100
point=size(model,1)+1;
model{point,1}=lambda;
model{point,2}=estimulo;
model{point,3}=alfan;
model{point,4}=sig_ephrin;
model{point,5}=sig_SOM;
model{point,6}=good_ori;
model{point,7}=calidad;

good_ori=0;counter=0;calidad=[];


end
end
end
save model
