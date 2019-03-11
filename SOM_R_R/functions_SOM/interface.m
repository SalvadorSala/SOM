function varargout = interface(varargin)
% INTERFACE MATLAB code for interface.fig
%      INTERFACE, by itself, creates a new INTERFACE or raises the existing
%      singleton*.
%
%      H = INTERFACE returns the handle to a new INTERFACE or the handle to
%      the existing singleton*.
%
%      INTERFACE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INTERFACE.M with the given input arguments.
%
%      INTERFACE('Property','Value',...) creates a new INTERFACE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before interface_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to interface_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help interface

% Last Modified by GUIDE v2.5 11-Jun-2018 12:03:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @interface_OpeningFcn, ...
                   'gui_OutputFcn',  @interface_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before interface is made visible.
function interface_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to interface (see VARARGIN)

% Choose default command line output for interface
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes interface wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = interface_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%Initial conections%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
load('correct_ori.mat')

iter=100; %number of iterations
counter=0;
sponta=100; %number of similiar initial localiced random inputs
Tau=0.00018; %decaying rate of lambda
tot_trial=12000; %total number of trials

stim={'ret_wave','sincro_ini','sincro_azar','azar'}; %'azar','sincro_azar','sincro_ini',ret_wave
weights='ephrin_gradient'; %'azar', ephrin_gradient

model=cell(1,6);model{1}='lambda';model{1,2}='coor ret';model{1,3}='alfa'; model{1,4}='sig_ini';
model{1,5}='sig'; model{1,6}='Good ori';  model{1,7}='quality'
calidad=[];

for alfan=0:0.03:0.21
for sig_ephrin=2%[0.1,1:10]
    for input=1:length(stim)
for z=1:iter
     
clearvars -except model calidad alfan z  stim tot_trial estimulo Tau nstim sponta ret_wave1 ret_wave2 z cum_stim1...
    cum_stim2 R_c1_cor L_c1_cor input R_c2_cor L_c2_cor input ini_stim counter iter weights alfan sig_ephrin ;...

estimulo=stim(input);
if strcmp('ret_wave',estimulo)
      
angle=0:0.001:2*pi; %angle of retinal wave
r_wave=0:0.01:0.3;  % radius od retinal wave
[ret_wave1]=ret_wave_simple2(r_wave,50,angle,0.5,0.4);
[ret_wave2]=ret_wave_simple2(r_wave,50,angle,0.5,0.4);
nstim=size(ret_wave1,2); %n stimulos que forman la onda el resto es espontaneo
elseif strcmp('azar',estimulo)
nstim=0;
elseif strcmp('sincro_azar',estimulo)
nstim=0;
elseif strcmp('sincro_ini',estimulo)
a = 0.4;
b = 0.5;
L_ini_stim = [(b-a).*rand(1000,1) + a,(b-a).*rand(1000,1) + a ]';
R_ini_stim = [(b-a).*rand(1000,1) + a,(b-a).*rand(1000,1) + a ]';
nstim=0;
end


nn=11; lambda_aux=0.1; ntrial=0; alfa=1;  sig1=1/(2*sig_ephrin^2); grid=[1:nn]./10-0.1;
sig_SOM=2; sig2=1/(2*sig_SOM^2); %alfan=0.12 sig_ephrin=1

[Rxmat_pos,Rymat_pos]=meshgrid(grid,grid); Rxmat_pos_nois=Rxmat_pos+(alfa*(normrnd(0,alfan,nn,nn))); Rymat_pos_nois=Rymat_pos+(alfa*(normrnd(0,alfan,nn,nn)));
[Rxmat_pre,Rymat_pre]=meshgrid(grid,grid); Rxmat_pre_nois=Rxmat_pre+(alfa*(normrnd(0,alfan,nn,nn))); Rymat_pre_nois=Rymat_pre+(alfa*(normrnd(0,alfan,nn,nn)));

[Lxmat_pos,Lymat_pos]=meshgrid(grid,grid); Lxmat_pos_nois=Lxmat_pos+(alfa*(normrnd(0,alfan,nn,nn))); Lymat_pos_nois=Lymat_pos+(alfa*(normrnd(0,alfan,nn,nn)));
[Lxmat_pre,Lymat_pre]=meshgrid(grid,grid); Lxmat_pre_nois=Lxmat_pre+(alfa*(normrnd(0,alfan,nn,nn))); Lymat_pre_nois=Lymat_pre+(alfa*(normrnd(0,alfan,nn,nn)));


if strcmp('ephrin_gradient',weights)
pesos=zeros(nn,nn,nn^2);

for i=1:nn^2
Rx=Rxmat_pos_nois(i);Lx=Lxmat_pos_nois(i);        Ry=Rymat_pos_nois(i);Ly=Lymat_pos_nois(i);

gauss = @(x,y) exp(-((Rxmat_pre_nois-x).^2+(Rymat_pre_nois-y).^2)*sig1);
R_pesos(:,:,i)=gauss(Rx,Ry);

gauss = @(x,y) exp(-((Lxmat_pre_nois-x).^2+(Lymat_pre_nois-y).^2)*sig1);
L_pesos(:,:,i)=gauss(Lx,Ly);

end
else
R_pesos=rand(nn,nn,nn^2);
L_pesos=rand(nn,nn,nn^2);
end


%%%Calcular peso incial de cada neurona postsinaptica, centro de masas
[Rpesoi]=peso_ini(Rxmat_pre,Rymat_pre,R_pesos);
[Lpesoi]=peso_ini(Lxmat_pre,Lymat_pre,L_pesos);

% Initial centres of prefered features:

R_c1=reshape(Rpesoi(:,1),[nn,nn]);
R_c2=reshape(Rpesoi(:,2),[nn,nn]);

L_c1=reshape(Lpesoi(:,1),[nn,nn]);
L_c2=reshape(Lpesoi(:,2),[nn,nn]);


[X,Y]=meshgrid(1:nn,1:nn); ntrial=0; 


   %% training session 
   while ntrial<tot_trial
                   if(mod(ntrial,100)==0)
% 
       clf; 
   hold on;
   axis square;
   axis([0 1 0 1]);
   subplot(1,2,1)
    plot(L_c1,L_c2,'k'); hold on;plot(L_c1',L_c2','k');axis([0 1 0 1]); axis square
       hold on;
   plot(L_c1(1,1),L_c2(1,1),'Marker','x','color' , 'k');  axis([0 1 0 1]); axis square
    plot(L_c1(11,11),L_c2(11,11),'Marker','x','color' , 'b');  axis([0 1 0 1]); axis square
     plot(L_c1(11,1),L_c2(11,1),'Marker','x','color' , 'r');  axis([0 1 0 1]); axis square
      plot(L_c1(1,11),L_c2(1,11),'Marker','x','color' , 'y');  axis([0 1 0 1]); axis square
      tstring=[int2str(ntrial) ' examples']; title(tstring);
    subplot(1,2,2)
     plot(R_c1,R_c2,'k');hold on; plot(R_c1',R_c2','k'); hold on; axis([0 1 0 1]); axis square
     hold on
   plot(R_c1(1,1),R_c2(1,1),'Marker','x','color' , 'k');  axis([0 1 0 1]); axis square
    plot(R_c1(11,11),R_c2(11,11),'Marker','x','color' , 'b');  axis([0 1 0 1]); axis square
     plot(R_c1(11,1),R_c2(11,1),'Marker','x','color' , 'r');  axis([0 1 0 1]); axis square
      plot(R_c1(1,11),R_c2(1,11),'Marker','x','color' , 'y');  axis([0 1 0 1]); axis square
 movegui('north')
   tstring=[int2str(ntrial) ' examples']; title(tstring);
%             end

lambda=exp(-Tau*ntrial)*lambda_aux;


if ntrial<nstim
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
      
    
if ntrial<nstim;
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

       ntrial=ntrial+1; 
   end
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%Cuantificar la oreintacion correcta%%%%%%%%%%%%%%%%%%%%%%%%%


if    (corr2(R_c1,L_c1)>0.95 &&  corr2(R_c1_cor,L_c1)>0.95)   &&    (corr2(R_c2,L_c2)>0.95 && corr2(R_c2_cor,L_c2)>0.95)
    
counter=counter+1;

end

Lq = quality(nn, nn^2, Y,X,L_c1,L_c2);
Rq = quality(nn, nn^2, Y,X,R_c1,R_c2);

calidad=[calidad ;Lq,Rq, (corr2(R_c1,L_c1)>0.95 &&  corr2(R_c1_cor,L_c1)>0.95)   &&    (corr2(R_c2,L_c2)>0.95 && corr2(R_c2_cor,L_c2)>0.95)];
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

% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called




% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');   
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)

tic
load('correct_ori.mat')

iter=100; %number of iterations
counter=0;
sponta=100; %number of similiar initial localiced random inputs
Tau=0.00018; %decaying rate of lambda
tot_trial=12000; %total number of trials

stim={'ret_wave','sincro_ini','sincro_azar','azar'}; %'azar','sincro_azar','sincro_ini',ret_wave
weights='ephrin_gradient'; %'azar', ephrin_gradient

model=cell(1,6);model{1}='lambda';model{1,2}='coor ret';model{1,3}='alfa'; model{1,4}='sig_ini';
model{1,5}='sig'; model{1,6}='Good ori';  model{1,7}='quality'
calidad=[];

for alfan=0:0.03:0.21
for sig_ephrin=2%[0.1,1:10]
    for input=1:length(stim)
for z=1:iter
     
clearvars -except model calidad alfan z  stim tot_trial estimulo Tau nstim sponta ret_wave1 ret_wave2 z cum_stim1...
    cum_stim2 R_c1_cor L_c1_cor input R_c2_cor L_c2_cor input ini_stim counter iter weights alfan sig_ephrin ;...

estimulo=stim(input);
if strcmp('ret_wave',estimulo)
      
angle=0:0.001:2*pi; %angle of retinal wave
r_wave=0:0.01:0.3;  % radius od retinal wave
[ret_wave1]=ret_wave_simple2(r_wave,50,angle,0.5,0.4);
[ret_wave2]=ret_wave_simple2(r_wave,50,angle,0.5,0.4);
nstim=size(ret_wave1,2); %n stimulos que forman la onda el resto es espontaneo
elseif strcmp('azar',estimulo)
nstim=0;
elseif strcmp('sincro_azar',estimulo)
nstim=0;
elseif strcmp('sincro_ini',estimulo)
a = 0.4;
b = 0.5;
L_ini_stim = [(b-a).*rand(1000,1) + a,(b-a).*rand(1000,1) + a ]';
R_ini_stim = [(b-a).*rand(1000,1) + a,(b-a).*rand(1000,1) + a ]';
nstim=0;
end


nn=11; lambda_aux=0.1; ntrial=0; alfa=1;  sig1=1/(2*sig_ephrin^2); grid=[1:nn]./10-0.1;
sig_SOM=2; sig2=1/(2*sig_SOM^2); %alfan=0.12 sig_ephrin=1

[Rxmat_pos,Rymat_pos]=meshgrid(grid,grid); Rxmat_pos_nois=Rxmat_pos+(alfa*(normrnd(0,alfan,nn,nn))); Rymat_pos_nois=Rymat_pos+(alfa*(normrnd(0,alfan,nn,nn)));
[Rxmat_pre,Rymat_pre]=meshgrid(grid,grid); Rxmat_pre_nois=Rxmat_pre+(alfa*(normrnd(0,alfan,nn,nn))); Rymat_pre_nois=Rymat_pre+(alfa*(normrnd(0,alfan,nn,nn)));

[Lxmat_pos,Lymat_pos]=meshgrid(grid,grid); Lxmat_pos_nois=Lxmat_pos+(alfa*(normrnd(0,alfan,nn,nn))); Lymat_pos_nois=Lymat_pos+(alfa*(normrnd(0,alfan,nn,nn)));
[Lxmat_pre,Lymat_pre]=meshgrid(grid,grid); Lxmat_pre_nois=Lxmat_pre+(alfa*(normrnd(0,alfan,nn,nn))); Lymat_pre_nois=Lymat_pre+(alfa*(normrnd(0,alfan,nn,nn)));


if strcmp('ephrin_gradient',weights)
pesos=zeros(nn,nn,nn^2);

for i=1:nn^2
Rx=Rxmat_pos_nois(i);Lx=Lxmat_pos_nois(i);        Ry=Rymat_pos_nois(i);Ly=Lymat_pos_nois(i);

gauss = @(x,y) exp(-((Rxmat_pre_nois-x).^2+(Rymat_pre_nois-y).^2)*sig1);
R_pesos(:,:,i)=gauss(Rx,Ry);

gauss = @(x,y) exp(-((Lxmat_pre_nois-x).^2+(Lymat_pre_nois-y).^2)*sig1);
L_pesos(:,:,i)=gauss(Lx,Ly);

end
else
R_pesos=rand(nn,nn,nn^2);
L_pesos=rand(nn,nn,nn^2);
end


%%%Calcular peso incial de cada neurona postsinaptica, centro de masas
[Rpesoi]=peso_ini(Rxmat_pre,Rymat_pre,R_pesos);
[Lpesoi]=peso_ini(Lxmat_pre,Lymat_pre,L_pesos);

% Initial centres of prefered features:

R_c1=reshape(Rpesoi(:,1),[nn,nn]);
R_c2=reshape(Rpesoi(:,2),[nn,nn]);

L_c1=reshape(Lpesoi(:,1),[nn,nn]);
L_c2=reshape(Lpesoi(:,2),[nn,nn]);


[X,Y]=meshgrid(1:nn,1:nn); ntrial=0; 


   %% training session 
   while ntrial<tot_trial
                   if(mod(ntrial,100)==0)
% 
       clf; 
   hold on;
   axis square;
   axis([0 1 0 1]);
   subplot(1,2,1)
    plot(L_c1,L_c2,'k'); hold on;plot(L_c1',L_c2','k');axis([0 1 0 1]); axis square
       hold on;
   plot(L_c1(1,1),L_c2(1,1),'Marker','x','color' , 'k');  axis([0 1 0 1]); axis square
    plot(L_c1(11,11),L_c2(11,11),'Marker','x','color' , 'b');  axis([0 1 0 1]); axis square
     plot(L_c1(11,1),L_c2(11,1),'Marker','x','color' , 'r');  axis([0 1 0 1]); axis square
      plot(L_c1(1,11),L_c2(1,11),'Marker','x','color' , 'y');  axis([0 1 0 1]); axis square
      tstring=[int2str(ntrial) ' examples']; title(tstring);
    subplot(1,2,2)
     plot(R_c1,R_c2,'k');hold on; plot(R_c1',R_c2','k'); hold on; axis([0 1 0 1]); axis square
     hold on
   plot(R_c1(1,1),R_c2(1,1),'Marker','x','color' , 'k');  axis([0 1 0 1]); axis square
    plot(R_c1(11,11),R_c2(11,11),'Marker','x','color' , 'b');  axis([0 1 0 1]); axis square
     plot(R_c1(11,1),R_c2(11,1),'Marker','x','color' , 'r');  axis([0 1 0 1]); axis square
      plot(R_c1(1,11),R_c2(1,11),'Marker','x','color' , 'y');  axis([0 1 0 1]); axis square
 movegui('north')
   tstring=[int2str(ntrial) ' examples']; title(tstring);
%             end

lambda=exp(-Tau*ntrial)*lambda_aux;


if ntrial<nstim
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
      
    
if ntrial<nstim;
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

       ntrial=ntrial+1; 
   end
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%Cuantificar la oreintacion correcta%%%%%%%%%%%%%%%%%%%%%%%%%


if    (corr2(R_c1,L_c1)>0.95 &&  corr2(R_c1_cor,L_c1)>0.95)   &&    (corr2(R_c2,L_c2)>0.95 && corr2(R_c2_cor,L_c2)>0.95)
    
counter=counter+1;

end

Lq = quality(nn, nn^2, Y,X,L_c1,L_c2);
Rq = quality(nn, nn^2, Y,X,R_c1,R_c2);

calidad=[calidad ;Lq,Rq, (corr2(R_c1,L_c1)>0.95 &&  corr2(R_c1_cor,L_c1)>0.95)   &&    (corr2(R_c2,L_c2)>0.95 && corr2(R_c2_cor,L_c2)>0.95)];
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

% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
