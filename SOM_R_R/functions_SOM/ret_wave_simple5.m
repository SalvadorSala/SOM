function[work_stim]=ret_wave_simple5(r,points,angle,x0,y0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generates the positions for a retinal wave type of stimuly
% as described in[3]
%Parameters:
% r............................................radius of the ret wave
% points.......................................n points in each step
% angle........................................angle for each set of points
% x0,y0........................................Initial position of the wave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

work_stim=[];


for i=1:length(r)

 x=r(i)*cos(angle)+x0;
 y=r(i)*sin(angle)+y0;

pick=randi(numel(x),points,1);

x=x(pick)+(normrnd(0,0.04,1,length(x(pick))));
y=y(pick)+(normrnd(0,0.04,1,length(y(pick))));
 index_1=(x<=1 & y<=1) & (x>=0 & y>=0) ;
 x=x(index_1);
 y=y(index_1);
 
work_stim=[work_stim, [x;y]];

% hold on,
% plot(x',y','.','MarkerSize',10)
% pause(0.1)
%  axis ([0 1 0 1])


end