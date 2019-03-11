function[work_stim]=ret_wave_simple2(r,points,angle,x0,y0)

work_stim=[];


for i=1:length(r)


    if r(i)>0.9
    



 x=r(i)*cos(angle)+x0;
 y=r(i)*sin(angle)+y0;
 index_1=(x<=1 & y<=1) & (x>=0 & y>=0) ;
 x=x(index_1);
 y=y(index_1);
 
pick=randi(numel(x),points,1);
 

 
if angle(1)==0
x=x(pick)-(-.2*(rand(1,points)-1));
y=y(pick)-(-.1*(rand(1,points)-1));
else
x=x(pick)+(-.2*(rand(1,points)-1));
y=y(pick)-(-.2*(rand(1,points)-1));
end


    else
        


angl=angle(randi(numel(angle),points,1));

if angle(1)==0
x=r(i)*cos(angl)+x0-(-.2*(rand(1,points)-1));
y=r(i)*sin(angl)+y0-(-.2*(rand(1,points)-1));
else
x=r(i)*cos(angl)+x0+(-.2*(rand(1,points)-1));
y=r(i)*sin(angl)+y0-(-.2*(rand(1,points)-1));
end

    end

work_stim=[work_stim, [x;y]];

% 
% 
% hold on,
% plot(x',y','.')
% pause(0.5)
%  axis ([0 1 0 1])
    
    
    
    
    
    
    
    
% 
%     
% 
% angle=0:0.001:pi/2;
% 
% 
% angl=angle(randi(numel(angle),points,1))
% 
%  
% 
% x=r(i)*cos(angl)+x0-(-.11*(rand(1,points)-1));
% 
% y=r(i)*sin(angl)+y0-(-.11*(rand(1,points)-1));
% 
% 
% 
% work_stim=[work_stim, [x;y]];
% 
% 
% 
% hold on,
% plot(x',y','.')
% pause(0.5)
% %  axis ([0 1 0 1])

end