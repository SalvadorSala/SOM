function[work_stim]=ret_wave_simple(r,points)

work_stim=[];
x0=0.1;
y0=0;  


for i=1:length(r)


    if r(i)>0.9
    
    
angle=0:0.001:pi/2;



 x=r(i)*cos(angle)+x0;
 y=r(i)*sin(angle)+y0;
 index_1=x<=1 & y<=1;
 x=x(index_1);
 y=y(index_1);
 
 
pick=randi(numel(x),points,1);
 

 

x=x(pick)-(-.1*(rand(1,points)-1));

y=y(pick)-(-.1*(rand(1,points)-1));





    else
        
angle=0:0.001:pi/2;


angl=angle(randi(numel(angle),points,1));

 

x=r(i)*cos(angl)+x0-(-.1*(rand(1,points)-1));

y=r(i)*sin(angl)+y0-(-.1*(rand(1,points)-1));
    end
    
    
    
%     if i>1
        
work_stim=[work_stim, [x;y]];
hold on,
plot(x',y','.')
pause(0.5)

%     end


    
    
    
    
    
    
    
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