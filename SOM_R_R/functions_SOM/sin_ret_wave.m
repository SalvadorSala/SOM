function[work_stim]=sin_ret_wave(r,points)

work_stim=[];




for i=1:length(r)

recta=@(x) -1*x+r(i)

if r(i)>1
    
% x=r(i)-1:0.001:1;
x=linspace(r(i)-1,1,points)
else
% x=0:0.001:r(i);
x=linspace(0,r(i),points)

end  

angle=recta(x);



x=x-(-.1*(rand(1,length(x))-1));

y=angle-(-.1*(rand(1,length(x))-1));


work_stim=[work_stim, [x;y]];
% 
% 
% 
% hold on,
% plot(x',y','.')
% pause(0.5)
% axis([0 1 0 1 ])



end