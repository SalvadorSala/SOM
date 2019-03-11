function[stim]=ret_wave_STIM(nstim,ret_wave)

size_wave=size(ret_wave,2);

for q=1:size(ret_wave,3)
coord=[];

for i=1:nstim
  
   row= randi([1,size_wave],1);
   col=randi([1,size_wave],1);
    
  if ret_wave(col,row,q)>rand(1)
     coord=[coord, [col;row]];
     
  else
      display('nostim')
      
  end 
end

 hold on, plot(coord(2,:),coord(1,:),'o'), axis([1 length(0:0.001:1) 1 length(0:0.001:1)])
      pause(1)

end





end