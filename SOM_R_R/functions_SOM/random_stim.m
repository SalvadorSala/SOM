function[cum_stim]= random_stim()
 cum_stim=[];    
rd=[rand(1,10000);rand(1,10000)]';
b=0.1;

point=[0.1:0.1:2;0.1:0.1:2]';
for r = [0.3]

    
    [idx1, dist1] = rangesearch(rd,point ,  0);
    [idx2, dist2] = rangesearch(rd, point, r+b);
 

    work=[0];
    for i=1:length(idx2)

  nce = setdiff([idx2{i}(1,2:end)], [idx1{i}(1,2:end)]); % Find non common elements
  
  pause(0.1)
   plot(rd(nce,1),rd(nce,2),'.'), axis([0,1 0,1])
  cum_stim=[cum_stim; rd(nce,1),rd(nce,2)];
    end

end




