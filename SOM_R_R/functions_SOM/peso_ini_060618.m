function[pesoi]=peso_ini_060618(xmat_pos,ymat_pos,pesos)

for i=1:length(pesos)

 pesoi(i,:)=[sum((xmat_pos(:)'.*reshape(pesos(:,:,i),[1,length(pesos)])))./sum(reshape(pesos(:,:,i),[1,length(pesos)])),...
  sum((ymat_pos(:)'.*reshape(pesos(:,:,i),[1,length(pesos)])))./sum(reshape(pesos(:,:,i),[1,length(pesos)]))];

end
end
