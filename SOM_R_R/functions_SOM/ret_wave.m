function[ret_wave]=ret_wave(vec_wave,size_wave,Rx,Ry)

sig1=0.1:0.1:0.6;sig2=0.2:0.1:0.7;
ret_wave=zeros(size_wave,size_wave,length(sig1));

Rxmat_pre=repmat(vec_wave,[size_wave,1]);
Lymat_pre=repmat(vec_wave',[1,size_wave]);



for zz=1:length(sig1)
sig1_aux=1/(2*sig1(zz)^2); 
 sig2_aux=1/(2*(sig2(zz))^2); 
%     sig2_aux=sig1_aux-2;
    
gauss1 = @(x,y) exp(-((Rxmat_pre-x).^2+(Lymat_pre-y).^2)*sig1_aux);
gauss2 = @(x,y) exp(-((Rxmat_pre-x).^2+(Lymat_pre-y).^2)*sig2_aux);
 
    
Base1=gauss1(Rx,Ry);
Base2=gauss2(Rx,Ry);
ret_wave(:,:,zz)=mat2gray(Base2-Base1);


end
% figure, 
% for u=1:size(ret_wave,3)
% imagesc(ret_wave(:,:,u))
% pause(0.2)
% end
% % 















end