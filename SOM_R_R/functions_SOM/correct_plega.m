function[plega,cruce,girado]=correct_plega(c1,c2,nn,umbral)

%correct_plega calculates if the map is correctly undolded
%     
% INPUT
%     c1 - x coordinates of the neurons RF centers
%     c2 - y coordinates of the neurons RF centers
%     nn - number of neurons
%     umbral - threshold to determine an error
%
% OUTPUT, this determines the types of error the map can have
% plega=      1 correct unfolding    or  0 incorrect unfoldfing
% cruce=      1    map is twisted    or  0 map is not twisted
% girado=     1    map is rotated    or  0 map is aligned with the edeges
%
%
% AUTHOR: SSP & AJV
% ------------------------------------------------------------------------
index_1=c1>umbral;
a1=index_1(1,1);
b1=index_1(nn,1);
c1=index_1(1,nn);
d1=index_1(nn,nn);
reg1=a1+b1+c1+d1;

index_2=c2>umbral;
a2=index_2(1,1);
b2=index_2(nn,1);
c2=index_2(1,nn);
d2=index_2(nn,nn);
reg2=a2+b2+c2+d2;

% plega=(a1==d1 || a2==d2)==0
plega=(a1==d1 || a2==d2)==0 && (reg1~=2 || reg2~=2)==0 ;

cruce=(a1==d1 || a2==d2)==1 && (reg1~=2 || reg2~=2)==0 ;
girado=(reg1~=2 || reg2~=2);


end