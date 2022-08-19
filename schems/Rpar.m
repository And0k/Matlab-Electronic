function R= Rpar(Rin)
Rnom=[1
    1.2
    1.6
    1.8
    2
    2.2
    2.4
    2.7
    3
    3.3
    3.6
    3.9
    4.3
    4.7
    5.1
    5.6
    6.2
    6.8
    7.5
    8.2
    9.1
    10
    12
    15
    18
    20
    22
    24
    34
    30  
    36
    47
    51
    56
    68
    75
    82
    91
    100
    110
    120
    150
    180
    220
    240
    270
    300
    330
    360
    430
    470
    510
    620
    750
    820
    1000
    1100
    3300
    10000];
ind= find(Rnom>=Rin, 1, 'first'); Rnom= Rnom(ind:end);
R2ex=Rin.*Rnom./(Rnom-Rin);
n= 1;
for k=1:size(Rnom,1)
    ind= find(Rnom<R2ex(k),1,'last');
    if ~isempty(ind) && (ind<size(Rnom,1))
        R(n)= Rnom(ind);
        R(n+1)= Rnom(ind+1);        
    else
        R(n)= Rnom(end);
        R(n+1)= Rnom(end);
    end
    n= n+2;
end
%R2ex= repmat(Rnom', size(Rnom)) - repmat(R2ex,size(Rnom'));
%[qality ind]= min(abs(R2ex));
%R= [Rnom Rnom(ind)];
RnomX2= reshape([Rnom Rnom]',[1, 2*size(Rnom,1)]);
qality= (RnomX2.*R./(RnomX2+R)) - Rin;
b_ind= (qality<abs(Rin-Rnom(1)));
R= [RnomX2(b_ind); R(b_ind)]';
Rnom= RnomX2';
for k=1:size(ind,2)
    i_eq= find(Rnom(k)==R(:,2));
    R(i_eq(Rnom(i_eq)==R(k,2)), 2)= NaN;
end
b_ind= ~isnan(R(:,2));
R= R(b_ind,:);
[qality ind]= sort((R(:,1).*R(:,2)./(R(:,1)+R(:,2))) - Rin);
b_ind= (abs(qality)<(Rnom(1)-Rin));
R= [R(ind(b_ind),:); [Rnom(1) inf]];
R= [R [qality(b_ind); (Rnom(1)-Rin)]];