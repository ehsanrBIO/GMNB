load('../result_GMNB_GSE52260.mat')

%%FLNA
k = 15270;
%k = 566;
%k = 6138;
%k = 6675;
%k = 6138;
%k = 2111;
%k = 11479;
%k = 5377;
x= [0 .5 1 2 4 6 12 24 48 72]; 


y1=result{k}.est_r(2,:);
y2=result{k}.est_r(3,:);
p1 = result{k}.PJT{2};
p2 = result{k}.PJT{3};
p1=1-p1;              
p2=1-p2;

v1=[]; for i=1:1000; v1=[v1; (nbinrnd(repmat(y1,3,1), p1))]; end
v2=[]; for i=1:1000; v2=[v2; (nbinrnd(repmat(y2,3,1), p2))]; end

mv1=mean(v1);
mv2=mean(v2);

e1=var(v1);
e2=var(v2);

lo1=mv1-2.58*sqrt(e1);                 
hi1=mv1+2.58*sqrt(e1);
lo2=mv2-2.58*sqrt(e2);
hi2=mv2+2.58*sqrt(e2);

in1=[x,x(end:-1:1)];
out1=[lo1,hi1(end:-1:1)];  
out2=[lo2,hi2(end:-1:1)];  

%figure; hold on;
patch('XData',in1,'YData',out1,'FaceColor','blue','FaceAlpha',.3);
plot(x,mv1,'b')

patch('XData',in1,'YData',out2,'FaceColor','red','FaceAlpha',.3);
plot(x,mv2,'r')

load('../data/GSE52260_rc_human_timeSeriesRNAseq.mat')

n1 = th0_data{k};
n2 = th17_data{k};

th0=plot(x,n1(1,:), 'bo', 'MarkerEdgeColor','b', 'MarkerFaceColor', 'b')
plot(x,n1(2:3,:), 'bo', 'MarkerEdgeColor','b', 'MarkerFaceColor', 'b')
th17=plot(x,n2(1,:), 'rd', 'MarkerEdgeColor','r', 'MarkerFaceColor','r')
plot(x,n2(2:3,:), 'rd', 'MarkerEdgeColor','r', 'MarkerFaceColor', 'r')

legend([th0 th17], 'Th0', 'Th17')
