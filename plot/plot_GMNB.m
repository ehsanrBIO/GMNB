load('../result_GMNB_GSE52260.mat')

%%FLNA
k = 15270;
%%EGR1
%k = 566;
%%NR4A1
%k = 6138;
%%ACTB
%k = 6675;
%%MYC
%k = 2111;
%%PKM2
%k = 11479;
%%EGR2
%k = 14014;
%%IL6ST
%k = 14202;
%%SEPT5
%k = 5377;
%%BATF3
%k = 6431;
%%COL1A2
%k = 10103;
%%ENO2
%k = 12694;

x= [0 .5 1 2 4 6 12 24 48 72]; 
y1=result{k}.est_r(2,:);
y2=result{k}.est_r(3,:);

e1=result{k}.var_r(2,:);
e2=result{k}.var_r(3,:);

lo1=y1-2.58*sqrt(e1);                 
hi1=y1+2.58*sqrt(e1);
lo2=y2-2.58*sqrt(e2);
hi2=y2+2.58*sqrt(e2);

in1=[x,x(end:-1:1)];
out1=[lo1,hi1(end:-1:1)];  
out2=[lo2,hi2(end:-1:1)];  

figure; hold on;
patch('XData',in1,'YData',out1,'FaceColor','blue','FaceAlpha',.3);
plot(x,y1,'b')

patch('XData',in1,'YData',out2,'FaceColor','red','FaceAlpha',.3);
plot(x,y2,'r')

load('../data/GSE52260_rc_human_timeSeriesRNAseq.mat') 

p1 = result{k}.PJT{2};
p2 = result{k}.PJT{3};

q1 = p1./(1-p1);
q2 = p2./(1-p2);

n1 = th0_data{k} ./q1;
n2 = th17_data{k} ./q2;

th0=plot(x,n1(1,:), 'bo', 'MarkerEdgeColor','b', 'MarkerFaceColor', 'b')
plot(x,n1(2:3,:), 'bo', 'MarkerEdgeColor','b', 'MarkerFaceColor', 'b')
th17=plot(x,n2(1,:), 'rd', 'MarkerEdgeColor','r', 'MarkerFaceColor','r')
plot(x,n2(2:3,:), 'rd', 'MarkerEdgeColor','r', 'MarkerFaceColor', 'r')

legend([th0 th17], 'Th0', 'Th17')
