set(0,'defaultfigurecolor','w');
clear
clc
count=5; 
maxtimestep=600;
button=1; %1 for standard deviation of Q,2 for optimal smoothing length h
tt=count/24:count/24:maxtimestep/24; 

results1=load('a00_dx1.5_var0.2.mat');
results2=load('a00_dx1.5_var0.4.mat');
results3=load('a00_dx1.5_var0.8.mat');
qoi_a1=results1.a;
qoi_a2=results2.a;
qoi_a3=results3.a;

results4 = load('results_dx1.5_var0.2_GI_q_realcov.mat');
results5 = load('results_dx1.5_var0.4_GI_q_realcov.mat');
results6 = load('results_dx1.5_var0.8_GI_q_realcov.mat');
qoi_f1=results4.sh;
qoi_f2=results5.sh;
qoi_f3=results6.sh;


for i=1:maxtimestep/count

    s1(i)= std(qoi_a1(i,:,4));
    s2(i)= std(qoi_a2(i,:,4));
    s3(i)= std(qoi_a3(i,:,4));

    h1(i) = mean(qoi_f1(i,:,8));
    h2(i) = mean(qoi_f2(i,:,8));
    h3(i) = mean(qoi_f3(i,:,8));
end
 nstart=1;
 v1= 0.889; %var0.2
 v2= 0.876; %var0.4
 v3= 0.851; %var0.8
tt1 = nstart*count/24*v1:count/24*v1:600/24*v1; %
tt2 = nstart*count/24*v2:count/24*v2:600/24*v2; %
tt3 = nstart*count/24*v3:count/24*v3:600/24*v3; %
maker_idx = 1:6:length(tt1);
if(button==1)
 p1=plot(tt1,s1(nstart:120),'--o','linewidth',1.8,'MarkerSize',5,'MarkerIndices',maker_idx);
 hold on
 p2=plot(tt2,s2(nstart:120),'--^','linewidth',1.8,'MarkerSize',5,'MarkerIndices',maker_idx);
 hold on
 p3=plot(tt3,s3(nstart:120),'--x','linewidth',1.8,'MarkerSize',7,'MarkerIndices',maker_idx);
 hold on
else
 p1=plot([0,tt1],[0.1875,h1],'--o','linewidth',1.8,'MarkerSize',5,'MarkerIndices',maker_idx);
 hold on
 p2=plot([0,tt2],[0.1875,h2],'--^','linewidth',1.8,'MarkerSize',5,'MarkerIndices',maker_idx);
 hold on
 p3=plot([0,tt3],[0.1875,h3],'--x','linewidth',1.8,'MarkerSize',7,'MarkerIndices',maker_idx);
 hold on
end

xlabel(' T = U t / I_{lnK}','FontSize',12,'FontName','Times New Roman');
if(button==1)
 ylabel('Standard deviation of Q_i','FontSize',10,'FontName','Times New Roman');
 ylim([0.02,0.12]);
else
ylabel('h_{GI}(m)','FontSize',10,'FontName','Times New Roman');
ylim([0.17,0.36]);
end
lgd1=legend([p1;p2;p3],'\sigma_{lnk}^2=0.2','\sigma_{lnk}^2=0.4','\sigma_{lnk}^2=0.8');
 set(lgd1,'FontSize',20,'FontName','Times New Roman','box','off');
xlim([0,23]);

set(gca, 'FontSize',20,'FontName','Times New Roman ','FontWeight','normal','linewidth',1.5);
set(gcf,'Units','centimeters','Position',[1 1 20 15]);



