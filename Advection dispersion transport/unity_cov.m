clc 
clear
set(0,'defaultfigurecolor','w');
button=1; %0for unity£¬1 for analysis
count=10; 
dt=3600;
maxtimestep=1200;
ntotal=19200; %160801,40401,10201,2601,76800
p=6;

if (button==0) %//////////////////////////////////////////////////////////////////////////
count_d=1*count;
for i=1:1:maxtimestep/count
    
    cov_temp=load(['a' num2str(count_d*i) '.dat']);
    cov(i,:,:)=cov_temp;
end
disp(i)
save('cov_var0.8_real.mat','cov');
disp('done')

else %//////////////////////////////////////////////////////////////////////////////////
results1=load('cov_var0.2_real.mat');
results2=load('cov_var0.4_real.mat');
results3=load('cov_var0.8_real.mat');
qoi_cov1=results1.cov;
qoi_cov2=results2.cov;
qoi_cov3=results3.cov;

for i=1:maxtimestep/count

    s1(1,i)= qoi_cov1(i,1);
    s2(1,i)= qoi_cov2(i,1);
    s3(1,i)= qoi_cov3(i,1);

    s1(2,i)= qoi_cov1(i,2);
    s2(2,i)= qoi_cov2(i,2);
    s3(2,i)= qoi_cov3(i,2);
end
nstart=1;
count=5;
 v1= 0.889; %var0.2
 v2= 0.876; %var0.4
 v3= 0.851; %var0.8
  tt=10:10:1200;
tt1 = nstart*count/24*v1:count/24*v1:600/24*v1; %
tt2 = nstart*count/24*v2:count/24*v2:600/24*v2; %
tt3 = nstart*count/24*v3:count/24*v3:600/24*v3; %

maker_idx = 1:6:length(tt1);

p1=plot(tt1,s1(1,nstart:120),'-o','color',[0,114,189]/255,'linewidth',1.8,'MarkerSize',5,'MarkerIndices',maker_idx);
hold on
p2=plot(tt2,s2(1,nstart:120),'-^','color',[217,83,25]/255,'linewidth',1.8,'MarkerSize',5,'MarkerIndices',maker_idx);
hold on
p3=plot(tt3,s3(1,nstart:120),'-x','color',[239,187,60]/255,'linewidth',1.8,'MarkerSize',7,'MarkerIndices',maker_idx);
hold on

p4=plot(tt1,s1(2,nstart:120),'--o','color',[0,114,189]/255,'linewidth',1.8,'MarkerSize',5,'MarkerIndices',maker_idx);
hold on
p5=plot(tt2,s2(2,nstart:120),'--^','color',[217,83,25]/255,'linewidth',1.8,'MarkerSize',5,'MarkerIndices',maker_idx);
hold on
p6=plot(tt3,s3(2,nstart:120),'--x','color',[239,187,60]/255,'linewidth',1.8,'MarkerSize',7,'MarkerIndices',maker_idx);
hold on

xlim([0,23]);
ylim([15,46]);
xlabel(' T = U t / I_{lnK}','FontSize',12,'FontName','Times New Roman');
ylabel('Covariance_{xx,yy}','FontSize',10,'FontName','Times New Roman');
lgd1=legend([p1;p2;p3;p4;p5;p6],'\Sigma_{xx,\sigma_{lnk}^2=0.2}','\Sigma_{xx,\sigma_{lnk}^2=0.4}','\Sigma_{xx,\sigma_{lnk}^2=0.8}','\Sigma_{yy,\sigma_{lnk}^2=0.2}','\Sigma_{yy,\sigma_{lnk}^2=0.4}','\Sigma_{yy,\sigma_{lnk}^2=0.8}');
set(lgd1,'FontSize',20,'FontName','Times New Roman','NumColumns',2,'box','off');
set(gca, 'FontSize',20,'FontName','Times New Roman ','FontWeight','normal','linewidth',1.5);
set(gcf,'Units','centimeters','Position',[1 1 20 15]);

end
  