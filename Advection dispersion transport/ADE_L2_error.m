set(0,'defaultfigurecolor','w');
clear
clc
count=5; 
maxtimestep=600;

% tt=0.03708*count:0.03708*count:0.03708*maxtimestep; 
tt=count/24:count/24:maxtimestep/24; 
 sa1=load('sa_dx6.0_var0.4_fdn.mat');


results1=load('results_dx1.5_var0.4_GI_q_realcov.mat');
results2=load('results_dx1.5_var0.4_smax.mat');

qoi_sa1=sa1.sa;
% qoi_sa2=sa2.sa;

qoi_s1=results1.s;
qoi_f1=results1.sh;
% % qoi_s2=results2.s;
qoi_f2=results2.sh;


dxq=0.25;
xq=0.5*dxq:dxq:40-0.5*dxq;
yq=0.5*dxq:dxq:15-0.5*dxq;
[X,Y]=meshgrid(xq,yq);
ngrid=9600;    

for i=1:maxtimestep/count
Fa1= scatteredInterpolant(qoi_sa1(i,:,2)',qoi_sa1(i,:,3)',qoi_sa1(i,:,6)','natural','nearest');
% Fa2= scatteredInterpolant(qoi_sa2(i,:,2)',qoi_sa2(i,:,3)',qoi_sa2(i,:,6)','natural','nearest');

Fs1 = scatteredInterpolant(qoi_s1(i,:,2)',qoi_s1(i,:,3)',qoi_s1(i,:,6)','natural','nearest');
Ff1 = scatteredInterpolant(qoi_f1(i,:,2)',qoi_f1(i,:,3)',qoi_f1(i,:,6)','natural','nearest');
Ff2 = scatteredInterpolant(qoi_f2(i,:,2)',qoi_f2(i,:,3)',qoi_f2(i,:,6)','natural','nearest');

% 
za1(i,:,:)=Fa1(X,Y);
% za2(i,:,:)=Fa2(X,Y);

zs1(i,:,:)=Fs1(X,Y);
zf1(i,:,:)=Ff1(X,Y);
zf2(i,:,:)=Ff2(X,Y);

end
%norm error
for i=1:maxtimestep/count
     ea(:,:)=za1(i,:,:); 
    
    es1(:,:)= zs1(i,:,:);
    ef1(:,:)= zf1(i,:,:);
    ef2(:,:)= zf2(i,:,:);
    
    e1_a(i)=norm(ea(:,:),1)/ngrid;
    e1_s1(i)=norm(es1(:,:)-ea(:,:),1)/ngrid;
    e1_f1(i)=norm(ef1(:,:)-ea(:,:),1)/ngrid; 
    e1_f2(i)=norm(ef2(:,:)-ea(:,:),1)/ngrid; 

    
    e2_a(i)=norm(ea(:,:),2)/ngrid;
    e2_s1(i)=norm(es1(:,:)-ea(:,:),2)/ngrid;
    e2_f1(i)=norm(ef1(:,:)-ea(:,:),2)/ngrid;
    e2_f2(i)=norm(ef2(:,:)-ea(:,:),2)/ngrid;

    RMSE_a(i)=norm(ea(:,:),2)/sqrt(ngrid);
    RMSE_s1(i)=norm(es1(:,:)-ea(:,:),2)/sqrt(ngrid);
    RMSE_f1(i)=norm(ef1(:,:)-ea(:,:),2)/sqrt(ngrid);
    
%     einf_s(i)=norm(es(:,:)-ea(:,:),inf);
%     einf_f(i)=norm(ef(:,:)-ea(:,:),inf);
    
%     pr(i)=e2_f(i)/e2_s(i);
end
%%
nstart=1;
 % v= 0.889; %var0.2
   v= 0.876; %var0.4
  % v= 0.851; %var0.8
tt = nstart*count/24*v:count/24*v:600/24*v; 
% tt = nstart*count:count:600; 
maker_idx = 1:4:length(tt);
p1=plot(tt,e2_s1(nstart:120),'-.+','color',[44,77,192]/255,'linewidth',1.8,'MarkerSize',5,'MarkerIndices',maker_idx);
hold on
p2=plot(tt,e2_f1(nstart:120),'-.o','color',[80,143,137]/255,'linewidth',1.8,'MarkerSize',5,'MarkerIndices',maker_idx);
hold on
p3=plot(tt,e2_f2(nstart:120),'-.x','color',[195,110,125]/255,'linewidth',1.8,'MarkerSize',5,'MarkerIndices',maker_idx);

xlim([0,23]);
ylim([0,6e-5]);

ylabel('L_2','FontSize',10,'FontName','Times New Roman');
xlabel(' T = U t / I_{lnK}','FontSize',12,'FontName','Times New Roman');
lgd1=legend([p1;p3;p2],'SPH-h_{min}','SPH-h_{max}','ADP-SPH');
set(lgd1,'FontSize',20,'FontName','Times New Roman','box','off');
set(gca, 'FontSize',20,'FontName','Times New Roman ','FontWeight','normal','linewidth',1.5);
set(gcf,'Units','centimeters','Position',[1 1 20 15]);
 
nstep=120;
e2_s = e2_s1(nstep);
e2_sh = e2_f1(nstep);
% e2_f2 = e2_f2(nstep);
% 
% RMSE_s = RMSE_s(120);
% RMSE_f = RMSE_f(120);
% RMSE_f_n = RMSE_f_n(120);

% e1_s= e1_s(120)/e1_a(120);
% e1_f = e1_f(120)/e1_a(120);
% e1_fd = e1_fd(120)/e1_a(120);
% e1_fd_n = e1_fd_n(120)/e1_a(120);

% e2_s_r = e2_s1(nstep)/e2_a(nstep);
% e2_sh_r = e2_f1(nstep)/e2_a(nstep);
% % e2_fd = e2_fd(120)/e2_a(120);
% % % e2_fd_n = e2_fd_n(120)/e2_a(120);
% 
% RMSE_s = RMSE_s(120)/RMSE_a(120);
% RMSE_f= RMSE_f(120)/RMSE_a(120);
% RMSE_fd= RMSE_fd(120)/RMSE_a(120);
% RMSE_fd_n = RMSE_fd_n(120)/RMSE_a(120);
%%