clc 
clear all
set(0,'defaultfigurecolor','w');
select=3; %1for maximum C£¬2 for L2 norm error, 3for negative C
button=2; %1for isotropic,2for anisotropic
count=10; 
dt=3600;
maxtimestep=1000;
if(button==1)
%isotropic
Dxx=1e-6;
Dxy=0;
Dyy=1e-6;
w=2;
else
%anisotropic
Dxx=1e-6;
Dxy=0;
Dyy=1e-7;
w=1.2;
end
x0=30;
y0=30;
c0=1;
C2=Dxx*Dyy-Dxy^2;
C3=Dxx+Dyy;
e_min=1e-100;


tt=0.018*count:0.018*count:0.018*maxtimestep; %time
tt=count:count:maxtimestep; %time
results1=load('results_dx0.125_var1_h_d0.01.mat');
results2=load('results_dx1.8_dx0.125_var1_GI_q_d0.01.mat');
ntotal=160801; %160801,103041,40401,10201
maker_idx = 1:5:length(tt);

%%
qoi_s=results1.s1;
qoi_f1=results1.s2;
qoi_f2=results1.s3;
for i=1:maxtimestep/count
    t=count*dt*i;
    A1=2*t*Dyy+w^2; 
    A2=2*t*Dxx+w^2;
    A3=t*Dxy;
    C4=sqrt(4*t^2*C2+2*t*w^2*C3+w^4);
    c(i,:)=c0*w^2/C4*exp((-(qoi_s(i,:,2)-x0).^2*A1 -(qoi_s(i,:,3)-y0).^2*A2  +4.*(qoi_s(i,:,2)-x0).*(qoi_s(i,:,3)-y0)*A3)/(8*t^2*C2+4*w^2*t*C3+2*w^4));
    
    
    cmax_num_s(i)=max(qoi_s(i,:,4));
    cmax_num_f1(i)=max(qoi_f1(i,:,4));
    cmax_num_f2(i)=max(qoi_f2(i,:,4));

    cmin_num_s(i)=min(qoi_s(i,:,4));
    cmin_num_f1(i)=min(qoi_f1(i,:,4));
    cmin_num_f2(i)=min(qoi_f2(i,:,4));

    if(cmin_num_s(i)>=0)  
       cmin_num_s(i)=e_min;
    end
    if(cmin_num_f1(i)>=0)  
       cmin_num_f1(i)=e_min;
    end
    if(cmin_num_f2(i)>=0)  
        cmin_num_f2(i)=e_min;
    end
    e1_s(i)=norm(qoi_s(i,:,4)-c(i,:),1)/ntotal;
    e1_f1(i)=norm(qoi_f1(i,:,4)-c(i,:),1)/ntotal; 
    e1_f2(i)=norm(qoi_f2(i,:,4)-c(i,:),1)/ntotal;

    e2_s(i)=norm(qoi_s(i,:,4)-c(i,:),2)/ntotal;
    e2_f1(i)=norm(qoi_f1(i,:,4)-c(i,:),2)/ntotal;
    e2_f2(i)=norm(qoi_f2(i,:,4)-c(i,:),2)/ntotal;
%     e2_fm(i)=norm(qoi_fm(i,:,4)-c(i,:),2)/ntotal;

    einf_s(i)=norm(qoi_s(i,:,4)-c(i,:),inf);
    einf_f1(i)=norm(qoi_f1(i,:,4)-c(i,:),inf);
    einf_f2(i)=norm(qoi_f2(i,:,4)-c(i,:),inf);
    [M,I] = max(c(i,:));     
    cpeak_num_s(i)=qoi_s(i,I,4); 
    cpeak_num_f1(i)=qoi_f1(i,I,4);
    cpeak_num_f2(i)=qoi_f2(i,I,4);

%     cpeak_ana(i)=c0*w^2./C4;
    cpeak_ana(i)=M;
    ecpeak_s(i)=cpeak_num_s(i)-cpeak_ana(i);
    ecpeak_f1(i)=cpeak_num_f1(i)-cpeak_ana(i);
    ecpeak_f2(i)=cpeak_num_f2(i)-cpeak_ana(i);

%     ecpeak_s(i)=cmax_num_s(i)-cpeak_ana(i);;
%     ecpeak_f(i)=cmax_num_f(i)-cpeak_ana(i);;
end
if(select==1)
p1=plot(tt,ecpeak_s,'--s','color',[50,127,205]/255,'linewidth',1.2,'markersize',5,'MarkerIndices',maker_idx);
hold on
p2=plot(tt,ecpeak_f1,'--^','color',[45,141,154]/255,'linewidth',1.2,'markersize',3,'MarkerIndices',maker_idx);
hold on
p3=plot(tt,ecpeak_f2,'--o','color',[62,199,157]/255,'linewidth',1.2,'markersize',3,'MarkerIndices',maker_idx);
hold on
elseif(select==2)
%e2
p1=plot(tt,e2_s,'--s','color',[50,127,205]/255,'linewidth',1.2,'markersize',5,'MarkerIndices',maker_idx);
hold on
p2=plot(tt,e2_f1,'--^','color',[45,141,154]/255,'linewidth',1.2,'markersize',3,'MarkerIndices',maker_idx);
hold on
p3=plot(tt,e2_f2,'--o','color',[62,199,157]/255,'linewidth',1.2,'markersize',3,'MarkerIndices',maker_idx);
hold on
else   %negative concentration 
p1=semilogy(tt,abs(cmin_num_s),'--s','color',[50,127,205]/255,'linewidth',1.2,'markersize',5,'MarkerIndices',maker_idx);
hold on
p2=semilogy(tt,abs(cmin_num_f1),'--^','color',[45,141,154]/255,'linewidth',1.2,'markersize',3,'MarkerIndices',maker_idx);
hold on
p3=semilogy(tt,abs(cmin_num_f2),'--o','color',[62,199,157]/255,'linewidth',1.2,'markersize',3,'MarkerIndices',maker_idx);
hold on
end
%%
%--------------------------------------------------------------------------
qoi_s=results1.s4;
qoi_f1=results1.s5;
qoi_f2=results1.s6;
for i=1:maxtimestep/count
    t=count*dt*i;
    A1=2*t*Dyy+w^2; 
    A2=2*t*Dxx+w^2;
    A3=t*Dxy;
    C4=sqrt(4*t^2*C2+2*t*w^2*C3+w^4);
    c(i,:)=c0*w^2/C4*exp((-(qoi_s(i,:,2)-x0).^2*A1 -(qoi_s(i,:,3)-y0).^2*A2  +4.*(qoi_s(i,:,2)-x0).*(qoi_s(i,:,3)-y0)*A3)/(8*t^2*C2+4*w^2*t*C3+2*w^4));
    
    cmax_num_s(i)=max(qoi_s(i,:,4));
    cmax_num_f1(i)=max(qoi_f1(i,:,4));
    cmax_num_f2(i)=max(qoi_f2(i,:,4));

    cmin_num_s(i)=min(qoi_s(i,:,4));
    cmin_num_f1(i)=min(qoi_f1(i,:,4));
    cmin_num_f2(i)=min(qoi_f2(i,:,4));
    if(cmin_num_s(i)>=0)  
       cmin_num_s(i)=e_min;
    end
    if(cmin_num_f1(i)>=0)  
       cmin_num_f1(i)=e_min;
    end
    if(cmin_num_f2(i)>=0)  
      cmin_num_f2(i)=e_min;
    end

    e1_s(i)=norm(qoi_s(i,:,4)-c(i,:),1)/ntotal;
    e1_f1(i)=norm(qoi_f1(i,:,4)-c(i,:),1)/ntotal; 
    e1_f2(i)=norm(qoi_f2(i,:,4)-c(i,:),1)/ntotal;

    e2_s(i)=norm(qoi_s(i,:,4)-c(i,:),2)/ntotal;
    e2_f1(i)=norm(qoi_f1(i,:,4)-c(i,:),2)/ntotal;
    e2_f2(i)=norm(qoi_f2(i,:,4)-c(i,:),2)/ntotal;
%     e2_fm(i)=norm(qoi_fm(i,:,4)-c(i,:),2)/ntotal;

    einf_s(i)=norm(qoi_s(i,:,4)-c(i,:),inf);
    einf_f1(i)=norm(qoi_f1(i,:,4)-c(i,:),inf);
    einf_f2(i)=norm(qoi_f2(i,:,4)-c(i,:),inf);
    [M,I] = max(c(i,:));    
    cpeak_num_s(i)=qoi_s(i,I,4); 
    cpeak_num_f1(i)=qoi_f1(i,I,4);
    cpeak_num_f2(i)=qoi_f2(i,I,4);

    cpeak_ana(i)=c0*w^2./C4;
    ecpeak_s(i)=cpeak_num_s(i)-cpeak_ana(i);
    ecpeak_f1(i)=cpeak_num_f1(i)-cpeak_ana(i);
    ecpeak_f2(i)=cpeak_num_f2(i)-cpeak_ana(i);
%     ecpeak_s(i)=cmax_num_s(i)-cpeak_ana(i);;
%     ecpeak_f(i)=cmax_num_f(i)-cpeak_ana(i);;
end
if(select==1)
%×î´óÅ¨¶ÈÎó²î

p4=plot(tt,ecpeak_s,'s--','color',[249,201,102]/255,'linewidth',1.2,'markersize',5,'MarkerIndices',maker_idx);
hold on
p5=plot(tt,ecpeak_f1,'--^','color',[247,153,117]/255,'linewidth',1.2,'markersize',3,'MarkerIndices',maker_idx);
hold on
p6=plot(tt,ecpeak_f2,'--o','color',[247,116,146]/255,'linewidth',1.2,'markersize',3,'MarkerIndices',maker_idx);
hold on

elseif(select==2)
%e2
p4=plot(tt,e2_s,'s--','color',[249,201,102]/255,'linewidth',1.2,'markersize',5,'MarkerIndices',maker_idx);
hold on
p5=plot(tt,e2_f1,'--^','color',[247,153,117]/255,'linewidth',1.2,'markersize',3,'MarkerIndices',maker_idx);
hold on
p6=plot(tt,e2_f2,'--o','color',[247,116,146]/255,'linewidth',1.2,'markersize',3,'MarkerIndices',maker_idx);
hold on
else   %negative concentration
p4=semilogy(tt,abs(cmin_num_s),'s--','color',[249,201,102]/255,'linewidth',1.2,'markersize',5,'MarkerIndices',maker_idx);
hold on
p5=semilogy(tt,abs(cmin_num_f1),'--^','color',[247,153,117]/255,'linewidth',1.2,'markersize',3,'MarkerIndices',maker_idx);
hold on
p6=semilogy(tt,abs(cmin_num_f2),'--o','color',[247,116,146]/255,'linewidth',1.2,'markersize',3,'MarkerIndices',maker_idx);
hold on
end
%%
qoi_s=results2.s;
qoi_f1=results2.sh;
qoi_f2=results2.sh;
for i=1:maxtimestep/count
   t=count*dt*i;
    A1=2*t*Dyy+w^2; 
    A2=2*t*Dxx+w^2;
    A3=t*Dxy;
    C4=sqrt(4*t^2*C2+2*t*w^2*C3+w^4);
    c(i,:)=c0*w^2/C4*exp((-(qoi_s(i,:,2)-x0).^2*A1 -(qoi_s(i,:,3)-y0).^2*A2  +4.*(qoi_s(i,:,2)-x0).*(qoi_s(i,:,3)-y0)*A3)/(8*t^2*C2+4*w^2*t*C3+2*w^4));
 
    cmax_num_s(i)=max(qoi_s(i,:,4));
    cmax_num_f1(i)=max(qoi_f1(i,:,4));
    cmax_num_f2(i)=max(qoi_f2(i,:,4));

    cmin_num_s(i)=min(qoi_s(i,:,4));
    cmin_num_f1(i)=min(qoi_f1(i,:,4));
    cmin_num_f2(i)=min(qoi_f2(i,:,4));

    if(cmin_num_s(i)>=0)  
       cmin_num_s(i)=e_min;
    end
    if(cmin_num_f1(i)>=0)  
       cmin_num_f1(i)=e_min;
    end
    if(cmin_num_f2(i)>=0)  
      cmin_num_f2(i)=e_min;
    end

    e1_s(i)=norm(qoi_s(i,:,4)-c(i,:),1)/ntotal;
    e1_f1(i)=norm(qoi_f1(i,:,4)-c(i,:),1)/ntotal; 
    e1_f2(i)=norm(qoi_f2(i,:,4)-c(i,:),1)/ntotal;

    e2_s(i)=norm(qoi_s(i,:,4)-c(i,:),2)/ntotal;
    e2_f1(i)=norm(qoi_f1(i,:,4)-c(i,:),2)/ntotal;
    e2_f2(i)=norm(qoi_f2(i,:,4)-c(i,:),2)/ntotal;

    einf_s(i)=norm(qoi_s(i,:,4)-c(i,:),inf);
    einf_f1(i)=norm(qoi_f1(i,:,4)-c(i,:),inf);
    einf_f2(i)=norm(qoi_f2(i,:,4)-c(i,:),inf);
    [M,I] = max(c(i,:));     
    cpeak_num_s(i)=qoi_s(i,I,4); 
    cpeak_num_f1(i)=qoi_f1(i,I,4);
    cpeak_num_f2(i)=qoi_f2(i,I,4);

    cpeak_ana(i)=c0*w^2./C4;
    ecpeak_s(i)=cpeak_num_s(i)-cpeak_ana(i);
    ecpeak_f1(i)=cpeak_num_f1(i)-cpeak_ana(i);
    ecpeak_f2(i)=cpeak_num_f2(i)-cpeak_ana(i);
%     ecpeak_s(i)=cmax_num_s(i)-cpeak_ana(i);;
%     ecpeak_f(i)=cmax_num_f(i)-cpeak_ana(i);;
end
if(select==1)
% maximum concentartion 
p7 = plot(tt,ecpeak_s,'x--','color',[51,95,114]/255,'linewidth',1.2,'markersize',5,'MarkerIndices',maker_idx);
hold on
p8 = plot(tt,ecpeak_f1,'-*b','linewidth',1.5,'markersize',5,'MarkerIndices',maker_idx);
hold on
elseif(select==2) %e2
p7 = plot(tt,e2_s,'--x','color',[51,95,114]/255,'linewidth',1.2,'markersize',5,'MarkerIndices',maker_idx);
hold on
p8 = plot(tt,e2_f1,'-*b','linewidth',1.5,'markersize',5,'MarkerIndices',maker_idx);
hold on
else   %negative concentration 
   p7=semilogy(tt,abs(cmin_num_s),'--x','color',[51,95,114]/255,'linewidth',1.2,'markersize',5,'MarkerIndices',maker_idx);
   hold on
   p8=semilogy(tt,abs(cmin_num_f1),'-*b','linewidth',1.5,'markersize',5,'MarkerIndices',maker_idx);
   hold on 
end
% ---------------------------------------------------------------------------
if(select==1) 
     ylim([0,0.05]);
     % ylim([0,0.03]);
      % ylim([0,0.025]);
elseif(select==2)
    ylim([0,6.2e-6]);
elseif(select==3)
    ylim([1e-15,1]);
end

lgd1=legend([p7;p1;p2;p3;p4;p5;p6;p8],'h=1.8\Deltax','h=2\Deltax','h=3\Deltax','h=4\Deltax','h=5\Deltax','h=6\Deltax','h=7\Deltax','h_{GI}=5.02\Deltax');  %5.02,4.44, 3.38
set(lgd1,'FontSize',14,'FontName','Times New Roman','orientation','horizontal','NumColumns',2,'box','off');


if(select==1)
    ylabel('(C_{N,max}-C_{A,max})/C_0','FontSize',14,'FontName','Times New Roman ');
elseif(select==2)
    ylabel('L_2','FontSize',14,'FontName','Times New Roman ');
else
    ylabel('|C_{min}/C_0|','FontSize',14,'FontName','Times New Roman ');
end
xlabel('Time(hours)','FontSize',14,'FontName','Times New Roman ');
set(gca, 'FontSize',16,'FontName','Times New Roman ','FontWeight','normal','linewidth',1.5);
set(gcf,'Units','centimeters','Position',[1 1 14 10]);