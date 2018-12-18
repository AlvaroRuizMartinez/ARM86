function InvitroFtsZModels

clear all

% KronD subroutine can be downloaded from: https://www.mathworks.com/matlabcentral/fileexchange/26485-kronecker-delta

% Ruiz-Martinez et al., 2016 (Parameters)

k1=0.38;
k1m=0.01;
k2=0.79;
kemb=199.8;
kep=6.6;
deltat=4.05;
deltam=deltat*2;
kannp=6.6;
kem=kemb*exp(-deltat);
kannm=kemb*exp(-deltam);
kbp=3.5981;
kbm=199.8221;
Ub=0.175;
kmtb=2.7288;
khyd1=0.6681;
khyd2=0.143;
khyd3=0.112;

% Ruiz-Martinez et al., 2018 (Parameters)

x_F=log(sqrt(3/2));
x_B2=log(sqrt(3/(2*2)));
x_B3=log(sqrt(3/(2*3)));
F_F=1.0304+0.0193*(x_F)+0.06229*(x_F)^2+0.00476*(x_F)^3+0.00166*(x_F)^4+2.66e-6*(x_F)^7;
F_B2=1.0304+0.0193*(x_B2)+0.06229*(x_B2)^2+0.00476*(x_B2)^3+0.00166*(x_B2)^4+2.66e-6*(x_B2)^7;
F_B3=1.0304+0.0193*(x_B3)+0.06229*(x_B3)^2+0.00476*(x_B3)^3+0.00166*(x_B3)^4+2.66e-6*(x_B3)^7;

kbp11=4.0955;
kbp12=kbp11*(1/F_F+1/((2)^(1/3)*F_B2))*F_F/2;
kbp13=kbp11*(1/F_F+1/((3)^(1/3)*F_B3))*F_F/2;

Ub2=0.175;
Ubb2=(0.175+0.0405);
kbmi=199.9704;
kmtbi=2.1957;
khyd1i=0.6998;

%----

Cc=0.7;

CTOTAL=[1 2.5 5 7.5 10];

M=1;
%comptime=zeros(length(CTOTAL),M);
%avcomptime=zeros(length(CTOTAL),1);
%numbeq=zeros(length(CTOTAL),1);
Ninit=2;
for k=1:length(CTOTAL)
    
CT=CTOTAL(k);

tmax=40;
tmax2=40;
t=[0,tmax];
t2=[0,tmax2];

if CT <0.7
N=1;
else
for N=Ninit:1000

numbeq(k,1)=7+N;
numbeq2(k)=10;

%---

% Ruiz-Martinez et al., 2016

xinitb=zeros(N-1,1);
initx=[CT, 0, 0, 0, 0, 0, 0, 0, xinitb'];

for j=1:M

tic
[t,x]=ode45(@rhsFtsZ, t, [initx]);
 
end

if (x(length(t),N+7))>=1e-4
    toc;
    continue
else
    comptime(k,j)=toc;

L=zeros(length(t),1);
LL=zeros(length(t),1);
WW=zeros(length(t),1);
WB=zeros(length(t),1);
B=zeros(length(t),1);

for j=1:length(t)
    for i=1:N
        WB(j,1)=i*x(j,i+7)+WB(j,1);
    end
    L(j,1)=((CT-x(j,1)-x(j,2)-2*x(j,3)-3*x(j,4)-4*x(j,5)-5*x(j,6)-6*x(j,7))/(WB(j,1)));
    LL(j,1)=(2*x(j,3)+3*x(j,4)+4*x(j,5)+5*x(j,6)+6*x(j,7)+L(j,1)*WB(j,1))/(x(j,3)+x(j,4)+x(j,5)+x(j,6)+x(j,7)+WB(j,1));
end

for j=1:length(t)
    for i=1:N
        B(j,1)=x(j,i+7)+B(j,1);
    end
    W(j,1)=WB(j,1)/(B(j,1));
    WW(j,1)=(2*x(j,3)+3*x(j,4)+4*x(j,5)+5*x(j,6)+6*x(j,7)+L(j,1)*WB(j,1))/(2*x(j,3)+3*x(j,4)+4*x(j,5)+5*x(j,6)+6*x(j,7)+L(j,1)*B(j,1));
end

figure(1)
hold on
plot(t,x(:,1)+x(:,2),'b','Linewidth',2)

xlabel(['Time (s)'])
ylabel(['Concentration of monomers (',char(181),'M)'])
    set(gca,'FontSize',18,'fontWeight','bold')
    set(findall(gcf,'type','text'),'FontSize',18,'fontWeight','bold')
    legend('Ruiz-Martinez et al., 2016')
        set(gca,'box','off')
    set(gcf,'Position',[80 80 800 600])
%set(gcf,'Units','inches');
%screenposition = get(gcf,'Position');
%set(gcf,...
%    'PaperPosition',[0 0 screenposition(3:4)],...
%    'PaperSize',[screenposition(3:4)]);

figure(2)
hold on
plot(t,LL,'b','Linewidth',2)
xlabel(['Time (s)'])
ylabel(['Average length'])
    set(gca,'FontSize',18,'fontWeight','bold')
    set(findall(gcf,'type','text'),'FontSize',18,'fontWeight','bold')
    legend('Ruiz-Martinez et al., 2016')
        set(gca,'box','off')
    set(gcf,'Position',[80 80 800 600])
%set(gcf,'Units','inches');
%screenposition = get(gcf,'Position');
%set(gcf,...
%    'PaperPosition',[0 0 screenposition(3:4)],...
%    'PaperSize',[screenposition(3:4)]);
    
figure(3)
hold on
plot(t,WW,'b','Linewidth',2)
xlabel(['Time (s)'])
ylabel(['Average number of filaments per bundle'])
    set(gca,'FontSize',18,'fontWeight','bold')
    set(findall(gcf,'type','text'),'FontSize',18,'fontWeight','bold')
    legend('Ruiz-Martinez et al., 2016')
        set(gca,'box','off')
    set(gcf,'Position',[80 80 800 600])
%set(gcf,'Units','inches');
%screenposition = get(gcf,'Position');
%set(gcf,...
%    'PaperPosition',[0 0 screenposition(3:4)],...
%    'PaperSize',[screenposition(3:4)]);
    
CC(1)=x(length(t),3)+x(length(t),4)+x(length(t),5)+x(length(t),6)+x(length(t),7)+x(length(t),8);
for i=1:N
    NN(i)=i;    
end
for i=1:N-1
    CC(i+1)=x(length(t),i+8);
end

figure(100)
bar([1:N],CC/CC(1))
xlabel(['Time (s)'])
ylabel(['Number of filaments per bundle'])
    set(gca,'FontSize',18,'fontWeight','bold')
    set(findall(gcf,'type','text'),'FontSize',18,'fontWeight','bold')
    legend('Ruiz-Martinez et al., 2016')
        set(gca,'box','off')
    set(gcf,'Position',[80 80 800 600])
%set(gcf,'Units','inches');
%screenposition = get(gcf,'Position');
%set(gcf,...
%    'PaperPosition',[0 0 screenposition(3:4)],...
%    'PaperSize',[screenposition(3:4)]);
    
Ninit=N;
    break
end
end
end

avcomptime(k)=sum(comptime(k,:))/M;

% Ruiz-Martinez et al., 2018

initx2=[CT, 0, 0, 0, 0, 0, 0, 0, 0, 0];
for j=1:M
tic
[t2,x2]=ode45(@rhsFtsZimp, t2, [initx2]);
comptime2(k,j)=toc;
end
avcomptime2(k)=sum(comptime2(k,:))/M;

LL2=zeros(length(t2),1);
Length2=zeros(length(t2),1);
WW2=zeros(length(t2),1);
WWB=zeros(length(t2),1);


figure(10)
hold on
plot(t2,x2(:,1)+x2(:,2),'r','Linewidth',2)
xlabel(['Time (s)'])
ylabel(['Concentration of monomers (',char(181),'M)'])
    set(gca,'FontSize',18,'fontWeight','bold')
    set(findall(gcf,'type','text'),'FontSize',18,'fontWeight','bold')
    legend('Ruiz-Martinez et al., 2018')
        set(gca,'box','off')
    set(gcf,'Position',[80 80 800 600])
%set(gcf,'Units','inches');
%screenposition = get(gcf,'Position');
%set(gcf,...
%    'PaperPosition',[0 0 screenposition(3:4)],...
%    'PaperSize',[screenposition(3:4)]);

for j=1:length(t2)
    LL2(j,1)=x2(j,10)./(x2(j,5)+2*x2(j,6)+3*x2(j,7)+x2(j,9));
    Length2(j,1)=(2*x2(j,3)+3*x2(j,4)+x2(j,10))./(x2(j,3)+x2(j,4)+x2(j,5)+2*x2(j,6)+3*x2(j,7)+x2(j,9));
end

figure(20)
plot(t2,Length2,'r','Linewidth',2)
hold on
axis([0 tmax2 0 35])
xlabel(['Time (s)'])
ylabel(['Average length'])
    set(gca,'FontSize',18,'fontWeight','bold')
    set(findall(gcf,'type','text'),'FontSize',18,'fontWeight','bold')
    legend('Ruiz-Martinez et al., 2018')
        set(gca,'box','off')
    set(gcf,'Position',[80 80 800 600])
%set(gcf,'Units','inches');
%screenposition = get(gcf,'Position');
%set(gcf,...
%    'PaperPosition',[0 0 screenposition(3:4)],...
%    'PaperSize',[screenposition(3:4)]);

for j=1:length(t2)
    WW2(j,1)=(2*x2(j,3)+3*x2(j,4)+x2(j,10))./(2*x2(j,3)+3*x2(j,4)+x2(j,10).*((x2(j,5)+x2(j,6)+x2(j,7)+x2(j,8))./(x2(j,5)+2*x2(j,6)+3*x2(j,7)+x2(j,9))));
end
figure(30)
plot(t2,WW2,'r','Linewidth',2)
hold on
axis([0 tmax2 1 2])
xlabel(['Time (s)'])
ylabel(['Average number of filaments per bundle'])
    set(gca,'FontSize',18,'fontWeight','bold')
    set(findall(gcf,'type','text'),'FontSize',18,'fontWeight','bold')
    legend('Ruiz-Martinez et al., 2018')
        set(gca,'box','off')
    set(gcf,'Position',[80 80 800 600])
%set(gcf,'Units','inches');
%screenposition = get(gcf,'Position');
%set(gcf,...
%    'PaperPosition',[0 0 screenposition(3:4)],...
%    'PaperSize',[screenposition(3:4)]);

figure(40)
hold on
plot(t2,x2(:,9),'r','Linewidth',2)
axis([0 tmax2 0 0.15])
xlabel(['Time (s)'])
ylabel(['C_{wb}^f (',char(181),'M)'])
    set(gca,'FontSize',18,'fontWeight','bold')
    set(findall(gcf,'type','text'),'FontSize',18,'fontWeight','bold')
    legend('Ruiz-Martinez et al., 2018')
        set(gca,'box','off')
    set(gcf,'Position',[80 80 800 600])
%set(gcf,'Units','inches');
%screenposition = get(gcf,'Position');
%set(gcf,...
%    'PaperPosition',[0 0 screenposition(3:4)],...
%    'PaperSize',[screenposition(3:4)]);

figure(50)
hold on
plot(t2,x2(:,10),'r','Linewidth',2)
xlabel(['Time (s)'])
ylabel(['C_{fb}^m (',char(181),'M)'])
    set(gca,'FontSize',18,'fontWeight','bold')
    set(findall(gcf,'type','text'),'FontSize',18,'fontWeight','bold')
    legend('Ruiz-Martinez et al., 2018')
        set(gca,'box','off')
    set(gcf,'Position',[80 80 800 600])
%set(gcf,'Units','inches');
%screenposition = get(gcf,'Position');
%set(gcf,...
%    'PaperPosition',[0 0 screenposition(3:4)],...
%    'PaperSize',[screenposition(3:4)]);

for j=1:length(t2)
    WWB(j,1)=x2(j,9)./x2(j,8);
end

figure(60)
plot(t2,WWB,'r','Linewidth',2)
hold on
axis([0 tmax2 0 10])
xlabel(['Time (s)'])
ylabel(['Filaments per wide bundle'])
    set(gca,'FontSize',18,'fontWeight','bold')
    set(findall(gcf,'type','text'),'FontSize',18,'fontWeight','bold')
    legend('Ruiz-Martinez et al., 2018')
        set(gca,'box','off')
    set(gcf,'Position',[80 80 800 600])
%set(gcf,'Units','inches');
%screenposition = get(gcf,'Position');
%set(gcf,...
%    'PaperPosition',[0 0 screenposition(3:4)],...
%    'PaperSize',[screenposition(3:4)]);

end


figure(200)
hold on
plot(CTOTAL,avcomptime,'b',CTOTAL,avcomptime2,'r','Linewidth',2)
xlabel(['Total Concentration (',char(181),'M)'])
ylabel(['Computational Time (s)'])
    set(gca,'FontSize',18,'fontWeight','bold')
    set(findall(gcf,'type','text'),'FontSize',18,'fontWeight','bold')
    legend('Ruiz-Martinez et al., 2016', 'Ruiz-Martinez et al., 2018')
    set(gca,'box','off')
    set(gcf,'Position',[80 80 800 600])
%set(gcf,'Units','inches');
%screenposition = get(gcf,'Position');
%set(gcf,...
%    'PaperPosition',[0 0 screenposition(3:4)],...
%    'PaperSize',[screenposition(3:4)]);
%    print('CompTime', '-dpdf', '-r300');
    
figure(300)
hold on
plot(CTOTAL,numbeq,'b',CTOTAL,numbeq2,'r','Linewidth',2)

xlabel(['Total Concentration (',char(181),'M)'])
ylabel(['Number of ODEs'])
    set(gca,'FontSize',18,'fontWeight','bold')
    set(findall(gcf,'type','text'),'FontSize',18,'fontWeight','bold')
    legend('Ruiz-Martinez et al., 2016', 'Ruiz-Martinez et al., 2018')
    set(gca,'box','off')
    set(gcf,'Position',[80 80 800 600])
%set(gcf,'Units','inches');
%screenposition = get(gcf,'Position');
%set(gcf,...
%    'PaperPosition',[0 0 screenposition(3:4)],...
%    'PaperSize',[screenposition(3:4)]);
%    print('NumbODEs', '-dpdf', '-r300');


% Ruiz-Martinez et al., 2016 (ODEs)

    function dxdt=rhsFtsZ(t,x)
        
        WB=0;
        for i=1:N
            WB=i*x(i+7)+WB;
        end
        L=((CT-x(1)-x(2)-2*x(3)-3*x(4)-4*x(5)-5*x(6)-6*x(7))/(WB));
        if L>1 
            KBM=kbm*exp(-(L-1)*Ub);
        else
            KBM=kbm;
        end
        
        sum_hyd_3=0;
        sum_hyd_2a=0;
        sum_hyd_2b=0;
        if N>1
            sum_hyd_2a=x(9)+sum_hyd_2a;
        end  
        if N>2
            sum_hyd_2b=x(10)+sum_hyd_2b;
        end  
        for i=1:N-1
            sum_hyd_3=x(i+8)+sum_hyd_3;
        end
        dxdt(1)=-k1*x(1)+k1m*x(2)+((khyd1+khyd2)*x(8)+khyd2*(sum_hyd_2a+sum_hyd_2b)+khyd3*sum_hyd_3)*(CT-x(1)-x(2))/(CT-0.7);
        sum_el_olig_p=0;
        sum_el_olig_m=0;
        sum_mb=0;
        for i=1:6
            sum_el_olig_p=x(i+2)+sum_el_olig_p;
        end
        for i=1:5
            sum_el_olig_m=x(i+3)+sum_el_olig_m;
        end
        for i=1:N-1
            sum_mb=x(i+8)+sum_mb;
        end
        dxdt(2)=k1*x(1)-k1m*x(2)-2*k2*(x(2))^2+2*kemb*x(3)-kep*x(2)*sum_el_olig_p+kem*sum_el_olig_m-kmtb*x(2)*sum_mb;
        dxdt(3)=k2*(x(2))^2-kemb*x(3)-kep*x(2)*x(3)+kem*x(4);
        for j=1:3
            dxdt(j+3)=kep*x(2)*(x(j+2)-x(j+3))+kem*(x(j+4)-x(j+3));
        end
        dxdt(7)=kep*x(2)*(x(6)-x(7))+kem*(-x(7));
        sum_bun_p=0;
        sum_bun_m=0;
        for i=2:N-1
            sum_bun_p=x(i+7)+sum_bun_p;
        end
        for i=3:N
            sum_bun_m=x(i+7)+sum_bun_m;
        end
        BUND=0;
        if N>1
            BUND=-kbp*x(8)*(2*x(8)+sum_bun_p)+KBM*(2*x(9)+sum_bun_m)+BUND;
        end
        dxdt(8)=kep*x(2)*x(7)-kannp*(x(8))^2+kannm*x(8)+BUND+khyd2*x(8)*(CT-x(1)-x(2))/(CT-0.7);
        
        if N>1
        %EVEN i
        if N/2-floor(N/2)==0
            K=N-2;
        else
            K=N-1;
        end
        for i=2:2:K
            sum_bun_p2=0;
            sum_bun_p3=0;
            sum_bun_m2=0;
            for j=1:i-1
                sum_bun_p2=(1-KronD(j,i-j))/2*x(j+7)*x(i-j+7)+sum_bun_p2;
            end 
            for j=1:N-i
                sum_bun_p3=(x(j+7)+x(i+7)*KronD(i,j))+sum_bun_p3;
            end 
            for j=i+1:N
                sum_bun_m2=(1+KronD(2*i,j))*x(j+7)+sum_bun_m2;
            end 
            dxdt(i+7)=-i/2*KBM*x(i+7)+kbp*((x(i/2+7))^2+sum_bun_p2)-kbp*x(i+7)*sum_bun_p3+KBM*sum_bun_m2+khyd2*x(i+7)*KronD(i,2)*(CT-x(1)-x(2))/(CT-0.7);
        end
        %ODD i
        if N/2-floor(N/2)==0
            K=N-1;
        else
            K=N-2;
        end
        for i=3:2:K
            sum_bun_p4=0;
            sum_bun_p5=0;
            sum_bun_m3=0;
            for j=1:i-1
                sum_bun_p4=x(j+7)*x(i-j+7)+sum_bun_p4;
            end 
            for j=1:N-i
                sum_bun_p5=(x(j+7)+x(i+7)*KronD(i,j))+sum_bun_p5;
            end 
            for j=i+1:N
                sum_bun_m3=(1+KronD(2*i,j))*x(j+7)+sum_bun_m3;
            end 
            dxdt(i+7)=-(i-1)/2*KBM*x(i+7)+kbp/2*sum_bun_p4-kbp*x(i+7)*sum_bun_p5+KBM*sum_bun_m3+khyd2*x(i+7)*KronD(i,3)*(CT-x(1)-x(2))/(CT-0.7);
        end

        if N/2-floor(N/2)==0
            sum_bun_p6=0;
            for j=1:N-1
                sum_bun_p6=(1-KronD(j,N-j))/2*x(j+7)*x(N-j+7)+sum_bun_p6;
            end
            dxdt(N+7)=-N/2*KBM*x(N+7)+kbp*((x(N/2+7))^2+sum_bun_p6)+khyd2*x(N+7)*KronD(N,2)*(CT-x(1)-x(2))/(CT-0.7);
        else
            sum_bun_p7=0;
            for j=1:N-1
                sum_bun_p7=x(j+7)*x(N-j+7)+sum_bun_p7;
            end
            dxdt(N+7)=-(N-1)/2*KBM*x(N+7)+1/2*kbp*sum_bun_p7+khyd2*x(N+7)*KronD(N,3)*(CT-x(1)-x(2))/(CT-0.7);
        end
        end
        dxdt=dxdt';
    end

% Ruiz-Martinez et al., 2018 (ODEs)

    function dxdt2=rhsFtsZimp(t2,x2)
           
        L2=x2(10)/(x2(5)+2*x2(6)+3*x2(7)+x2(9));
        if L2>1 
            KBM2a=kbmi*exp(-(L2-1)*Ub2);
            KBM2b=kbmi*exp(-(L2-1)*Ubb2);
        else
            KBM2a=kbmi;
            KBM2b=kbmi;
        end
        W2=x2(9)/x2(8);
        xx=log(sqrt(3/(2*W2)));
        if W2>1.5e-4
            KBP1W=(1/F_F+1/(W2^(1/3)*(1.0304+0.0193*xx+0.06229*xx^2+0.00476*xx^3+0.00166*xx^4+2.66e-6*xx^7)))*F_F/2*kbp11;
            KBPWW=(2/(W2^(1/3)*(1.0304+0.0193*xx+0.06229*xx^2+0.00476*xx^3+0.00166*xx^4+2.66e-6*xx^7)))*F_F/2*kbp11;
            KBMW=kbmi*exp(-(26)*Ubb2);
            kmtbW=kmtbi;
            khyd3W=khyd3;
        else
            KBP1W=0;
            KBPWW=0;
            KBMW=0;
            kmtbW=0;
            khyd3W=0;
        end
            
           dxdt2(1)=-k1*x2(1)+k1m*x2(2)+(khyd1i*x2(5)+khyd2*(x2(5)+x2(6)+x2(7))+khyd3*(x2(6)+x2(7))+khyd3W*(x2(8)))*(CT-x2(1)-x2(2))/(CT-Cc);
           dxdt2(2)=k1*x2(1)-k1m*x2(2)-2*k2*(x2(2))^2+2*kemb*x2(3)-kep*x2(2)*(x2(3)+x2(4)+x2(5))+kem*(x2(4)+x2(5))-kmtbi*(x2(6)+x2(7))*x2(2)-kmtbW*x2(8)*x2(2);
           dxdt2(3)=k2*(x2(2))^2-kemb*(x2(3))-kep*x2(2)*(x2(3))+kem*x2(4);
           dxdt2(4)=-kem*x2(4)-kep*x2(2)*x2(4)+kep*x2(2)*(x2(3));
           dxdt2(5)=kep*x2(2)*x2(4)-kannp*(x2(5))^2+kannm*x2(5)-2*kbp11*(x2(5))^2-kbp12*x2(6)*x2(5)-kbp13*x2(7)*x2(5)-KBP1W*x2(8)*x2(5)+KBM2a*(2*x2(6))+KBM2b*exp(-0.0405)*(x2(7))+KBMW*exp(-0.0405)*x2(8)+khyd2*x2(5)*(CT-x2(1)-x2(2))/(CT-Cc);
           dxdt2(6)=-KBM2a*x2(6)+kbp11*(x2(5))^2-kbp12*x2(6)*(x2(5))+KBM2b*exp(-0.0405)*(x2(7))+khyd2*x2(6)*(CT-x2(1)-x2(2))/(CT-Cc);
           dxdt2(7)=-KBM2b*exp(-0.0405)*x2(7)+kbp12*x2(5)*x2(6)-kbp13*x2(7)*(x2(5))+khyd2*x2(7)*(CT-x2(1)-x2(2))/(CT-Cc);
           dxdt2(8)=kbp13*x2(7)*x2(5)-KBPWW*(x2(8))^2;
           dxdt2(9)=kbp13*x2(7)*x2(5)*4+((KBP1W*x2(5)*x2(8)-KBMW*exp(-0.0405)*x2(8)));
           dxdt2(10)=4*kep*x2(2)*x2(4)-1*((khyd1i+khyd2)*x2(5)+khyd3*(x2(6)+x2(7))+khyd3W*(x2(8))+khyd2*(x2(6)+x2(7)))*(CT-x2(1)-x2(2))/(CT-Cc)+kmtbi*(x2(6)+x2(7))*x2(2)+kmtbW*x2(8)*x2(2)+kep*x2(2)*x2(5)-kem*x2(5);
           dxdt2=dxdt2';
    end

end
