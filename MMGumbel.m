a=zeros(2000,6);
options = optimset('fmincon');
options = optimset(options,'Display','off','TolCon',10^-12,'TolFun',10^-9,'TolX',10^-12);
RMSEstep = zeros(1000,1);
lower = [-25;-25;-25;-25;-25;-25;-25;-25;0.00001;0.00001]; 
upper = [25;25;25;25;25;25;25;25;0.99999;0.99999];
datestr(now,31)

theta0 = [0.2;0;0;0.2;0;0;-0.3;-0.3;0.9;0.9];
RNMSEstepL= zeros(1000,1);
RNMSEstepU= zeros(1000,1);
SAVE3=zeros(1000,4);
LL15= zeros(1000,1);
LL16= zeros(1000,1);
kappa = [0.5,0.9];

parfor i=1:1000
n=5000;%nnnnnnnnnnnnnnnnnnnnnnnnnnn=1000 2000 5000
[u,v,thetaU,thetaL]= Simstep(n);
[ kappa11] = fmincon('MSJC',theta0,[],[],[],[],lower,upper,[],options,[u,v],kappa);
[LL15(i), QthetaU,QthetaL] = Copy_of_MSJC(kappa11,[u,v],kappa);
NMSEU=0;NMSEL=0;
for j=1:n
        NMSEU=(QthetaU(j)-thetaU(j)).^2+NMSEU;
        NMSEL=(QthetaL(j)-thetaL(j)).^2+NMSEL;
end
RNMSEstepU(i) = NMSEU/n;
RNMSEstepL(i) = NMSEL/n;

[ kappa12] = fmincon('otherMSJC',theta0,[],[],[],[],lower,upper,[],options,[u,v],kappa);
[LL16(i)] =  otherMSJC(kappa12,[u,v],kappa);
NMSEU=0;NMSEL=0;


% fprintf(' \n  %1f\t%1f\t%1f\t%1f\t%1f',i,RNMSEstepU(i) ,RNMSEstepL(i),LL15(i) ,LL16(i));
% fprintf(' \n 3  %1f	',i);
end

SAVE3(:,1) = RNMSEstepU(:,1);
SAVE3(:,2) = RNMSEstepL(:,1);
SAVE3(:,3) = LL15(:,1);
SAVE3(:,4) = LL16(:,1);
csvwrite('SAVE3.csv',SAVE3)


theta0 = [0.2;0;0;0.2;0;0;-0.3;-0.3;0.9;0.9];
RNMSEotherL= zeros(1000,1);
RNMSEotherU= zeros(1000,1);
otherTAL1= zeros(1000,1);
otherTAL2= zeros(1000,1);
otherTAU1= zeros(1000,1);
otherTAU2= zeros(1000,1);
otherLLLL= zeros(1000,1);
LL17= zeros(1000,1);
LL18= zeros(1000,1);
SAVE4=zeros(1000,6);
otherRNU= zeros(1000,1);
otherRNL= zeros(1000,1);
kappa = [0.5,0.9];
parfor i=1:1000
n=5000;%nnnnnnnnnnnnnnnnnnnnnnnnnnn=1000 2000 5000
[u,v,thetaU,thetaL,state]= Simother(n);
[ kappa11] = fmincon('MSJC',theta0,[],[],[],[],lower,upper,[],options,[u,v],kappa);
[LL17(i), QthetaU,QthetaL] = Copy_of_MSJC(kappa11,[u,v],kappa);
NMSEU=0;NMSEL=0;
for j=1:n
        NMSEU=(QthetaU(j)-thetaU(j)).^2+NMSEU;
        NMSEL=(QthetaL(j)-thetaL(j)).^2+NMSEL;
end
RNMSEotherL(i) = NMSEU/n;
RNMSEotherU(i) = NMSEL/n;

[ kappa12] = fmincon('otherMSJC',theta0,[],[],[],[],lower,upper,[],options,[u,v],kappa);
[LL18(i), otherTAL1,otherTAL2,otherTAU1,otherTAU2] =  CopyotherMSJC(kappa12,[u,v],kappa);
NMSEU=0;NMSEL=0;

if kappa12(1)>kappa12(7)
    
    otherLLLL=otherTAL1;
    otherTAL1=otherTAL2;
    otherTAL2=otherLLLL;
end
if kappa12(4)>kappa12(8)
    
    otherLLLL=otherTAU1;
    otherTAU1=otherTAU2;
    otherTAU2=otherLLLL;
end
otherQthetaU = zeros(n,1);
otherQthetaL = zeros(n,1);
for j=1:n
    if state(j)>0.5
        otherQthetaU(j) = otherTAU2(j);
        otherQthetaL(j) = otherTAL2(j);
        NMSEU=(otherQthetaU(j)-thetaU(j)).^2+NMSEU;
        NMSEL=(otherQthetaL(j)-thetaL(j)).^2+NMSEL;
    else
        otherQthetaU(j) = otherTAU1(j);
        otherQthetaL(j) = otherTAL1(j);
        NMSEU=(otherQthetaU(j)-thetaU(j)).^2+NMSEU;
        NMSEL=(otherQthetaL(j)-thetaL(j)).^2+NMSEL;
   
    end
end

otherRNU(i) = NMSEU/n;
otherRNL(i) = NMSEL/n;
% fprintf(' \n  %1f\t%1f\t%1f\t%1f\t%1f\t%1f\t%1f',i,RNMSEotherL(i) ,RNMSEotherU(i),otherRNU(i) ,otherRNL(i),LL17(i),LL18(i));
% fprintf(' \n 3  %1f	',i);
end

SAVE4(:,1) = RNMSEotherL(:,1);
SAVE4(:,2) = RNMSEotherU(:,1);
SAVE4(:,3) = otherRNU(:,1);
SAVE4(:,4) = otherRNL(:,1);
SAVE4(:,5) = LL17(:,1);
SAVE4(:,6) = LL18(:,1);

csvwrite('SAVE4.csv',SAVE4)


RNMSEAR1U= zeros(1000,1);
RNMSEAR1L= zeros(1000,1);
LL11= zeros(1000,1);
LL12= zeros(1000,1);
theta0 = [0.2;0.1;0;0.2;0;0;0.3;0.3;0.9;0.9];
kappa = [0.2,0.2];
SAVE1=zeros(1000,4);

% 
parfor i=1:1000
n=1000;%nnnnnnnnnnnnnnnnnnnnnnnnnnn=1000 2000 5000
[u,v,thetaU,thetaL]= SimAR1(n);
[ kappa11] = fmincon('MSJC',theta0,[],[],[],[],lower,upper,[],options,[u,v],kappa);

[LL11(i), QthetaU,QthetaL] = Copy_of_MSJC(kappa11,[u,v],kappa);
NMSEU=0;NMSEL=0;
for j=1:n
        NMSEU=(QthetaU(j)-thetaU(j)).^2+NMSEU;
        NMSEL=(QthetaL(j)-thetaL(j)).^2+NMSEL;
end

RNMSEAR1U(i) = NMSEU/n;
RNMSEAR1L(i) = NMSEL/n;
[ kappa12] = fmincon('otherMSJC',theta0,[],[],[],[],lower,upper,[],options,[u,v],kappa);
[LL12(i)] =  otherMSJC(kappa12,[u,v],kappa);


% fprintf(' \n  %1f\t%1f\t%1f\t%1f\t%1f',i,RNMSEAR1U(i) ,RNMSEAR1L(i),LL11(i) ,LL12(i));
% fprintf(' \n 1 %1f	',i);

end
SAVE1(:,1) = RNMSEAR1U(:,1);
SAVE1(:,2) = RNMSEAR1L(:,1);
SAVE1(:,3) = LL11(:,1);
SAVE1(:,4) = LL12(:,1);
csvwrite('SAVE9.csv',SAVE1)


theta0 = [0.2;0;0;0.2;0;0;-0.3;-0.3;0.9;0.9];
RNMSEsineL= zeros(1000,1);
RNMSEsineU= zeros(1000,1);
LL13= zeros(1000,1);
LL14= zeros(1000,1);
SAVE2=zeros(1000,4);
kappa = [0.5,0.9];
parfor i=1:1000
n=1000;%nnnnnnnnnnnnnnnnnnnnnnnnnnn=1000 2000 5000
[u,v,thetaU,thetaL]= Simsine(n);
[ kappa11] = fmincon('MSJC',theta0,[],[],[],[],lower,upper,[],options,[u,v],kappa);
[LL13(i), QthetaU,QthetaL] = Copy_of_MSJC(kappa11,[u,v],kappa);
NMSEU=0;NMSEL=0;
for j=1:n
        NMSEU=(QthetaU(j)-thetaU(j)).^2+NMSEU;
        NMSEL=(QthetaL(j)-thetaL(j)).^2+NMSEL;
end
RNMSEsineU(i) = NMSEU/n;
RNMSEsineL(i) = NMSEL/n;

[ kappa12] = fmincon('otherMSJC',theta0,[],[],[],[],lower,upper,[],options,[u,v],kappa);
[LL14(i)] =  otherMSJC(kappa12,[u,v],kappa);

% 
% fprintf(' \n  %1f\t%1f\t%1f\t%1f\t%1f',i,RNMSEsineU(i) ,RNMSEsineL(i),LL13(i) ,LL14(i));
% fprintf(' \n 2 %1f	',i);

end
SAVE2(:,1) = RNMSEsineU(:,1);
SAVE2(:,2) = RNMSEsineL(:,1);
SAVE2(:,3) = LL13(:,1);
SAVE2(:,4) = LL14(:,1);
csvwrite('SAVE10.csv',SAVE2)
% 

theta0 = [0.2;0;0;0.2;0;0;-0.3;-0.3;0.9;0.9];
RNMSEstepL= zeros(1000,1);
RNMSEstepU= zeros(1000,1);
SAVE3=zeros(1000,4);
LL15= zeros(1000,1);
LL16= zeros(1000,1);
kappa = [0.5,0.9];

parfor i=1:1000
n=1000;%nnnnnnnnnnnnnnnnnnnnnnnnnnn=1000 2000 5000
[u,v,thetaU,thetaL]= Simstep(n);
[ kappa11] = fmincon('MSJC',theta0,[],[],[],[],lower,upper,[],options,[u,v],kappa);
[LL15(i), QthetaU,QthetaL] = Copy_of_MSJC(kappa11,[u,v],kappa);
NMSEU=0;NMSEL=0;
for j=1:n
        NMSEU=(QthetaU(j)-thetaU(j)).^2+NMSEU;
        NMSEL=(QthetaL(j)-thetaL(j)).^2+NMSEL;
end
RNMSEstepU(i) = NMSEU/n;
RNMSEstepL(i) = NMSEL/n;

[ kappa12] = fmincon('otherMSJC',theta0,[],[],[],[],lower,upper,[],options,[u,v],kappa);
[LL16(i)] =  otherMSJC(kappa12,[u,v],kappa);
NMSEU=0;NMSEL=0;


% fprintf(' \n  %1f\t%1f\t%1f\t%1f\t%1f',i,RNMSEstepU(i) ,RNMSEstepL(i),LL15(i) ,LL16(i));
% fprintf(' \n 3  %1f	',i);
end

SAVE3(:,1) = RNMSEstepU(:,1);
SAVE3(:,2) = RNMSEstepL(:,1);
SAVE3(:,3) = LL15(:,1);
SAVE3(:,4) = LL16(:,1);
csvwrite('SAVE11.csv',SAVE3)


theta0 = [0.2;0;0;0.2;0;0;-0.3;-0.3;0.9;0.9];
RNMSEotherL= zeros(1000,1);
RNMSEotherU= zeros(1000,1);
otherTAL1= zeros(1000,1);
otherTAL2= zeros(1000,1);
otherTAU1= zeros(1000,1);
otherTAU2= zeros(1000,1);
otherLLLL= zeros(1000,1);
LL17= zeros(1000,1);
LL18= zeros(1000,1);
SAVE4=zeros(1000,6);
otherRNU= zeros(1000,1);
otherRNL= zeros(1000,1);
kappa = [0.5,0.9];
parfor i=1:1000
n=1000;%nnnnnnnnnnnnnnnnnnnnnnnnnnn=1000 2000 5000
[u,v,thetaU,thetaL,state]= Simother(n);
[ kappa11] = fmincon('MSJC',theta0,[],[],[],[],lower,upper,[],options,[u,v],kappa);
[LL17(i), QthetaU,QthetaL] = Copy_of_MSJC(kappa11,[u,v],kappa);
NMSEU=0;NMSEL=0;
for j=1:n
        NMSEU=(QthetaU(j)-thetaU(j)).^2+NMSEU;
        NMSEL=(QthetaL(j)-thetaL(j)).^2+NMSEL;
end
RNMSEotherL(i) = NMSEU/n;
RNMSEotherU(i) = NMSEL/n;

[ kappa12] = fmincon('otherMSJC',theta0,[],[],[],[],lower,upper,[],options,[u,v],kappa);
[LL18(i), otherTAL1,otherTAL2,otherTAU1,otherTAU2] =  CopyotherMSJC(kappa12,[u,v],kappa);
NMSEU=0;NMSEL=0;

if kappa12(1)>kappa12(7)
    
    otherLLLL=otherTAL1;
    otherTAL1=otherTAL2;
    otherTAL2=otherLLLL;
end
if kappa12(4)>kappa12(8)
    
    otherLLLL=otherTAU1;
    otherTAU1=otherTAU2;
    otherTAU2=otherLLLL;
end
otherQthetaU = zeros(n,1);
otherQthetaL = zeros(n,1);
for j=1:n
    if state(j)>0.5
        otherQthetaU(j) = otherTAU2(j);
        otherQthetaL(j) = otherTAL2(j);
        NMSEU=(otherQthetaU(j)-thetaU(j)).^2+NMSEU;
        NMSEL=(otherQthetaL(j)-thetaL(j)).^2+NMSEL;
    else
        otherQthetaU(j) = otherTAU1(j);
        otherQthetaL(j) = otherTAL1(j);
        NMSEU=(otherQthetaU(j)-thetaU(j)).^2+NMSEU;
        NMSEL=(otherQthetaL(j)-thetaL(j)).^2+NMSEL;
   
    end
end

otherRNU(i) = NMSEU/n;
otherRNL(i) = NMSEL/n;
% fprintf(' \n  %1f\t%1f\t%1f\t%1f\t%1f\t%1f\t%1f',i,RNMSEotherL(i) ,RNMSEotherU(i),otherRNU(i) ,otherRNL(i),LL17(i),LL18(i));
% fprintf(' \n 3  %1f	',i);
end

SAVE4(:,1) = RNMSEotherL(:,1);
SAVE4(:,2) = RNMSEotherU(:,1);
SAVE4(:,3) = otherRNU(:,1);
SAVE4(:,4) = otherRNL(:,1);
SAVE4(:,5) = LL17(:,1);
SAVE4(:,6) = LL18(:,1);

csvwrite('SAVE12.csv',SAVE4)



RNMSEAR1U= zeros(1000,1);
RNMSEAR1L= zeros(1000,1);
LL11= zeros(1000,1);
LL12= zeros(1000,1);
theta0 = [0.2;0.1;0;0.2;0;0;0.3;0.3;0.9;0.9];
kappa = [0.2,0.2];
SAVE1=zeros(1000,4);

% 
parfor i=1:1000
n=2000;%nnnnnnnnnnnnnnnnnnnnnnnnnnn=1000 2000 5000
[u,v,thetaU,thetaL]= SimAR1(n);
[ kappa11] = fmincon('MSJC',theta0,[],[],[],[],lower,upper,[],options,[u,v],kappa);

[LL11(i), QthetaU,QthetaL] = Copy_of_MSJC(kappa11,[u,v],kappa);
NMSEU=0;NMSEL=0;
for j=1:n
        NMSEU=(QthetaU(j)-thetaU(j)).^2+NMSEU;
        NMSEL=(QthetaL(j)-thetaL(j)).^2+NMSEL;
end

RNMSEAR1U(i) = NMSEU/n;
RNMSEAR1L(i) = NMSEL/n;
[ kappa12] = fmincon('otherMSJC',theta0,[],[],[],[],lower,upper,[],options,[u,v],kappa);
[LL12(i)] =  otherMSJC(kappa12,[u,v],kappa);


% fprintf(' \n  %1f\t%1f\t%1f\t%1f\t%1f',i,RNMSEAR1U(i) ,RNMSEAR1L(i),LL11(i) ,LL12(i));
% fprintf(' \n 1 %1f	',i);

end
SAVE1(:,1) = RNMSEAR1U(:,1);
SAVE1(:,2) = RNMSEAR1L(:,1);
SAVE1(:,3) = LL11(:,1);
SAVE1(:,4) = LL12(:,1);
csvwrite('SAVE5.csv',SAVE1)


theta0 = [0.2;0;0;0.2;0;0;-0.3;-0.3;0.9;0.9];
RNMSEsineL= zeros(1000,1);
RNMSEsineU= zeros(1000,1);
LL13= zeros(1000,1);
LL14= zeros(1000,1);
SAVE2=zeros(1000,4);
kappa = [0.5,0.9];
parfor i=1:1000
n=2000;%nnnnnnnnnnnnnnnnnnnnnnnnnnn=1000 2000 5000
[u,v,thetaU,thetaL]= Simsine(n);
[ kappa11] = fmincon('MSJC',theta0,[],[],[],[],lower,upper,[],options,[u,v],kappa);
[LL13(i), QthetaU,QthetaL] = Copy_of_MSJC(kappa11,[u,v],kappa);
NMSEU=0;NMSEL=0;
for j=1:n
        NMSEU=(QthetaU(j)-thetaU(j)).^2+NMSEU;
        NMSEL=(QthetaL(j)-thetaL(j)).^2+NMSEL;
end
RNMSEsineU(i) = NMSEU/n;
RNMSEsineL(i) = NMSEL/n;

[ kappa12] = fmincon('otherMSJC',theta0,[],[],[],[],lower,upper,[],options,[u,v],kappa);
[LL14(i)] =  otherMSJC(kappa12,[u,v],kappa);

% 
% fprintf(' \n  %1f\t%1f\t%1f\t%1f\t%1f',i,RNMSEsineU(i) ,RNMSEsineL(i),LL13(i) ,LL14(i));
% fprintf(' \n 2 %1f	',i);

end
SAVE2(:,1) = RNMSEsineU(:,1);
SAVE2(:,2) = RNMSEsineL(:,1);
SAVE2(:,3) = LL13(:,1);
SAVE2(:,4) = LL14(:,1);
csvwrite('SAVE6.csv',SAVE2)
% 

theta0 = [0.2;0;0;0.2;0;0;-0.3;-0.3;0.9;0.9];
RNMSEstepL= zeros(1000,1);
RNMSEstepU= zeros(1000,1);
SAVE3=zeros(1000,4);
LL15= zeros(1000,1);
LL16= zeros(1000,1);
kappa = [0.5,0.9];

parfor i=1:1000
n=2000;%nnnnnnnnnnnnnnnnnnnnnnnnnnn=1000 2000 5000
[u,v,thetaU,thetaL]= Simstep(n);
[ kappa11] = fmincon('MSJC',theta0,[],[],[],[],lower,upper,[],options,[u,v],kappa);
[LL15(i), QthetaU,QthetaL] = Copy_of_MSJC(kappa11,[u,v],kappa);
NMSEU=0;NMSEL=0;
for j=1:n
        NMSEU=(QthetaU(j)-thetaU(j)).^2+NMSEU;
        NMSEL=(QthetaL(j)-thetaL(j)).^2+NMSEL;
end
RNMSEstepU(i) = NMSEU/n;
RNMSEstepL(i) = NMSEL/n;

[ kappa12] = fmincon('otherMSJC',theta0,[],[],[],[],lower,upper,[],options,[u,v],kappa);
[LL16(i)] =  otherMSJC(kappa12,[u,v],kappa);
NMSEU=0;NMSEL=0;


% fprintf(' \n  %1f\t%1f\t%1f\t%1f\t%1f',i,RNMSEstepU(i) ,RNMSEstepL(i),LL15(i) ,LL16(i));
% fprintf(' \n 3  %1f	',i);
end

SAVE3(:,1) = RNMSEstepU(:,1);
SAVE3(:,2) = RNMSEstepL(:,1);
SAVE3(:,3) = LL15(:,1);
SAVE3(:,4) = LL16(:,1);
csvwrite('SAVE7.csv',SAVE3)


theta0 = [0.2;0;0;0.2;0;0;-0.3;-0.3;0.9;0.9];
RNMSEotherL= zeros(1000,1);
RNMSEotherU= zeros(1000,1);
otherTAL1= zeros(1000,1);
otherTAL2= zeros(1000,1);
otherTAU1= zeros(1000,1);
otherTAU2= zeros(1000,1);
otherLLLL= zeros(1000,1);
LL17= zeros(1000,1);
LL18= zeros(1000,1);
SAVE4=zeros(1000,6);
otherRNU= zeros(1000,1);
otherRNL= zeros(1000,1);
kappa = [0.5,0.9];
parfor i=1:1000
n=2000;%nnnnnnnnnnnnnnnnnnnnnnnnnnn=1000 2000 5000
[u,v,thetaU,thetaL,state]= Simother(n);
[ kappa11] = fmincon('MSJC',theta0,[],[],[],[],lower,upper,[],options,[u,v],kappa);
[LL17(i), QthetaU,QthetaL] = Copy_of_MSJC(kappa11,[u,v],kappa);
NMSEU=0;NMSEL=0;
for j=1:n
        NMSEU=(QthetaU(j)-thetaU(j)).^2+NMSEU;
        NMSEL=(QthetaL(j)-thetaL(j)).^2+NMSEL;
end
RNMSEotherL(i) = NMSEU/n;
RNMSEotherU(i) = NMSEL/n;

[ kappa12] = fmincon('otherMSJC',theta0,[],[],[],[],lower,upper,[],options,[u,v],kappa);
[LL18(i), otherTAL1,otherTAL2,otherTAU1,otherTAU2] =  CopyotherMSJC(kappa12,[u,v],kappa);
NMSEU=0;NMSEL=0;

if kappa12(1)>kappa12(7)
    
    otherLLLL=otherTAL1;
    otherTAL1=otherTAL2;
    otherTAL2=otherLLLL;
end
if kappa12(4)>kappa12(8)
    
    otherLLLL=otherTAU1;
    otherTAU1=otherTAU2;
    otherTAU2=otherLLLL;
end
otherQthetaU = zeros(n,1);
otherQthetaL = zeros(n,1);
for j=1:n
    if state(j)>0.5
        otherQthetaU(j) = otherTAU2(j);
        otherQthetaL(j) = otherTAL2(j);
        NMSEU=(otherQthetaU(j)-thetaU(j)).^2+NMSEU;
        NMSEL=(otherQthetaL(j)-thetaL(j)).^2+NMSEL;
    else
        otherQthetaU(j) = otherTAU1(j);
        otherQthetaL(j) = otherTAL1(j);
        NMSEU=(otherQthetaU(j)-thetaU(j)).^2+NMSEU;
        NMSEL=(otherQthetaL(j)-thetaL(j)).^2+NMSEL;
   
    end
end

otherRNU(i) = NMSEU/n;
otherRNL(i) = NMSEL/n;
% fprintf(' \n  %1f\t%1f\t%1f\t%1f\t%1f\t%1f\t%1f',i,RNMSEotherL(i) ,RNMSEotherU(i),otherRNU(i) ,otherRNL(i),LL17(i),LL18(i));
% fprintf(' \n 3  %1f	',i);
end

SAVE4(:,1) = RNMSEotherL(:,1);
SAVE4(:,2) = RNMSEotherU(:,1);
SAVE4(:,3) = otherRNU(:,1);
SAVE4(:,4) = otherRNL(:,1);
SAVE4(:,5) = LL17(:,1);
SAVE4(:,6) = LL18(:,1);

csvwrite('SAVE8.csv',SAVE4)



fprintf(1,'\n finish  ');
% datestr(now,31)
% figure(i),plot((1:n)',otherQthetaU,(1:n)',thetaU,'m--'),legend('simulation','real'),title('otherU');
% figure(i+4),plot((1:n)',QthetaU,(1:n)',thetaU,'m--'),legend('simulation','real'),title('weU');
% figure(i+8),plot((1:n)',otherQthetaL,(1:n)',thetaL,'m--'),legend('simulation','real'),title('otherL');
% figure(i+12),plot((1:n)',QthetaL,(1:n)',thetaL,'m--'),legend('simulation','real'),title('weL');

