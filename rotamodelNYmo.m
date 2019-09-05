%Age-structured model for rotavirus -- 1 strain w/ 3 stages of
%susceptibility, 3 stages of infectiousness (32 age classes)
clear global

global B wm wi1 wi2 u um beta b1 phi d1 d2 rr1 rr2 ri2 ri3 al v2 boost wA;  
%By setting our parameters as "global", we allow MATLAB to call them both
%from the Desktop environment and from within functions, which we'll need
%later on. 

age=[0:1/12:11/12 1:4 7.5 15 30 50 70];
al=length(age);
avgage=age;
agep=[(1/960)*ones(1,12) (1/80)*ones(1,4) 1/16 1/8 1/4 1/4 1/4];

Nstate_NY=17.99*1000000; %State population in 1990
Bstate_NY=[.0165; .0161; .0158; .0154; .0151; .0146; .0142; .0138; .0138; .0135; .0136; .0133; .0131; .0132; .0130; .0128; .0130; .0131; .0128; .0127; .0125; .0123; .0122; .0120; .0121]; %Annual crude birth rate, 1990-2014      
N0_nyc=7000000; %NYC pouplation in 1980

%Birth rate by month for NY state 
BstNY=[ones(12*17.5,1)*.0175; ones(60,1)*mean(Bstate_NY); ones(12,1)*Bstate_NY(1,1); ones(12,1)*Bstate_NY(2,1); ones(12,1)*Bstate_NY(3,1); 
    ones(12,1)*Bstate_NY(4,1); ones(12,1)*Bstate_NY(5,1); ones(12,1)*Bstate_NY(6,1); ones(12,1)*Bstate_NY(7,1); ones(12,1)*Bstate_NY(8,1); 
    ones(12,1)*Bstate_NY(9,1); ones(12,1)*Bstate_NY(10,1); ones(12,1)*Bstate_NY(11,1); ones(12,1)*Bstate_NY(12,1); 
    ones(12,1)*Bstate_NY(13,1); ones(12,1)*Bstate_NY(14,1); ones(12,1)*Bstate_NY(15,1); ones(12,1)*Bstate_NY(16,1); 
    ones(12,1)*Bstate_NY(17,1); ones(12,1)*Bstate_NY(18,1); ones(12,1)*Bstate_NY(19,1); ones(12,1)*Bstate_NY(20,1);
    ones(12,1)*Bstate_NY(21,1); ones(12,1)*Bstate_NY(22,1); ones(12,1)*Bstate_NY(23,1);  ones(12,1)*Bstate_NY(24,1);
    ones(12*10,1)*Bstate_NY(25,1)]; %ones(52*2,1)*Bstate_NY(25,1)]; 

tmax=length(BstNY); %length of simulation
t0=12*33.5; %length of "burn-in" period

N=N0_nyc*agep; %Nstate_NY*agep; %population size (in each age group)
B=[BstNY zeros(tmax,al-1)]; %number of new births per week into each age class
u=[ones(1,12) 1/12*ones(1,4) 1/(12*5) 1/120 1/240 1/240 1/240]; %rate for aging out of each age class (per person per week)
um=0; %mortality rate in individuals < Amax=70 yrs 

%%% Best-fit parameters estimated from age-specific hospitalization data (see Pitzer et al, Science, 2009)
parNY=[24.1842; .0492; .6776; .0369; 1.5704; 2.0822; 1.4958]; 

%%% CONTACT RATES %%%
cNY=[parNY(5,1)*ones(al,12) parNY(6,1)*ones(al,1) parNY(7,1)*ones(al,1) ones(al,al-14)]; %Age-related acquisition for <3 yr olds

%%% INFECTION PARAMETERS %%%
dur=1/4.3; %duration of infectiousness (1 week = 0.23 months)
d1=1/dur; %rate of recovery from primary infection (per month)
d2=2*d1; %rate of recovery from subsequent infection (per month)
rr1=0.62;  %relative risk of second infection 
rr2=0.35; %relative risk of third infection 
ri2=0.5;   %relative infectiousness of secondary infection
ri3=0.1;   %relative infectiousness of asymptomatic infection 
wm=1/3; %rate of waning maternal immunity (avg duration of immunity = 3 months)
wi1=1/9; %rate of waning immunity following primary infection (avg duration of immunity = 9 months)
wi2=1/12; %rate of waning immunity following 2nd infection (avg duration of immunity = 12 months)
wA=0; %rate of waning of immunity to symptomatic infections (S2->S1)

%%% INPUT PARAMETERS %%%
ptrans=parNY(1)/dur; %probability of transmission given contact
beta=ptrans/sum(N)*cNY; %transmission matrix 

b1=parNY(2); %seasonal forcing 
phi=parNY(3); %seasonal offset

h=parNY(4); %proportion of severe diarrhea cases hospitalized
hosp1=0.13*h*ones(1,al); %proportion of primary infections with severe diarrhea who are hospitalized
hosp2=0.03*h*ones(1,al); %proportion of secondary infections with severe diarrhea who are hospitalized
hosp3=zeros(1,al); %proportion of subsequent infections with severe diarrhea who are hospitalized
delta1=0.45*ones(1,al); %proportion of primary infections that are symptomatic (any diarrhea)
delta2=0.25*ones(1,al); %proportion of secondary infections that are symptomatic (any diarrhea)
delta3=zeros(1,al); %rate of detection of subsequent infections

R0=max(eig(dur*beta.*(ones(al,1)*N)));

%%% INITIALIZE OUTCOME VARIABLES %%%
H_nyc1=zeros(tmax-t0,4);
H_nyc2=zeros(tmax-t0,4);
H_nyc3=zeros(tmax-t0,4);
D_nyc1=zeros(tmax-t0,4);
D_nyc2=zeros(tmax-t0,4);
D_nyc3=zeros(tmax-t0,4);
SSE_nyc=zeros(16,1);

%%% VACCINATION PARAMETERS %%%
tvacc=48; %month prior to vaccine introduction (=January 2006)
vdate=datevec(datenum([2006 1 15 0 0 0]):30:datenum([2015 12 31 0 0 0]));
vdateNY=[2006*ones(12,1) (1:12)' 15*ones(12,1) zeros(12,3); 2007*ones(12,1) (1:12)' 15*ones(12,1) zeros(12,3); 2008*ones(12,1) (1:12)' 15*ones(12,1) zeros(12,3); 2009*ones(12,1) (1:12)' 15*ones(12,1) zeros(12,3); 2010*ones(12,1) (1:12)' 15*ones(12,1) zeros(12,3); 2011*ones(12,1) (1:12)' 15*ones(12,1) zeros(12,3); 2012*ones(12,1) (1:12)' 15*ones(12,1) zeros(12,3); 2013*ones(12,1) (1:12)' 15*ones(12,1) zeros(12,3); 2014*ones(12,1) (1:12)' 15*ones(12,1) zeros(12,3); 2015*ones(12,1) (1:12)' 15*ones(12,1) zeros(12,3)];
vcovNYmo=csvread('monthly_vcov_nyc.csv');

for v=19 %1:21 %1:5  
ve1=(69+v)/100; %(65+5*v)/100; %relative effectiveness of 1-dose coverage estimates (compared to full coverage efficacy) 
ve2=.95; %relative effectiveness of second dose of vaccine (based on seroconversion data)
v2=zeros(tmax,al); %initialize vaccination rate across all ages
boost=zeros(tmax,al); %boosting effect of second dose of vaccine
for scenario=1:3
if scenario==1
v2(:,3)=[zeros(t0+tvacc,1); vcovNYmo(1:120,1); vcovNYmo(121,1)*ones(24,1); .932*ones(tmax-t0-tvacc-120-24,1)]*ve1; %proportion vaccinated at 2 months of age
boost(:,6)=[zeros(t0+tvacc,1); vcovNYmo(1:120,2); vcovNYmo(121,2)*ones(24,1); .7*ones(tmax-t0-tvacc-120-24,1)]*ve2; %proportion vaccinated with second dose at 4 months of age
elseif scenario==2    
v2(:,3)=[zeros(t0+tvacc,1); vcovNYmo(1:120,1); vcovNYmo(121,1)*ones(24,1); .932*ones(tmax-t0-tvacc-120-24,1)]*ve1; %proportion vaccinated at 2 months of age
boost(:,6)=[zeros(t0+tvacc,1); vcovNYmo(1:120,2); vcovNYmo(121,2)*ones(24,1); .932*ones(tmax-t0-tvacc-120-24,1)]*ve2; %proportion vaccinated with second dose at 4 months of age
elseif scenario==3    
v2(:,3)=[zeros(t0+tvacc,1); vcovNYmo(1:120,1); vcovNYmo(121,1)*ones(36,1); .932*ones(tmax-t0-tvacc-120-36,1)]*ve1; %proportion vaccinated at 2 months of age
boost(:,6)=[zeros(t0+tvacc,1); vcovNYmo(1:120,2); vcovNYmo(121,2)*ones(36,1); .932*ones(tmax-t0-tvacc-120-36,1)]*ve2; %proportion vaccinated with second dose at 4 months of age
end
vcov1(:,scenario)=100*v2(t0+1:end,3)/ve1;
vcov2(:,scenario)=100*boost(t0+1:end,6)/ve2;

%Initialize vector to keep track of the number of people in each state
St0=[];
St0(1:al,1)=[N(1) zeros(1,al-1)]; %Maternal immunity
St0(al+1:2*al,1)=[0 N(2:al)-ones(1,al-1)]; %Susceptible_0
St0(2*al+1:3*al,1)=[0 ones(1,al-1)]; %Infectious_1 (primary) 
St0(3*al+1:4*al,1)=zeros(1,al); %Recovered_1
St0(4*al+1:5*al,1)=zeros(1,al); %Susceptible_1
St0(5*al+1:6*al,1)=zeros(1,al); %Infectious_2 (2nd time)
St0(6*al+1:7*al,1)=zeros(1,al); %Recovered_2
St0(7*al+1:8*al,1)=zeros(1,al); %Susceptible-Resistant
St0(8*al+1:9*al,1)=zeros(1,al); %Asymptomatic Infectious_3 (subsequent)
St0(9*al+1:10*al,1)=zeros(1,al); %Temp Resistant

clear St lambda H %clear outcome variables which may be in memory

options=odeset('NonNegative',1:length(St0)); %force solutions to differential equations to be non-negative
[time St]=ode45('rasisAM',1:tmax,St0,options); %solve differential equations defined in 'rasisAW'

time(1:t0,:)=[]; %delete output from from burn-in period
St(1:t0,:)=[]; 

%Initialize outcome variables
lambda=zeros(tmax-t0,al); H=zeros(tmax-t0,al); D=zeros(tmax-t0,al);

for t=1:length(St) %calculate force of infection at time t
    lambda(t,:)=(1+b1*cos(2*pi*(time(t)-phi*12)/12))*((St(t,2*al+1:3*al)+ri2*St(t,5*al+1:6*al)+ri3*St(t,8*al+1:9*al))*beta);
end

for i=1:al %calculate number of rotavirus hospitalizations (H) and diarrhea cases (D) in each age group across time
    H(:,i)=hosp1(i)*St(:,al+i).*lambda(:,i)+hosp2(i)*rr1*St(:,4*al+i).*lambda(:,i)+hosp3(i)*rr2*St(:,7*al+i).*lambda(:,i);
    D(:,i)=delta1(i)*St(:,al+i).*lambda(:,i)+delta2(i)*rr1*St(:,4*al+i).*lambda(:,i)+delta3(i)*rr2*St(:,7*al+i).*lambda(:,i);
end

pop=St(:,1:al);
for j=1:9
    pop=pop+St(:,j*al+1:(j+1)*al);
end

H_nyc1(:,v)=sum(H(1:end,1:13),2); %./sum(pop(1:end,1:13),2);
H_nyc2(:,v)=sum(H(1:end,14:16),2); %./sum(pop(1:end,14:16),2);
H_nyc3(:,v)=sum(H(1:end,17),2); %./sum(pop(1:end,17:21),2);
D_nyc1(:,v)=sum(D(1:end,1:13),2)./sum(pop(1:end,1:13),2);
D_nyc2(:,v)=sum(D(1:end,14:16),2)./sum(pop(1:end,14:16),2);
D_nyc3(:,v)=sum(D(1:end,17),2)./sum(pop(1:end,17:21),2);

if v==19
Hnyc1_vacc(:,scenario)=H_nyc1(:,v);    
Hnyc2_vacc(:,scenario)=H_nyc2(:,v);    
Hnyc3_vacc(:,scenario)=H_nyc3(:,v);    

Hnyc_vacc=Hnyc1_vacc+Hnyc2_vacc+Hnyc3_vacc;

agedist_prevacc=[sum(H_nyc1(7:54,v)) sum(H_nyc2(7:54,v)) sum(H_nyc3(7:54,v))]/sum(Hnyc_vacc(7:54,1));
for i=1:14
    cases_byyr(i,:)=[sum(H_nyc1(12*(i-1)+7:12*i+6,v)) sum(H_nyc2(12*(i-1)+7:12*i+6,v)) sum(H_nyc3(12*(i-1)+7:12*i+6,v)) sum(Hnyc_vacc(12*(i-1)+7:12*i+6,1))];
    agedist_byyr(i,:)=[sum(H_nyc1(12*(i-1)+7:12*i+6,v)) sum(H_nyc2(12*(i-1)+7:12*i+6,v)) sum(H_nyc3(12*(i-1)+7:12*i+6,v))]/sum(Hnyc_vacc(12*(i-1)+7:12*i+6,1));
end
end


rotaNYC_hosp=csvread('monthly_RVhosp_nyc.csv');
rotaNYC_lab=csvread('monthly_RVlab_nyc.csv');
labmo_nyc=[2008*ones(10,1) (3:12)' ones(10,1); 2009*ones(12,1) (1:12)' ones(12,1); 2010*ones(12,1) (1:12)' ones(12,1); 2011*ones(12,1) (1:12)' ones(12,1); 2012*ones(12,1) (1:12)' ones(12,1); 2013*ones(12,1) (1:12)' ones(12,1); 2014*ones(12,1) (1:12)' ones(12,1); 2015*ones(12,1) (1:12)' ones(12,1); 2016*ones(12,1) (1:12)' ones(12,1)];
hospmo_nyc=[2002*ones(12,1) (1:12)' ones(12,1); 2003*ones(12,1) (1:12)' ones(12,1); 2004*ones(12,1) (1:12)' ones(12,1); 2005*ones(12,1) (1:12)' ones(12,1); 2006*ones(12,1) (1:12)' ones(12,1); 2007*ones(12,1) (1:12)' ones(12,1); 2008*ones(2,1) (1:2)' ones(2,1); labmo_nyc];

[corr_hosplab,corr_pvalue]=corr(sum(rotaNYC_hosp(75:end,:),2),sum(rotaNYC_lab(1:106,:),2));

agedist_hosp_prevacc=[sum(rotaNYC_hosp(7:54,2)) sum(rotaNYC_hosp(7:54,3)) sum(sum(rotaNYC_hosp(7:54,4:5)))]/sum(rotaNYC_hosp(7:54,1));
for i=1:14
    rotaNYC_hosp_byyr(i,:)=[sum(rotaNYC_hosp(12*(i-1)+7:12*i+6,2)) sum(rotaNYC_hosp(12*(i-1)+7:12*i+6,3)) sum(sum(rotaNYC_hosp(12*(i-1)+7:12*i+6,4:5))) sum(rotaNYC_hosp(12*(i-1)+7:12*i+6,1))];
    agedist_hosp_byyr(i,:)=[sum(rotaNYC_hosp(12*(i-1)+7:12*i+6,2)) sum(rotaNYC_hosp(12*(i-1)+7:12*i+6,3)) sum(sum(rotaNYC_hosp(12*(i-1)+7:12*i+6,4:5)))]/sum(rotaNYC_hosp(12*(i-1)+7:12*i+6,1));
end

rotaNYC_lab_byyr(1,:)=[sum(rotaNYC_lab(1:4,1)) sum(rotaNYC_lab(1:4,2)) sum(sum(rotaNYC_lab(1:4,3:4))) sum(rotaNYC_lab(1:4,5))];
agedist_lab_byyr(1,:)=[sum(rotaNYC_lab(1:4,1)) sum(rotaNYC_lab(1:4,2)) sum(sum(rotaNYC_lab(1:4,3:4)))]/sum(rotaNYC_lab(1:4,5));
for i=1:8
    rotaNYC_lab_byyr(i+1,:)=[sum(rotaNYC_lab(12*(i-1)+5:12*i+4,1)) sum(rotaNYC_lab(12*(i-1)+5:12*i+4,2)) sum(sum(rotaNYC_lab(12*(i-1)+5:12*i+4,3:4))) sum(rotaNYC_lab(12*(i-1)+5:12*i+4,5))];
    agedist_lab_byyr(i+1,:)=[sum(rotaNYC_lab(12*(i-1)+5:12*i+4,1)) sum(rotaNYC_lab(12*(i-1)+5:12*i+4,2)) sum(sum(rotaNYC_lab(12*(i-1)+5:12*i+4,3:4)))]/sum(rotaNYC_lab(12*(i-1)+5:12*i+4,5));
end

SSE_nyc_lab(v,1)=sum((rotaNYC_lab(1:106,1)-H_nyc1(75:180,v)).^2 + (rotaNYC_lab(1:106,2)-H_nyc2(75:180,v)).^2 + (sum(rotaNYC_lab(1:106,3:4),2)-H_nyc3(75:180,v)).^2);
SSE_nyc_hosp(v,1)=sum((rotaNYC_hosp(73:180,2)-H_nyc1(73:180,v)).^2 + (rotaNYC_hosp(73:180,3)-H_nyc2(73:180,v)).^2 + (sum(rotaNYC_hosp(73:180,4:5),2)-H_nyc3(73:180,v)).^2);
SSE_nyc=SSE_nyc_lab+SSE_nyc_hosp;

SSTot_lab=sum((rotaNYC_lab(1:106,1)-mean(rotaNYC_lab(1:106,1))).^2 + (rotaNYC_lab(1:106,2)-mean(rotaNYC_lab(1:106,2))).^2 + (sum(rotaNYC_lab(1:106,3:4),2)-mean(sum(rotaNYC_lab(1:106,3:4),2))).^2);
SSTot_hosp=sum((rotaNYC_hosp(73:180,2)-mean(rotaNYC_hosp(73:180,2))).^2 + (rotaNYC_hosp(73:180,3)-mean(rotaNYC_hosp(73:180,3))).^2 + (sum(rotaNYC_hosp(73:180,4:5),2)-mean(sum(rotaNYC_hosp(73:180,4:5),2))).^2);

Rsquared_lab=1-SSE_nyc_lab(v,1)/SSTot_lab;
Rsquared_hosp=1-SSE_nyc_hosp(v,1)/SSTot_hosp;

end
end

%%
%v=19;

figure
subplot(4,1,1); hold on
plot(H_nyc1(:,v),'b'); 
plot(75:182,rotaNYC_lab(:,1),'Color',[.9 .7 .2])
plot(1:180,rotaNYC_hosp(:,2),'Color','r')
set(gca,'XLim',[0 276],'XTick',0:24:264,'XTickLabel',2002:2:2024,'FontSize',9);
title('<2 year of age')
%legend('75% relative efficacy','80% relative efficacy','85% relative efficacy','90% relative efficacy')

subplot(4,1,2); hold on
plot(H_nyc2(:,v),'b'); 
plot(1:180,rotaNYC_hosp(:,3),'Color','r')
plot(75:182,rotaNYC_lab(:,2),'Color',[.9 .7 .2])
ylabel('Predicted incidence of severe RVGE (per month)')
set(gca,'XLim',[0 276],'XTick',0:24:264,'XTickLabel',2002:2:2024,'FontSize',9);
title('2-4 years of age')

subplot(4,1,3); hold on
plot(H_nyc3(:,v),'b'); 
plot(1:180,rotaNYC_hosp(:,4),'Color','r')
plot(75:182,sum(rotaNYC_lab(:,3:4),2),'Color',[.9 .7 .2])
set(gca,'XLim',[0 276],'XTick',0:24:264,'XTickLabel',2002:2:2024,'FontSize',9);
title('5+ years of age')

subplot(4,1,4)
hold on; box on
plot(100*v2(t0+1:end,3)/ve1,'b'); 
plot(100*boost(t0+1:end,6)/ve2,'Color',[0 .7 0]); 
set(gca,'XLim',[0 276],'XTick',0:24:264,'XTickLabel',2002:2:2024,'FontSize',9);
title('Vaccine coverage')
ylabel('% vaccinated')
legend('\geq1 dose','full schedule')

%%
%legend('Model prediction','Observed (lab reports)')
%legend('75% relative efficacy','80% relative efficacy','85% relative efficacy','90% relative efficacy')
%gtext('\bfA')
%gtext('\bfB')
%gtext('\bfC')
%gtext('\bfD')

%%
%figure
%hold on
%plot(datenum(hospmo_nyc),H_nyc1(1:180,:)+H_nyc2(1:180,:)+H_nyc3(1:180,:))
%plot(datenum(hospmo_nyc),rotaNYC_hosp(:,1),'k')
%plot(datenum(labmo_nyc),rotaNYC_lab(1:106,5),'Color',[0 .7 0])
%datetick('x','mmm-yy')
%legend('70% relative efficacy','80% relative efficacy','90% relative efficacy','100% relative efficacy')

%%
%figure 
%subplot(2,1,1); plot(sum(H,2)); set(gca,'XLim',[0 180],'XTick',12:24:180,'XTickLabel',2003:2:2017,'FontSize',9);
%subplot(2,1,2); plot(sum(pop,2)); set(gca,'XLim',[0 180],'XTick',12:24:180,'XTickLabel',2003:2:2017,'FontSize',9);

%%
%figure
%hold on
%plot(Hnyc_vacc)
%plot([tvacc+144 tvacc+144],[0 350],'--b')
%plot([tvacc+156 tvacc+156],[0 350],'--r')
%plot(Hnyc_vacc(:,1),'b')
%set(gca,'XLim',[0 276],'XTick',0:24:264,'XTickLabel',2002:2:2024,'FontSize',9);
%ylabel({'Predicted incidence of RVGE';'(per month)'})
