% Plots the Sea Level Equivalent (SLE) from the 240 000 yr integration 
% calculated in CalcDAIS. The bars in the subplots are model 
%targets from data reconstructions. Total sea level in plotted as black dashed
%lines in the lower two subplots. Like Figure 6 in the DAIS publication
 
load ('OutDais')

figure(1)

%240 kyr BP to 2010 AD hindcast of AIS SLE
box on
subplot(2,2,1)

set(gca,'xtick',[-200000 -150000 -100000 -50000  0])
set(gca,'xticklabel',[ -200 -150 -100 -50 0],'fontsize',14)
xlabel('Date (kyr BP)','fontsize',16) 
set(gca,'ytick',[-20 -15 -10 -5 0 5 10],'yaxislocation', 'left')
set(gca,'yticklabel',[-20 -15 -10 -5 0 5 10], 'fontsize',14)
ylabel('SLE (m)','fontsize',16)
axis([-245000 5000 -19 7])
hold on
plot(date(1:240010),SLE(1:240010)-mean(SLE(239961:239990)),'b','LineWidth',2);hold on
line([date(110000) date(130000)],[2.5 2.5],'Color', 'k','LineWidth',2) 
line([date(110000) date(130000)],[6 6],'Color', 'k','LineWidth',2) 
line([date(120000) date(120000)],[2.5 6],'Color', 'k','LineWidth',2) 
line([date(1) date(240010)],[0 0],'Color', 'k','LineWidth',0.5) 
text(-237000, 5.3,'a','fontweight','bold','fontsize',14)

%25 kyr BP to 2010 AD hindcast of AIS SLE
box on
subplot(2,2,2)
set(gca,'xtick',[-25000 -20000 -15000 -10000 -5000 0])
set(gca,'xticklabel',[-25 -20 -15 -10 -5 0],'fontsize',14)
xlabel('Date (kyr BP)','fontsize',16) 
set(gca,'ytick',[-20 -15 -10 -5 0 5 10],'yaxislocation', 'left')
set(gca,'yticklabel',[-20 -15 -10 -5 0 5 10], 'fontsize',14)
ylabel('SLE (m)','fontsize',16)
axis([-26000 1000 -19 7])
hold on
plot(date(215000:240010),SLE(215000:240010)-mean(SLE(239961:239990)),'b','LineWidth',2);hold on
line([date(21500) date(240010)],[0 0],'Color', 'k','LineWidth',0.5) 
line([date(220000) date(220000)],[-8 -17],'Color', 'k','LineWidth',2) 
line([date(219000) date(221000)],[-8 -8],'Color', 'k','LineWidth',2) 
line([date(219000) date(221000)],[-17 -17],'Color', 'k','LineWidth',2) 
text(-25200, 5.3,'b','fontweight','bold','fontsize',14)

%6 kyr BP to 2010 AD hindcast of AIS SLE
subplot(2,2,3)
box on
set(gca,'xtick',[-6000 -4000 -2000 0])
set(gca,'xticklabel',[-6000 -4000 -2000 0],'fontsize',14)
xlabel('Date (kyr BP)','fontsize',16) 
set(gca,'ytick',[-6 -5 -4 -3 -2 -1 0 1],'yaxislocation', 'left')
set(gca,'yticklabel',[-6 -5 -4 -3 -2 -1 0 1], 'fontsize',14)
ylabel('SLE (m)','fontsize',16)
axis([-6500 300 -6 1])
hold on
plot(date(234000:240010),SLE(234000:240010)-mean(SLE(239961:239990)),'b','LineWidth',2);hold on
plot(date(234000:240010),SL(234000:240010), 'k--','LineWidth',2)
line([date(21500) date(240010)],[0 0],'Color', 'k','LineWidth',0.5) 
line([date(234000) date(234000)],[-2 -4],'Color', 'k','LineWidth',2) 
line([date(233800) date(234200)],[-2 -2],'Color', 'k','LineWidth',2) 
line([date(233800) date(234200)],[-4 -4],'Color', 'k','LineWidth',2) 
text(-6350, 0.6,'c','fontweight','bold','fontsize',14)

%1800 to 2010 AD hindcast of AIS SLE
box on
subplot(2,2,4)
set(gca,'xtick',[-200 -100 0])
set(gca,'xticklabel',[1800 1900 2000],'fontsize',14)
xlabel('Date (yr AD)','fontsize',16) 
set(gca,'ytick',[-0.2 -0.1 0 0.1],'yaxislocation', 'left')
set(gca,'yticklabel',[-0.2 -0.1 0 0.1], 'fontsize',14)
ylabel('SLE (m)','fontsize',16)
axis([-220 30 -0.21 0.1])
hold on
plot(date(239800:240010),SLE(239800:240010)-mean(SLE(239961:239990)),'b','LineWidth',2);hold on
plot(date(239800:240010),SL(239800:240010), 'k--','LineWidth',2)
line([date(239800) date(240010)],[0 0],'Color', 'k','LineWidth',0.5) 
text(-214, 0.08,'d','fontweight','bold','fontsize',14)

