%(N=3)
% n, offline, online
n_timing =[3, 3.851396904, 2.139412408;
4, 5.744539697, 2.548018473;
5, 8.317854901, 3.243462204;
6, 11.053912582, 4.027396481;
7, 14.663377918, 4.909351186;
8, 18.303002231, 6.009072964;
9, 22.711086936, 7.470598776;
10, 27.666474141, 9.156292053;
15, 70.623225393, 22.868706809]
%25, 249.361241122, 30.500706989];

%(n=5)
N_timing =[2, 0.3096907;
3, 1.891401442;
4, 24.632279412;
5, 352.611930983];
%(n=5)
N_timing2 = [2, 0.112225146;
3, 2.057751418;
4, 27.154251206;
5, 374.472727385];

N_avg =[
    2.0000,    0.2110;
    3.0000,    1.9746;
    4.0000,   25.8933;
    5.0000,  363.5423];


data = csvread('203_project_data.csv');
iter_timing = csvread('203_project_timing.csv');

iter = data(:,1);
feas = data(:,2);
rand_index = data(:,3);
state = data(:,4:5);
normstate = state(:,1).^2 + state(:,2).^2;
control = data(:,6:7);
r  = data(:,8);
cost = data(:,9);
gamma = data(:,10);
expected_cost = data(:,11);


% Defaults for this blog post
width = 5;     % Width in inches
height = 5;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 11;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize


% Runtime versus number of states - N = 3
figure(1);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
plot(n_timing(:,1),n_timing(:,2),'bo-',n_timing(:,1), n_timing(:,3),'rx--','LineWidth',lw,'MarkerSize',msz); %<- Specify plot properites
h_legend =legend('Offline', 'Online', 'Location', 'NorthWest');
set(h_legend,'FontSize',12,'FontName','Times');
xlabel('Number of states','FontSize',12,'FontName','Times');
ylabel('Runtime (s)','FontSize',12,'FontName','Times')
%title('Runtime versus number of states (look-ahead = 3)','FontSize',12,'FontName','Times');
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
print('203_runtime_n','-dpng','-r300');

% Runtime versus N, (m = 5)
figure(2);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
plot(N_avg(:,1),N_avg(:,2),'ko-','LineWidth',lw,'MarkerSize',msz); %<- Specify plot properites
xlim([1 6]);
xlabel('Look-ahead distance','FontSize',12,'FontName','Times');
ylabel('Runtime (s)','FontSize',12,'FontName','Times')
%title('Runtime versus Look-ahead (number of states = 5)','FontSize',12,'FontName','Times');
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
print('203_runtime_N','-dpng','-r300');

% Error versus iteration
figure(3);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
plot(iter(1:49),sqrt(normstate(1:49)),'-','LineWidth',lw,'MarkerSize',msz); %<- Specify plot properites
xlabel('Iteration','FontSize',12,'FontName','Times');
ylabel('||X||','FontSize',12,'FontName','Times')
%title('Runtime versus Look-ahead (number of states = 5)','FontSize',12,'FontName','Times');
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
print('203_convergence','-dpng','-r300');

% Time per iteration
figure(4);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
plot(iter_timing(1:49),'ko-','LineWidth',lw,'MarkerSize',msz); %<- Specify plot properites
xlabel('Iteration','FontSize',12,'FontName','Times');
ylabel('Runtime (s)','FontSize',12,'FontName','Times')
%title('Runtime versus Look-ahead (number of states = 5)','FontSize',12,'FontName','Times');
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
print('203_runtime_i','-dpng','-r300');

