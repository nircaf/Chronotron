
% TODO: Student1, Nir Cafri 313160194
% TODO: Student2, Aviv Azmanov 204805469

% Introduction to Computation and Cognition
% See PDF document for instructions

%%%%%%%%%%%%%%%%%%%%%%%
clear;close all;clc;
tic;
% loading data
% load first_test_data;

load train_data;
%% define parameters:
eta = 0.005;
tau_m = 0.02; %[s]
tau_s = (tau_m/4); %[s]
tau_q = tau_m; %[s]
theta = 30; % mVolt
gamma_r = tau_q;
t = (0:0.0001:0.5); %Time from 0s to 0.5s by (1*10^-4)s
t = round(t,4);
VPD = 0;
W=rand(1,N)*(theta/2); % Random W starting weights
epochs = 50; % How many times chronotron repeats to learn

%% Define the (normalized) kernel function
alpha   = (tau_m/tau_s);
kappa   = alpha^(alpha/(alpha-1))/(alpha-1);
K       = @(t) (t > 0).*kappa.*(exp(-t/tau_m) - exp(-t/tau_s));
sigma       = @(X) ((X^2)/2);
%% Learning
for (epo = 1:epochs)
    samp_order = randperm(length(Samples)); % Choose a random sample
    for (samp = 1:length(Samples))
        example = Samples(samp_order(samp));
        input_times = example.times;
        input_neurons = example.neurons;
        y0 = example.y0;
        [deltaw] = chronotron_learn(input_times, input_neurons, y0, t, W, K, tau_m, theta, tau_q, eta, gamma_r);
        W = W+deltaw; % Update the weigths
        if (epochs==epo)
            [V, spk_times] = IF_sim(input_times, input_neurons, t, W, K, tau_m, theta); %Find the spike times
            F_source = spk_times;
            F_target = y0;
            [F_add, F_rm, F_mv, d] = VP_spike_classify(F_source, F_target, tau_q, sigma);
            VPD = [VPD , d]; % Update victor Purpura distance
        end
    end
end

%% Plotting
for (j =1:3) % 3 Samples
    Vgraph2 = zeros(3,length(t));
    randsamp = randi([1 length(Samples)],1);
    input_neurons = Samples(randsamp).neurons;
    input_times = Samples(randsamp).times;
    y0 = Samples(randsamp).y0;
    F_target = y0;
    Tgraph2 = [];
    [Vgraph2(j,:), spk] = IF_sim(input_times, input_neurons, t, W, K, tau_m, theta); % V/T
    F_source = spk;
    
    [F_add, F_rm, F_mv, d] = VP_spike_classify(F_source, F_target, tau_q, sigma); % VP distance
    subplot(3,1,j);
    hold on;
    %% V post learning
    plot(t, Vgraph2(j,:),'Color',[0.5,0.5,0.9],'DisplayName','Voltage');
    
    %% Spikes post learning
    if (~isempty(spk))
        for s = 1:length(spk)
            Tgraph2(s) = find(t == spk(s)); %
        end
        scatter(spk,Vgraph2(j,Tgraph2),70,'MarkerEdgeColor',[0.1 0.1 0.1],...
            'MarkerFaceColor',[.1 0.9 .1],...
            'LineWidth',0.5,'DisplayName','Chronotron');
    end
    %% Teacher
    if (~isempty(y0))
        yy0 = linspace(theta, theta, length(y0));
        scat = scatter(y0,yy0,15,'MarkerEdgeColor',[0 0 0],...
            'MarkerFaceColor',[.1 0.1 .9],...
            'LineWidth',0.1,'DisplayName','Theacher');
    end
    %% miscelaneous
    ylabel('Voltage [mV]');
    title(['V post learning', num2str(j)]);
    str = {['VP distance: ',num2str(d)]};
    
    xlabel('Time [s]');
    
    yline(theta,'-.r','Theta','DisplayName','Theta');
    
    text(0,0,str);
    legend
    hold off;
end

VPD_new = VPD(2:1:end);
VPmean=mean(VPD_new);
standarddeviation=std(VPD_new);
fprintf ('Mean distance(s): %g \n Standard deviation(s): %g \n Learning time(min): %g \n',VPmean,standarddeviation,toc/60);