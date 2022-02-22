function [deltaW] = chronotron_learn(input_times, input_neurons, y0, t, W, K, tau_m,theta, tau_q, eta, gamma_r)
%% Chronotron Learning Function
N = size(W,2);
F_target = y0;
sigma       = @(X) ((X^2)/2);
[v_spk,spk] = IF_sim(input_times, input_neurons, t, W, K, tau_m, theta);
F_source = spk; 
[F_add, F_rm, F_mv, d] = VP_spike_classify(F_source, F_target, tau_q, sigma);
[lambda] = IF_get_lambdas(N, input_times, input_neurons, t, K);
if ~isempty(F_rm)  % neurons to remove
        remove = F_source(F_rm); 
    for j=1:length(remove)
        remove(j) = find(t==remove(j)); 
    end
    rmSum = sum(lambda(:,remove),2);
else% no neurons to remove
    rmSum = 0;

end
if ~isempty(F_add) % neurons to add
     F_target = round(F_target,4);
    Spikeadd = F_target(F_add);
    for j=1:length(Spikeadd)
        Spikeadd(j) = find(t == Spikeadd(j)); 
    end
    addSum = sum(lambda(:,Spikeadd),2);
else % no neurons to add
      addSum = 0;

end

if ~isempty(F_mv) % neurons to move
       spk_mv_A = F_source(F_mv(1,:)); % Spike point A
    sp_mv_B = F_target(F_mv(2,:)); % Spike point B
    for j=1:length(sp_mv_B)
        mv(j) = find(t==spk_mv_A(j));
    end
    mvSum = sum((spk_mv_A - sp_mv_B).*(lambda(:,mv).^2),2);
else % no neurons to move
    mvSum = 0;
end
deltaW = eta*(addSum - rmSum + (gamma_r/(tau_q)^2).* mvSum);
deltaW = deltaW';

end

