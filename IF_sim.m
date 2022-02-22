function [V, spk_times] = IF_sim(input_times, input_neurons, t, W, K, tau_m, theta)
%% Spike Response Model for I&F V(t):
N = size(W,2);
lambdas = IF_get_lambdas(N,input_times, input_neurons, t, K);
V = zeros(size(t));
N = size(W,2);
spk_times = [];
for i = 1:length(t)
    if isempty(spk_times)
        V(i) = (W*lambdas(:,i));
    else
        V(i) = ((W*lambdas(:,i)) - theta*sum(exp(-(t(i) - spk_times)/tau_m))); %Identity from introduction
    end
    if V(i) >= theta %If Action Potential has occured...
        spk_times = [spk_times , t(i)]; %...mark them
    end
end

end

