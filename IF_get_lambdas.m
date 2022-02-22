function [lambdas] = IF_get_lambdas(N, input_times, input_neurons, t, K)
%% Lambda(t) for I&F model:

lambdas = zeros(N,length(t));

% for n=1 : N
    for i=1: length(input_times)
        lambdas(input_neurons(i),:) = lambdas(input_neurons(i),:)+ K(t-input_times(i));
%         for j=1 : length(t)
%             lambdas(j) = lambdas(j)+ K(t(j)-input_times(i)); %Lambda(t) as shown in introduction
%             
%         end
    end
% end

end
