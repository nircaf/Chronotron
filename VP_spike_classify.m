function [Fadd, Frm, Fmv, d] = VP_spike_classify(F_source, F_target, tau_q, sigma)
%% Victor-Purpura Spike Classifier Algorithm as shown in project introduction:
n_source = length(F_source); %Input
n_target = length(F_target); %Output
D = zeros(n_source+1, n_target+1); %D Matrix
sadd = cell(n_source+1, n_target+1); %Spikes to add
srm = cell(n_source+1, n_target+1); %Spikes to remove
smv = cell(n_source+1, n_target+1); %Spikes to move

%% Row #1
for j = 1:(n_target+1) 
    D(1,j) = j-1; 
    sadd{1,j} = [1:j-1];
end

%% Column #1
for i = 2:(n_source+1) 
    D(i,1) = i-1; 
    srm{i,1} = [srm{i-1,1}, i-1];
end

%% Classifier for next J-1 columns:
for i = 2:(n_source+1) 
    for j = 2:(n_target+1) 
        a1 = D(i-1, j) + 1;
        a2 = D(i, j-1) + 1;
        a3 = D(i-1, j-1) + sigma(abs(F_source(i-1) - F_target(j-1))/tau_q);
        if (a1<=a2) && (a1<=a3)
            D(i,j) = a1;
            sadd{i,j} = sadd{i-1,j};
            srm{i,j} = [srm{i-1,j}, i-1];
            smv{i,j} = smv{i-1,j};
        elseif a2<=a3
            D(i,j) = a2;
            sadd{i,j} = [sadd{i,j-1},j-1];
            srm{i,j} = srm{i,j-1};
            smv{i,j} = smv{i,j-1};
        else
            D(i,j) = a3;
            sadd{i,j} = sadd{i-1,j-1};
            srm{i,j} = srm{i-1,j-1};
            smv{i,j} = [smv{i-1,j-1},[(i-1);(j-1)]];
        end
    end 
end
Fmv = smv{n_source+1,n_target+1};
Frm = srm{n_source+1,n_target+1};
Fadd = sadd{n_source+1,n_target+1};
d = D(n_source+1,n_target+1); %The last cubicle 'd' of 'D' is the VP-Distance

end

