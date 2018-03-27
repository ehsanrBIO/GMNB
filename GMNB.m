function result = GMNB(data1, data2)

%GMNB   The GMNB is a dynamic differential expression analysis method
%       for temporal RNA-seq count data.
%   
%       data1           is a data matrix consisting the raw red counts of
%                       the first condition
%                       (replicates x timepoints)
%       data2           is a data matrix consisting the raw read counts of
%                       the second condition
%                       (replicates x timepoints)
%
%   Author: Ehsan Hajiramezanali <ehsanr@tamu.edu>
%   Paper:  Differential expression analysis of dynamical sequencing count data with a gamma Markov chain


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% model parameters
ezero=1.0; fzero=1.0; 
rzeroInt=100; 
hzero = 1.0;
marginal_likelihood = zeros(3,1);
a0 = .01; b0 = .01;

%% Gibbs sampling specific parameters
burnin  = 1000; collection = 1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for combinationIdx=1:3
  % shared model
    if combinationIdx == 1
        for k=1:length(data1)
          n_j_t_k{k} = [data1{k}; data2{k}];
        end
  % individual model (data1)
    elseif combinationIdx == 2
        for k=1:length(data1)
          n_j_t_k{k} = data1{k};
        end
  % individual model (data2)
    elseif combinationIdx == 3
        for k=1:length(data1)
          n_j_t_k{k} = data2{k};
        end
    end 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% for the proposed model
  [J,T]    = size(n_j_t_k{1});
  if ~exist('est_r','var')
    est_r = repmat({zeros(3,T)}, 1, length(n_j_t_k));
  end
  r = repmat({ones(T,1)}, 1, length(n_j_t_k));
  rzero = repmat({rzeroInt}, 1, length(n_j_t_k));
  p_j_t = 1/2 * ones(J,T);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Samples = repmat({[]}, 1, length(n_j_t_k));
  likelihood_value = repmat({[]}, 1, length(n_j_t_k));
  
  for iter=1:burnin + collection
    iter
    parfor k = 1:length(n_j_t_k)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% sample c's first; T-specific variables
      c_k{k} = gamrnd(ezero + rzero{k} + sum(r{k}(1:end-1)),1./(fzero+sum(r{k})));
      L_j_t{k} = zeros(J,T);
      L_t{k} = zeros(1,T);
      U_t{k} = zeros(1,T);
      q_t{k} = zeros(1,T);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% backward sampling; only the L's which are T-specific variables
      for t=T:-1:2
          [~,L_j_t{k}(:,t)] = CRT(n_j_t_k{k}(:,t),r{k}(t));
          L_t{k}(t) = sum(L_j_t{k}(:,t));
          [~,U_t{k}(t-1)] = CRT(L_t{k}(t)+U_t{k}(t),r{k}(t-1));
          q_t{k}(t-1)     = ((- log(max((1-q_t{k}(t)),eps)) - sum(log(max((1-p_j_t(:,t)),eps)))) ./ max(c_k{k} - log(max((1-q_t{k}(t)),eps)) - sum(log(max((1-p_j_t(:,t)),eps))),eps));
      end
      [~,L_j_t{k}(:,1)] = CRT(n_j_t_k{k}(:,1),r{k}(1,:));
      L_t{k}(1) = sum(L_j_t{k}(:,1));
      [~,U_zero{k}] = CRT(L_t{k}(1)+U_t{k}(1),rzero{k});
      q_zero{k} = ((- log(max((1-q_t{k}(1)),eps)) - sum(log(max((1-p_j_t(:,1)),eps)))) ./ max(c_k{k} - log(max((1-q_t{k}(1)),eps)) - sum(log(max((1-p_j_t(:,1)),eps))),eps));
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% forward sampling; only the theta's which are T-specific variables
%      c_0       = 1;
      rzero{k}  = gamrnd(hzero + U_zero{k},1./max(c_k{k}-log(max(1-q_zero{k},eps)),eps));
%      rzero{k}  = gamrnd(hzero + U_zero{k},1./max(c_0-log(max(1-q_zero{k},eps)),eps));
      r{k}(1) = gamrnd(rzero{k} + U_t{k}(1) + L_t{k}(1), 1./max(c_k{k} - sum(log(max((1-p_j_t(:,1)),eps))) - log(max(1-q_t{k}(1),eps)),eps));
      for t=2:T-1
          param1     = r{k}(t-1) + U_t{k}(t) + L_t{k}(t);
          param2     = (1./max(c_k{k} - sum(log(max((1-p_j_t(:,t)),eps))) - log(max((1-q_t{k}(t)),eps)),eps));
          r{k}(t) = gamrnd(param1,param2);
      end
      param1     = r{k}(T-1) + L_t{k}(T);
      param2     = (1./(max(c_k{k} - sum(log(max(1-p_j_t(:,T),eps))),eps)));
      r{k}(T) = gamrnd(param1,param2);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if(iter>burnin)
   
          %% for proposed model
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          if (mod(iter,10)==0)
            likelihood_value{k} = [likelihood_value{k} get_NB_loglikelihood(n_j_t_k{k}, r{k}, p_j_t)];
            Samples{k} = [Samples{k} r{k}];
          end
      end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p_j_t = betarnd(a0 + sum(cat(3, n_j_t_k{:}), 3), b0 + repmat(sum(cell2mat(r), 2)', J, 1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end
  
  for k = 1:length(n_j_t_k)
    marginal_log_likelihood{k}(combinationIdx) = log(sum(exp(likelihood_value{k} - max(likelihood_value{k})))) + max(likelihood_value{k}) - log(length(likelihood_value{k}));
    est_r{k}(combinationIdx,:) = mean(Samples{k},2);
    var_r{k}(combinationIdx,:) = var(Samples{k},0,2);
  end
  PJT{combinationIdx}=p_j_t;
end

% output variables
for k = 1:length(n_j_t_k)
  result{k}.log_bayes_factor = sum(marginal_log_likelihood{k}(2:3)) - marginal_log_likelihood{k}(1);
  result{k}.bayes_factor = exp(result{k}.log_bayes_factor);
  result{k}.est_r = est_r{k};
  result{k}.var_r = var_r{k};
  result{k}.PJT = PJT;
  result{k}.marginal = marginal_log_likelihood{k};
end

end
