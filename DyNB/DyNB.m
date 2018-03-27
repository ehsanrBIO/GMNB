function result = DyNB(data1,data2,t,t_star,data1_sizefactors,data2_sizefactors,coefficients)
%DyNB  The DyNB method is a method for quantifying differential
%       time-course dynamics of time-course RNA-seq datasets.
%       This implementation analyzed RNA-seq data without
%       taking into account the possibility of time scaling.
%   
%   usage: result = DyNB(data1,data2,t,t_star,data1_sizefactors,data2_sizefactors,coefficients)
%
%   where:
%       result          is the result vector, where
%		 .bayes_factor  the estimated Bayes factor for differential expression
%                       (individual models / shared model)
%        .data1_normalized
%                       normalized read counts (using the size factors) of
%                       the first condition
%        .data2_normalized
%                       normalized read counts (using the size factors) of
%                       the second condition
%		 .f
%                       posterior means of different f's evaluated at t,
%                       that is, f{1} is from the shared model (data1 and data2),
%                       f{2} is from the individual model of data1, and f{3} is
%                       from the individual model of data2
%		 .f_star
%                       posterior means of different f's evaluated at t_star,
%                       that is, f_star{1} is from the shared model (data1 and data2),
%                       f_star{2} is from the individual model of data1, and f_star{3}
%                       is from the individual model of data2
%		 .variance
%						variances of the posterior means of different f's evaluated at t
%                       (corresponds to f, indexing is the same)
%		 .variance_star
%						variances of the posterior means of different f's evaluated at t_star
%                       (corresponds to f_star, indexing is the same)
%
%       data1           is a data matrix consisting the raw red counts of
%                       the first condition
%                       (R replicates x N timepoints)
%       data2           is a data matrix consisting the raw read counts of
%                       the second condition
%                       (R replicates x N timepoints)
%       t     			is a timepoint vector (corresponds to data1 
%                       and data2)
%                       (N timepoints x 1)
%       t_star     		is a timepoint vector for evaluating Gaussian
%                       processes, can be more dense that the measurement
%                       grid but has to hold the actual measurement timepoints
%                       (M timepoints x 1)
%       data1_sizefactors
%                       is a size factor matrix of data1 (size factors
%                       per data set are given by DESeq)
%                       (R replicates x N timepoints)
%       data2_sizefactors
%                       is a size factor matrix of data2 (size factors
%                       per data set are given by DESeq)
%                       (R replicates x N timepoints)
%       coefficients    is a vector holding the coefficients from
%                       the dispersion estimation (1 x 3)
%                       the experimental conditions
%
%   Author: Tarmo ?ij? <tarmo.aijo@aalto.fi>
%
%   Updated: 08/11/2014

% let us check that the dimensions of the input variables make sense
t_indices = arrayfun(@(x) find(x == t_star),t);
if ~issame(t,t_star(t_indices))
	disp('ERROR: t should be included in t_star!')
	return
end
if size(data1,2) ~= length(t) || size(data2,2) ~= length(t)
	disp('ERROR: size(data1,2) and size(data2,2) should be equal to length(t)')
	return
end
if ~issame(size(data1),size(data1_sizefactors)) || ~issame(size(data2),size(data2_sizefactors))
	disp('ERROR: size(data1) should be equal to size(data1_sizefactors) and size(data2) should be equal to size(data2_sizefactors)!')
	return
end

% let us skip the genes with zero reads in one condition or in both conditions
if (sum(data1(:) == 0) == numel(data1)) || (sum(data2(:) == 0) == numel(data2))
    disp('WARNING: not enough reads!')
    return
end

% number of MCMC iterations 
n_iterations = 1e5;

% initialize output variables
accepted_f = cell(3,1);
marginal_likelihood = zeros(3,1);
for combinationIdx=1:3
	accepted_f{combinationIdx} = zeros(length(t),n_iterations);
	accepted_f_star{combinationIdx} = zeros(length(t_star),n_iterations);
end

for combinationIdx=1:3
	% shared model
    if combinationIdx == 1
        k = [data1;data2];
        sizefactors = [data1_sizefactors;data2_sizefactors];
	% individual model (data1)
    elseif combinationIdx == 2
        k = data1;
        sizefactors = data1_sizefactors;
	% individual model (data2)
    elseif combinationIdx == 3
        k = data2;
        sizefactors = data2_sizefactors;
    end
    
    % parameter intervals (l)
    parameter_intervals = [5e-1 1e0];

	% let us define the prior for f
    GP_prior_mean = mean([min(k(:)),max(k(:))])*ones(size(t));
    GP_prior_K = get_covariance_matrix(t,t,[500*mean([min(k(:)),max(k(:))]),0.75]);
    
    % proposal distribution, standard deviations (l)
    s = max(parameter_intervals)/100;
    
    % get a hyperparameters sample (sigma, l, mu)
    hyperparameter_sample = [10*sqrt(var(k(:)));1e0;median(k(:))];
    
    % calculate the covariance matrix
    K = get_covariance_matrix(t,t,hyperparameter_sample(1:2));
    K_star = get_covariance_matrix(t_star,t_star,hyperparameter_sample(1:2));
    
    % get an initial f sample 
    if median(k(:)) ~= 0
        foobar = median(k(:))*ones(size(t_star));
    else % if the median of read counts is zero, then the f sample is the vector of ones
        foobar = ones(size(t_star));
    end
    GP_sample_star = get_GP_sample(t_star,K_star,foobar);
    GP_sample = GP_sample_star(t_indices);
      
    % calculate the value of the likelihood given the sampled f
    old_likelihood_value = get_NB_likelihood(GP_sample,k,coefficients,sizefactors);

    % initialize the result vectors
    accepted_samples = zeros(1,n_iterations);

	% MCMC loop
    mcmcIdx = 1;
    while (mcmcIdx <= n_iterations)
        u = rand(1);

        % get a hyperparameters sample (sigma, l, mu)
        hyperparameter_sample_new = [10*sqrt(var(k(:)));get_new_hyperparameter_sample(hyperparameter_sample,parameter_intervals,s);median(k(:))];

        % calculate the covariance matrices 
        K_new = get_covariance_matrix(t,t,hyperparameter_sample_new(1:2));
        K_new_star = get_covariance_matrix(t_star,t_star,hyperparameter_sample_new(1:2));

        % sample a new f sample
        GP_sample_new_star = get_GP_sample(t_star,K_new_star,GP_sample_star);
        GP_sample_new = GP_sample_new_star(t_indices);
        
        % calculate the probability of getting the sampled hyperparameters sample
        % given the old hyperparameters and vice versa
        [prob_from_new_to_old_hp prob_from_old_to_new_hp] = get_hyperparameter_jump_probability(hyperparameter_sample_new,hyperparameter_sample,parameter_intervals,s);
        
        % calculate the probability of the new f sample given the old
        % hyperparameters sample        
        prob_new_GP_w_old_hyperparameters = mvnpdf(GP_sample_new,GP_sample,K);
        
        % calculate the probability of the old f sample given the new
        % hyperparameters sample
        prob_old_GP_w_new_hyperparameters = mvnpdf(GP_sample,GP_sample_new,K_new);
        
        % calculate the value of the likelihood given the new f sample
        new_likelihood_value = get_NB_likelihood(GP_sample_new,k,coefficients,sizefactors);
        
		% if we have encountered a numerical issue, then let us start from the beginning
        if isnan(new_likelihood_value) | isinf(new_likelihood_value)
            mcmcIdx = 1;
            % get a hyperparameters sample
            hyperparameter_sample = [10*sqrt(var(k(:)));1e0;median(k(:))];

            % calculate the covariance matrix
            K = get_covariance_matrix(t,t,hyperparameter_sample(1:2));
            K_star = get_covariance_matrix(t_star,t_star,hyperparameter_sample(1:2));

            % get an initial f sample 
            if median(k(:)) ~= 0
                foobar = median(k(:))*ones(size(t_star));
            else % if the median of read counts is zero, then the f sample is the vector of ones
                foobar = ones(size(t_star));
            end

        	% sample a new f sample
            GP_sample_star = get_GP_sample(t_star,K_star,foobar);
            GP_sample = GP_sample_star(t_indices);

            % calculate the value of likelihood given the sampled f
            old_likelihood_value = get_NB_likelihood(GP_sample,k,coefficients,sizefactors);
            continue
        end

        prior_ratio = mvnpdf(GP_sample_new,GP_prior_mean,GP_prior_K)/mvnpdf(GP_sample,GP_prior_mean,GP_prior_K);
         
	    % let us calculate the value of the metropolis-hasting criteria
        if u < min(1,prior_ratio*(prob_old_GP_w_new_hyperparameters*prob_from_new_to_old_hp*new_likelihood_value)/(prob_new_GP_w_old_hyperparameters*prob_from_old_to_new_hp*old_likelihood_value))
			% store the Gaussian process samples
			accepted_f{combinationIdx}(:,mcmcIdx) = GP_sample_new;
			accepted_f_star{combinationIdx}(:,mcmcIdx) = GP_sample_new_star;
            
            % take the value of the integrand given the sample
            accepted_samples(mcmcIdx) = new_likelihood_value;
            
            % update the sample variables
            GP_sample_star = GP_sample_new_star;
            GP_sample = GP_sample_new;
            hyperparameter_sample = hyperparameter_sample_new;
            K = K_new;
            K_star = K_new_star;
            old_likelihood_value = new_likelihood_value;  
        else
			accepted_f{combinationIdx}(:,mcmcIdx) = GP_sample;
			accepted_f_star{combinationIdx}(:,mcmcIdx) = GP_sample_star;
            accepted_samples(mcmcIdx) = old_likelihood_value;
        end
    	mcmcIdx = mcmcIdx + 1;
    end
    marginal_likelihood(combinationIdx) = 1/(1/length(accepted_samples)*sum(1./accepted_samples(1:end)));
end

% output variables
result.bayes_factor = prod(marginal_likelihood(2:3))/marginal_likelihood(1);

result.data1_normalized = data1./data1_sizefactors;
result.data2_normalized = data2./data2_sizefactors;

for combinationIdx=1:3
	result.f{combinationIdx} = mean(accepted_f{combinationIdx},2);
	result.f_star{combinationIdx} = mean(accepted_f_star{combinationIdx},2);
	result.variance{combinationIdx} = var(accepted_f{combinationIdx},0,2);
	result.variance_star{combinationIdx} = var(accepted_f_star{combinationIdx},0,2);
end
end

function [y1 y2] = get_hyperparameter_jump_probability(x1,x2,parameter_intervals,s)
y1 = truncated_normal_distribution_pdf(x2(2),x1(2),s,parameter_intervals(1),parameter_intervals(2));
y2 = truncated_normal_distribution_pdf(x1(2),x2(2),s,parameter_intervals(1),parameter_intervals(2));
end

function y = truncated_normal_distribution_pdf(x,mu,sigma,a,b)
y = normpdf((x-mu)/sigma) / (sigma*(normcdf((b-mu)/sigma)-normcdf((a-mu)/sigma)));
end

function y = get_new_hyperparameter_sample(x,parameter_intervals,s)
y = zeros(1,1);
y(1) = randraw('normaltrunc',[parameter_intervals(1), parameter_intervals(2), x(2), s]);
end

function y = get_GP_sample(t,K,mean_vector)
chol_K = chol(K);
% rejection sampling
while true
    y = mean_vector+(randn(1,length(t))*chol_K)';
    if sum(y < 0) == 0
        break
    end
end 
end

function y = get_covariance_matrix(t,tp,hyperparameters)
t = repmat(t,[1 size(tp,1)]);
tp = repmat(tp',[size(t,1) 1]);
y = hyperparameters(1)*exp(-1/hyperparameters(2)^2*(t-tp).^2);
% for improving numerical stability
y = y + hyperparameters(1)*1e-12*eye(size(y));
end

function y = get_NB_likelihood(mu,k,coefficients,sizefactors)
k = k';
y = 1;
for replicateIdx=1:size(k,2)
    X = [ones(size(mu)) log10(mu) log10(mu).^2];
    dispersion = 2*10.^(X*coefficients);
    
    s2 = mu.*sizefactors(replicateIdx,:)'+dispersion;
    tmpR = (mu.*sizefactors(replicateIdx,:)').^2./(s2-(mu.*sizefactors(replicateIdx,:)'));
    tmpP = (s2-(mu.*sizefactors(replicateIdx,:)'))./(s2);
    
    if ~isreal(tmpR) || ~isreal(tmpP)
        y = Inf;
        return
    end
    
    tmp2 = (1-tmpP).^tmpR;
    if sum(k(:,replicateIdx) == 0) > 0
        zero_indices = k(:,replicateIdx) == 0;
        tmp1 = gamma(tmpR);
        y = prod([y; (gamma(k(zero_indices,replicateIdx)+tmpR(zero_indices)) ./ (arrayfun(@(x) prod(1:x),k(zero_indices,replicateIdx)).*tmp1(zero_indices)) .* tmp2(zero_indices) .* tmpP(zero_indices).^(k(zero_indices,replicateIdx)))]);
        y = prod([y; (1./k(~zero_indices,replicateIdx) .* 1./beta(k(~zero_indices,replicateIdx),tmpR(~zero_indices)) .* tmp2(~zero_indices) .* tmpP(~zero_indices).^(k(~zero_indices,replicateIdx)))]);
    else
        y = prod([y; (1./k(:,replicateIdx) .* 1./beta(k(:,replicateIdx),tmpR) .* tmp2 .* tmpP.^(k(:,replicateIdx)))]);
    end
end
end
