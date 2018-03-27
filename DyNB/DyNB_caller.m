function result = DyNB_caller()

% let us read the data
% two conditions, with five timepoints and three replicates
% moreover, we have size factors for each individual data point
% (given by DESeq), and labels for the genes
load('../data/DyNB_example_data.mat')

% let us define time poinsts
% please scale the time points so that they are in the interval [0,1]
% measurement time points
t = [0 12 24 48 72]'/72;
% timepoints for evaluating Gaussian processes
t_star = [0:72]'/72;

% get the coefficients from the dispersion estimation
coefficients = estimate_dispersion(th0_data,th17_data,th0_sizefactors,th17_sizefactors);

% you should use distributed computing here because this will take some time
%for idx=1:length(th0_data)
%	result = DyNB(th0_data{idx},th17_data{idx},t,t_star,th0_sizefactors,th17_sizefactors,coefficients);
%end
result = DyNB(th0_data{15270},th17_data{15270},t,t_star,th0_sizefactors,th17_sizefactors,coefficients);
end

function coefficients = estimate_dispersion(data1,data2,data1_sizefactors,data2_sizefactors)
% estimate dispersion
x = [];
y = [];
for idx=1:length(data1)
	% NOTICE: the timepoint t=0 is shared here
	data = [reshape(data1{idx}./data1_sizefactors,[numel(data1{idx}) 1]); reshape(data2{idx}(:,2:end)./data2_sizefactors(:,2:end),[numel(data2{idx}(:,2:end)) 1])];
	q = mean(data)';
	w = mean(var([data1{idx}./data2_sizefactors data2{idx}(:,2:end)./data2_sizefactors(:,2:end)]))';
	x = [x; q];
	y = [y; w];
end
z = x*.1/length(data).*sum(sum(1./[data1_sizefactors data2_sizefactors(:,2:end)]));
tmp = y;
X = [ones(size(x(x > 0 & tmp > 0 & ~isinf(tmp)))) log10(x(x > 0 & tmp > 0 & ~isinf(tmp))) log10(x(x > 0 & tmp > 0 & ~isinf(tmp))).^2];
Y = log10(tmp(x > 0 & tmp > 0 & ~isinf(tmp)));
coefficients = (X'*X)\X'*Y;
end
