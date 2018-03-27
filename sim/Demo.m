% GMNB for the simulated data
if 1
  for idx=1:20
    load(sprintf('data_simAR_%d.mat', idx))
    result = DyNB_MG_Beta(th0_data,th17_data);
    save(sprintf('result_data_simAR_%d.mat', idx),'result')
  end
end
if 0
  for idx=1:20
    load(sprintf('data_simGP_%d.mat', idx))
    result = DyNB_MG_Beta(th0_data,th17_data);
    save(sprintf('result_data_simGP_%d.mat', idx),'result')
  end
end
if 0
  for idx=1:20
    load(sprintf('data_simNB_%d.mat', idx))
    result = DyNB_MG_Beta(th0_data,th17_data);
    save(sprintf('result_data_simNB_%d.mat', idx),'result')
  end
end

% Plot ROC and PR curves for one simulation
t=zeros(1000,1);
t(901:1000)=1;

if 1
  load('result_data_simAR_1.mat')
  result = DyNB_MG_Beta(th0_data,th17_data);
  for i=1:1000; bf(i)=result{i}.log_bayes_factor; end
  prec_rec(bf, t)
end
if 0
  load('result_data_simNB_1.mat')
  result = DyNB_MG_Beta(th0_data,th17_data);
  for i=1:1000; bf(i)=result{i}.log_bayes_factor; end
  prec_rec(bf, t)
end
if 0
  load('result_data_simGP_1.mat')
  result = DyNB_MG_Beta(th0_data,th17_data);
  for i=1:1000; bf(i)=result{i}.log_bayes_factor; end
  prec_rec(bf, t)
end

