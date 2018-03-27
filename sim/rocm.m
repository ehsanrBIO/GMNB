t = zeros(1000,1);
t(901:1000)=1;
PR=[];
ROC = [];

if 1
for idx=1:20
  load(sprintf('result_data_simGP_%d.mat', idx))
  for i=1:1000; bf(i)=result{i}.log_bayes_factor;end

  [X,Y,T,AUCPR] = perfcurve(t,bf,1, 'xCrit', 'reca', 'yCrit', 'prec');
  [X,Y,T,AUC] = perfcurve(t,bf,1);

  PR=[PR, AUCPR];
  ROC=[ROC, AUC];
end
prmean=mean(PR)
prvar=sqrt(var(PR))

rocmean=mean(ROC)
rocvar=sqrt(var(ROC))
end

if 0
for idx=1:20
  load(sprintf('result_data_simNB_%d.mat', idx))
  for i=1:1000; bf(i)=result{i}.log_bayes_factor;end

  [X,Y,T,AUCPR] = perfcurve(t,bf,1, 'xCrit', 'reca', 'yCrit', 'prec');
  [X,Y,T,AUC] = perfcurve(t,bf,1);

  PR=[PR, AUCPR];
  ROC=[ROC, AUC];
end
prmean=mean(PR)
prvar=sqrt(var(PR))

rocmean=mean(ROC)
rocvar=sqrt(var(ROC))
end

if 0
for idx=1:20
  load(sprintf('result_data_simAR_%d.mat', idx))
  for i=1:1000; bf(i)=result{i}.log_bayes_factor;end

  [X,Y,T,AUCPR] = perfcurve(t,bf,1, 'xCrit', 'reca', 'yCrit', 'prec');
  [X,Y,T,AUC] = perfcurve(t,bf,1);

  PR=[PR, AUCPR];
  ROC=[ROC, AUC];
end
prmean=mean(PR)
prvar=sqrt(var(PR))

rocmean=mean(ROC)
rocvar=sqrt(var(ROC))
end
