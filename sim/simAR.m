numsamples=4
T=5

for idx=1:20
  q1 = unifrnd(.8,1.2, numsamples, T);
  q2 = unifrnd(.8,1.2, numsamples, T);
  sigma = .05;
  for k=1:900
    K = gamrnd(2000,1);
    phi = unifrnd(.1, .9);
    betai=unifrnd(4.5,5.5);
    trueMu = normrnd(0, sigma/(1-phi^2));
    for i = 2:T
      trueMu(i) = normrnd(phi * trueMu(i-1), sigma);
    end
    mu = exp(trueMu+betai);
    data1{k} = nbinrnd(K, 1-(q1.*mu./(K + q1.*mu)));
    data2{k} = nbinrnd(K, 1-(q2.*mu./(K + q2.*mu)));
  end
  for k=901:1000
    K = gamrnd(2000,1);
    phi = unifrnd(.1, .9);
    betai=unifrnd(4.5,5.5);
    trueMu = normrnd(0, sigma/(1-phi^2));
    for i = 2:T
      trueMu(i) = normrnd(phi * trueMu(i-1), sigma);
    end
    mu = exp(trueMu +betai);
    data1{k} = nbinrnd(K, 1-(q1.*mu./(K + q1.*mu)));

    if (phi > .5) 
      phi = phi/1.5
    else
      phi = phi*1.5
    end 
    for i = 2:T
      trueMu(i) = normrnd(phi * trueMu(i-1), sigma);
    end
    mu = exp(trueMu+betai);
    data2{k} = nbinrnd(K, 1-(q2.*mu./(K + q2.*mu)));
  end
  save(sprintf('data_simAR_%d.mat', idx),'data1','data2')

  a1=cat(3, data1{:});
  x=[];
  for i = 1:5
    x=[x squeeze(a1(:, i, :))'];
  end
  dlmwrite(sprintf('simAR_sizeFactor_data1_%d.txt', idx), x);
  a1=cat(3, data2{:});
  x=[];
  for i = 1:5
    x=[x squeeze(a1(:, i, :))'];
  end
  dlmwrite(sprintf('simAR_sizeFactor_data2_%d.txt', idx), x);

end
