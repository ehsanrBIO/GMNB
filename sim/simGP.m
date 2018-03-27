numsamples=4
T=5
covt =[
         0   -0.1667   -0.3333   -0.6667   -1.0000
    0.1667         0   -0.1667   -0.5000   -0.8333;
    0.3333    0.1667         0   -0.3333   -0.6667;
    0.6667    0.5000    0.3333         0   -0.3333;
    1.0000    0.8333    0.6667    0.3333         0];
coefficients = [0.3155; 1.2515; 0.0843];

for idx=1:20
  q1 = unifrnd(.8,1.2, numsamples, T);
  q2 = unifrnd(.8,1.2, numsamples, T);
  for k=1:900
    ma = 25000;
    while(ma > 20000)
      theta_1 = unifrnd(100, 10000);
      theta_2 = unifrnd(.5, 1);
      mu = unifrnd(1000, 2000, 1, T);
      sigma = theta_1 * exp(-covt.^2/(theta_2^2)) + theta_1*1e-12*eye(5);
      r = mvnrnd(mu,sigma);
      while (nnz(r<0) > 0)
        r = mvnrnd(mu,sigma);
      end
      X = [ones(size(r')) log10(r') log10(r').^2];
      disper = 2*10.^(X*coefficients);
      S2 = q1 .* r + (q1 .^2) .* disper';
      data1 = nbinrnd((q1 .* r).^2./(S2 - q1 .* r), 1-((S2 - q1 .* r)./S2));


      r = mvnrnd(mu,sigma);
      while (nnz(r<0) > 0)
        r = mvnrnd(mu,sigma);
      end
      X = [ones(size(r')) log10(r') log10(r').^2];
      disper = 2*10.^(X*coefficients);
      S2 = q2 .* r + (q2 .^2) .* disper';
      data2 = nbinrnd((q2 .* r).^2./(S2 - q2 .* r), 1-((S2 - q2 .* r)./S2));
   
      ma = max(max(data2(:)), max(data1(:))); 
    end
    res1{k} = data1;
    res2{k} = data2;
  end
  for k=901:1000
    ma = 25000;
    while( ma > 20000)
      theta_1 = unifrnd(100, 10000);
      theta_2 = unifrnd(.5, 1);
      mu = unifrnd(1000, 2000, 1, T);
      sigma = theta_1 * exp(-covt.^2/(theta_2^2)) + theta_1*1e-12*eye(5);
      rn = randn(1,5);
      r = mu+(rn*chol(sigma));
      while (nnz(r<0) > 0)
        rn = randn(1,5);
        r = mu+(rn*chol(sigma));
      end
      r = mvnrnd(mu,sigma);
      while (nnz(r<0) > 0)
        r = mvnrnd(mu,sigma);
      end
      X = [ones(size(r')) log10(r') log10(r').^2];
      disper = 2*10.^(X*coefficients);
      
      S2 = q1 .* r + (q1 .^2) .* disper';
      data1 = nbinrnd((q1 .* r).^2./(S2 - q1 .* r), 1-((S2 - q1 .* r)./S2));
      
      theta_1 =  10 * theta_1;
      if (theta_2 > .75)
        theta_2 = theta_2 -.25;
      else
        theta_2 = theta_2 +.25;
      end 
      sigma = theta_1 * exp(-covt.^2/(theta_2^2)) + theta_1*1e-12*eye(5);
      mu = 1.5*mu;
      r = mvnrnd(mu,sigma);
      while (nnz(r<0) > 0)
        r = mvnrnd(mu,sigma);
      end
      X = [ones(size(r')) log10(r') log10(r').^2];
      disper = 2*10.^(X*coefficients);
      S2 = q2 .* r + (q2 .^2) .* disper';
      data2 = nbinrnd((q2 .* r).^2./(S2 - q2 .* r), 1-((S2 - q2 .* r)./S2));
   
      ma = max(max(data2(:)), max(data1(:))); 
    end
    res1{k} = data1;
    res2{k} = data2;
  end
  data1=res1;
  data2=res2;
  save(sprintf('data_simGP_%d.mat', idx),'data1','data2')

  a1=cat(3, data1{:});
  x=[];
  for i = 1:5
    x=[x squeeze(a1(:, i, :))'];
  end
  dlmwrite(sprintf('simGP_sizeFactor_data1_%d.txt', idx), x);
  a1=cat(3, data2{:});
  x=[];
  for i = 1:5
    x=[x squeeze(a1(:, i, :))'];
  end
  dlmwrite(sprintf('simGP_sizeFactor_data2_%d.txt', idx), x);
end
