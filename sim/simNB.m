numsamples=4
T=5

for idx=1:20
  q1 = unifrnd(.8,1.2, numsamples, T); 
  p1 = q1 ./ (1+q1);
  q2 = unifrnd(.8,1.2, numsamples, T); 
  p2 = q2 ./ (1+q2);
  for k=1:900
    e0 = unifrnd(30, 50);
    truelambda = gamrnd(e0,10);
    ctrue = unifrnd(0.8, 2); 
    for i = 2:T 
      truelambda(i) = gamrnd(truelambda(i-1), 1/(ctrue));
    end 
    data1{k} = nbinrnd(repmat(truelambda,numsamples,1), (1-p1));
    
    data2{k} = nbinrnd(repmat(truelambda,numsamples,1), (1-p2));
  end
  for k=901:1000
    e0 = unifrnd(30, 50);
    truelambda = gamrnd(e0,10);
    ctrue = unifrnd(0.8, 2);
    for i = 2:T
      truelambda(i) = gamrnd(truelambda(i-1), 1/(ctrue));
    end

    data1{k} = nbinrnd(repmat(truelambda,numsamples,1), (1-p1));
   
    if (ctrue < 1) 
      ctrue = ctrue + .02;
    else
      ctrue = ctrue - .02;
    end

    for i = 2:T
      truelambda(i) = gamrnd(truelambda(i-1), 1/(ctrue));
    end
    data2{k} = nbinrnd(repmat(truelambda,numsamples,1), (1-p2));
    end

    save(sprintf('data_simNB_%d.mat', idx),'data1','data2')
  
    a1=cat(3, data1{:});
    x=[];
    for i = 1:5
      x=[x squeeze(a1(:, i, :))'];
    end
    dlmwrite(sprintf('simNB_sizeFactor_data1_%d.txt', idx), x);
    a1=cat(3, data2{:});
    x=[];
    for i = 1:5
      x=[x squeeze(a1(:, i, :))'];
    end
    dlmwrite(sprintf('simNB_sizeFactor_data2_%d.txt', idx), x);

end
