if 1
  load('data/DyNB_example_data.mat')
  result = DyNB_MG_Beta(th0_data,th17_data);
  save('result_GMNB_DyNBexample.mat', 'result')
end
if 0
  load('data/GSE52260_rc_human_timeSeriesRNAseq.mat')
  result = DyNB_MG_Beta(th0_data,th17_data);
  save('result_GMNB_GSE52260.mat', 'result')
end
