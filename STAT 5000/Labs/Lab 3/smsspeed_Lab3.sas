/*  Speed of text messages between teenagers and adults. */
data SMS;
	infile '~/my_shared_file_links/u63538023/STAT5000_Fall2024_ISU/smsspeed.csv' 
		dlm=',' firstobs=2;
	input Age AgeGroup $ Own Control;
run;

title1 'Default: 95% confidence interval';
proc ttest data=SMS; 
  class AgeGroup;
  var Own;
run;


title2 'Modified: 99% confidence interval';
proc ttest data=SMS alpha=0.01; 
  class AgeGroup;
  var Own;
run;