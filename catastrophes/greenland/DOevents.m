plot(age_ngrip, deltaO18_ngrip, 'b-'); hold on; plot(event_age, -40*ones(18,1), 'rx'); 

O18data1=deltaO18_ngrip(1:floor(event_age(1)/20));
O18data2=deltaO18_ngrip(floor(event_age(1)/20):floor(event_age(2)/20));
O18data3=deltaO18_ngrip(floor(event_age(2)/20):floor(event_age(3)/20));
O18data4=deltaO18_ngrip(floor(event_age(3)/20):floor(event_age(4)/20));
O18data5=deltaO18_ngrip(floor(event_age(4)/20):floor(event_age(5)/20));
O18data6=deltaO18_ngrip(floor(event_age(5)/20):floor(event_age(6)/20));
O18data7=deltaO18_ngrip(floor(event_age(6)/20):floor(event_age(7)/20));
O18data8=deltaO18_ngrip(floor(event_age(7)/20):floor(event_age(8)/20));
O18data9=deltaO18_ngrip(floor(event_age(8)/20):floor(event_age(9)/20));
O18data10=deltaO18_ngrip(floor(event_age(9)/20):floor(event_age(10)/20));
O18data11=deltaO18_ngrip(floor(event_age(10)/20):floor(event_age(11)/20));
O18data12=deltaO18_ngrip(floor(event_age(11)/20):floor(event_age(12)/20));
O18data13=deltaO18_ngrip(floor(event_age(12)/20):floor(event_age(13)/20));
O18data14=deltaO18_ngrip(floor(event_age(13)/20):floor(event_age(14)/20));
O18data15=deltaO18_ngrip(floor(event_age(14)/20):floor(event_age(15)/20));
O18data16=deltaO18_ngrip(floor(event_age(15)/20):floor(event_age(16)/20));
O18data17=deltaO18_ngrip(floor(event_age(16)/20):floor(event_age(17)/20));
O18data18=deltaO18_ngrip(floor(event_age(17)/20):floor(event_age(18)/20));
O18data19=deltaO18_ngrip(floor(event_age(18)/20):end);

%reverse the data  (I am sure there is a faster way to do this)
clear A*;clear B*;
A=deltaO18_ngrip;
for k=1:length(A), B(k)=A(length(A)-k+1); end
RdeltaO18_ngrip=B';

clear A*;clear B*;
A=O18data1; for k=1:length(A), B(k)=A(length(A)-k+1); end
R_O18data1=B';

clear A*;clear B*;
A=O18data2; for k=1:length(A), B(k)=A(length(A)-k+1); end
R_O18data2=B';

clear A*;clear B*;
A=O18data3; for k=1:length(A), B(k)=A(length(A)-k+1); end
R_O18data3=B';

clear A*;clear B*;
A=O18data4; for k=1:length(A), B(k)=A(length(A)-k+1); end
R_O18data4=B';

clear A*;clear B*;
A=O18data5; for k=1:length(A), B(k)=A(length(A)-k+1); end
R_O18data5=B';

clear A*;clear B*;
A=O18data6; for k=1:length(A), B(k)=A(length(A)-k+1); end
R_O18data6=B';

clear A*;clear B*;
A=O18data7; for k=1:length(A), B(k)=A(length(A)-k+1); end
R_O18data7=B';

clear A*;clear B*;
A=O18data8; for k=1:length(A), B(k)=A(length(A)-k+1); end
R_O18data8=B';

clear A*;clear B*;
A=O18data9; for k=1:length(A), B(k)=A(length(A)-k+1); end
R_O18data9=B';

clear A*;clear B*;
A=O18data10; for k=1:length(A), B(k)=A(length(A)-k+1); end
R_O18data10=B';

clear A*;clear B*;
A=O18data11; for k=1:length(A), B(k)=A(length(A)-k+1); end
R_O18data11=B';

clear A*;clear B*;
A=O18data12; for k=1:length(A), B(k)=A(length(A)-k+1); end
R_O18data12=B';

clear A*;clear B*;
A=O18data13; for k=1:length(A), B(k)=A(length(A)-k+1); end
R_O18data13=B';

clear A*;clear B*;
A=O18data14; for k=1:length(A), B(k)=A(length(A)-k+1); end
R_O18data14=B';

clear A*;clear B*;
A=O18data15; for k=1:length(A), B(k)=A(length(A)-k+1); end
R_O18data15=B';

clear A*;clear B*;
A=O18data16; for k=1:length(A), B(k)=A(length(A)-k+1); end
R_O18data16=B';

clear A*;clear B*;
A=O18data17; for k=1:length(A), B(k)=A(length(A)-k+1); end
R_O18data17=B';

clear A*;clear B*;
A=O18data18; for k=1:length(A), B(k)=A(length(A)-k+1); end
R_O18data18=B';

clear A*;clear B*;
A=O18data19; for k=1:length(A), B(k)=A(length(A)-k+1); end
R_O18data19=B';