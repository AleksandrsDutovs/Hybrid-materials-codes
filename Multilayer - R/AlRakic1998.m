function [nm nk]=AlRakic1998()
  data =load("-ascii",'AlRakic.dat');
  nm=data(:,1)*1000;
  nk=data(:,2)+j*data(:,3);
endfunction