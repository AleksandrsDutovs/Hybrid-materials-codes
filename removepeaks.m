function rez=removepeaks(ind,w)
  k=[-w:w];
  h1=exp(-10*k.*k/w/w);
  h=h1/sum(h1);
  [L,W]=size(ind);
  nm=ind(:,1);
  intens=ind(:,2);
  filt=conv(intens,h,"same");
  diff=intens-filt;
  d2=diff.^2;
  max2=max(d2);
  m=1;
  for n=1:L
    if(d2(n)<max2)
      rez(m,:)=ind(n,1:2);
      m++;
    endif
  endfor
endfunction
