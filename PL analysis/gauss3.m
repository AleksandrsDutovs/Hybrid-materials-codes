function rez=gauss3(eV,p)
  rez=p(3)*exp(-((eV-p(1)).^2)/(2*p(2).^2)) ...
    + p(6)*exp(-((eV-p(4)).^2)/(2*p(5).^2)) ...
    + p(9)*exp(-((eV-p(7)).^2)/(2*p(8).^2));
endfunction