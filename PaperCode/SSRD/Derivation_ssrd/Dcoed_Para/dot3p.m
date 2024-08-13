function s = dot3p(x,y,s)
   shi = double(single(s));
   slo = s - shi;
   for k = 1:length(x)
      xhi = double(single(x(k)));
      xlo = x(k) - xhi;
      yhi = double(single(y(k)));
      ylo = y(k) - yhi;
      tmp = xhi*yhi;
      zhi = double(single(tmp));
      zlo = tmp - zhi + xhi*ylo + xlo*yhi + xlo*ylo;

      tmp = shi + zhi;
      del = tmp - shi - zhi;
      shi = double(single(tmp));
      slo = tmp - shi + slo + zlo - del;
   end
   s = shi +  slo;
