function useelement = findinavec(usecurinx,avec,i,j)

usecurinx2last=(usecurinx+4-j)/4;
avec_i=avec(i,:);
useelement=avec_i(usecurinx2last);

end

