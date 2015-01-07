function par_val = calc_PAPR(data)
% calc peak-to-average power ratio of input data (in vector form)


n = numel(data);

Re = real(data); 
Im = imag(data);
data = Re.*Re + Im.*Im;

par_val = n * max(data) / sum(data);



end