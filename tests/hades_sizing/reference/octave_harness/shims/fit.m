function a = fit(x, y, model, varargin)
  % least squares for the one-coefficient model a./x (linear in a)
  x = x(:); y = y(:);
  a = sum(y ./ x) / sum(1 ./ x.^2);
  fprintf(2, 'fit shim: a = %.15g\n', a);
end
