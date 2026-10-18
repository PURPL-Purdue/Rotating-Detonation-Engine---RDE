function sol = pdepe(m, pdefun, icfun, bcfun, x, t, varargin)
  x = x(:)'; nx = numel(x); xL = x(1:end-1); xR = x(2:end);
  if m == 0
    xi = (xL + xR)/2; zeta = xi;
  else
    xi = (xR - xL) ./ log(xR ./ xL); zeta = sqrt(xi .* (xR + xL)/2);
  end
  wR = (zeta.^(m+1) - xL.^(m+1))/(m+1); wL = (xR.^(m+1) - zeta.^(m+1))/(m+1);
  u0 = arrayfun(icfun, x)';
  opts = odeset('RelTol', 1e-3, 'AbsTol', 1e-6);
  [~, U] = ode15s(@(tt,u) rhs(tt, u), t, u0, opts);
  sol = U;
  function du = rhs(tt, u)
    u = u(:)'; duu = u(2:end) - u(1:end-1);
    if m == 0
      Ui = u(1:end-1) + duu .* (xi - xL) ./ (xR - xL); Ux = duu ./ (xR - xL);
    else
      lg = log(xR ./ xL); Ui = u(1:end-1) + duu .* log(xi ./ xL) ./ lg; Ux = duu ./ (xi .* lg);
    end
    c = zeros(1, nx-1); f = c; s = c;
    for i = 1:nx-1
      [c(i), f(i), s(i)] = pdefun(xi(i), tt, Ui(i), Ux(i));
    end
    [pl, ql, pr, qr] = bcfun(x(1), u(1), x(end), u(end), tt);
    fl = xi.^m .* f; du = zeros(nx, 1);
    du(1) = (fl(1) + x(1)^m * pl/ql + wR(1)*s(1)) / (wR(1)*c(1));
    du(2:nx-1) = (fl(2:end) - fl(1:end-1) + wL(1:end-1).*s(1:end-1) + wR(2:end).*s(2:end)) ./ (wL(1:end-1).*c(1:end-1) + wR(2:end).*c(2:end));
    du(nx) = (-x(end)^m * pr/qr - fl(end) + wL(end)*s(end)) / (wL(end)*c(end));
  end
end
