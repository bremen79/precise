function [lcblist,ucblist] = precise_a_co96(x,delta)
%PRECISE_A_CO96    Portfolio REgret for Confidence SEquences with
%                  Approximation using Cover and Ordentlich [1996].
%   [L,U] = PRECISE_A_CO96(X,delta) produces two matrices, of the same 
%   dimension as X and with lower and upper confidence sequences with
%   probability of error delta. X is numer-of-samples by
%   number-of-repetitions.

%  This algorithm is described in  Orabona and Jun, "Tight Concentrations
%  and Confidence Sequences from the Regret of Universal Portfolio", ArXiv
%  2021.

N = size(x,1);
n_try = size(x,2);
dt = delta;
tol = 1e-4;

reg = @(t)   log(sqrt(pi)) + gammaln(t+1) - gammaln(t+.5);

n_algo = 2;

lcblist = zeros(N,n_try);
ucblist = zeros(N,n_try);
for i_try=1:n_try
  rdata = zeros(N,n_algo,2);
  runningmax = zeros(1,n_algo);
  runningmin = ones(1,n_algo);

  for t=1:N
    data = x(1:t,i_try); %Y(1:t);
    me = mean(data);
    va = var(data, 1);
    lcb = zeros(1,n_algo);
    ucb = ones(1,n_algo);
    i_algo = 0;
  
    %--- fan
    i_algo = i_algo + 1;
    rhs = reg(t) + log(1/dt);
    lb = 0;
    ub = me;
    lcbmaxfn = @(m,me,va,t) max(max_logwealth_fan3_lcb(m, me, va, t), max_logwealth_kl(m, me, va, t));
    if (lb == ub || lcbmaxfn(lb,me,va,t) - rhs <= 0)
      lcb(i_algo) = 0.0;
    else
      [lcb(i_algo), uu] = bsearch(@(m) lcbmaxfn(m,me,va,t) - rhs, lb, ub, tol);
    end
    lb = me;
    ub = 1;
    ucbmaxfn = @(m,me,va,t) max(max_logwealth_fan3_ucb(m, me, va, t), max_logwealth_kl(m, me, va, t));
    if (lb == ub || ucbmaxfn(ub, me, va, t) - rhs <= 0)
      ucb(i_algo) = 1.0;
    else
      [ll, ucb(i_algo)] = bsearch(@(m) ucbmaxfn(m,me,va,t) - rhs, lb, ub, tol);
    end

    runningmax(:) = max(runningmax, lcb);
    runningmin(:) = min(runningmin, ucb);
    rdata(t,:,1) = runningmax;
    rdata(t,:,2) = runningmin;
  end
  lcblist(:,i_try) = rdata(:,1,1);
  ucblist(:,i_try) = rdata(:,1,2);
end

end

function max_logwealth = max_logwealth_fan3_lcb(m, mu_hat, var_hat, t)
  if (m==0.0)
    max_logwealth = (.5*(mu_hat^2)/(var_hat + mu_hat^2))*t;
  elseif m == mu_hat
    max_logwealth = 0.0;
  else
    A = ((mu_hat - m)/m);
    B = ((var_hat + (mu_hat - m)^2)/m^2);
    lam = A/(A+B);
    max_logwealth = (A*A/(A+B) - (-log(1-lam) - lam)*B)*t;
  end
end

function max_logwealth = max_logwealth_fan3_ucb(m, mu_hat, var_hat, t)
  if (m==1.0)
    max_logwealth = (.5*((1 - mu_hat)^2)/(var_hat + (1 - mu_hat)^2))*t ;
  elseif m == mu_hat
    max_logwealth = 0.0;
  else
    A = ((m - mu_hat)/(1-m));
    B = ((var_hat + (m - mu_hat)^2)/(1-m)^2);
    lam = A/(A+B);
    max_logwealth = (A*A/(A+B) - (-log(1-lam) - lam)*B)*t;    
  end
end

function val = kl(p,q)
  if (p == 0)
    val = (1-p)*log((1-p)/(1-q));
  elseif (p == 1)
    val = p*log(p/q);
  else
    val = p*log(p/q) + (1-p)*log((1-p)/(1-q));
  end
end

function max_logwealth = max_logwealth_kl(m, mu_hat, var_hat, t)
  max_logwealth = t*kl(mu_hat,m);
end

function [lb,ub] = bsearch(fn, lb, ub, tol)
  if (~exist('tol','var'))
    tol = 1e-6;
  end

  fnlb = fn(lb);
  fnub = fn(ub);
  assert(fnlb * fnub ~= 0.0)
  if ~(( fnlb <= 0 && fnub >= 0 ) || ( fnlb >= 0 && fnub <= 0))
    keyboard;
  end
  assert(( fnlb <= 0 && fnub >= 0 ) || ( fnlb >= 0 && fnub <= 0))
  sign_fnlb = 1.0;
  if (fnlb <= 0 && fnub >= 0)
    sign_fnlb = -1.0;
  end
  max_iter = 1000;
  i_iter = 0;
  while ub-lb >= tol && i_iter <= max_iter
    mid = (lb+ub)/2;
    val = fn(mid);
    if (val*sign_fnlb > 0.0) % when the sign of val is equal to that of fn(lb)
      lb = mid;
    else
      ub = mid;
    end
    i_iter = i_iter + 1;
  end
  if (i_iter >= max_iter)
    fprintf("WARNING: max_iter has reached\n");
  end
end
