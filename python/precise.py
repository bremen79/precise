from typing import Tuple
from math import lgamma as gammaln, exp, sqrt, pi, asin, pi, log as ln
from math import ceil, floor

import numpy as np

from opt_ctypes import find_mean_max_log_wealth_constrained

def regret_func(b,k,t):
  eps = np.finfo(float).eps  
  return k*ln(b+eps)+(t-k)*ln(1-b+eps)+gammaln(t+1)+2*gammaln(1/2)-gammaln(k+1/2)-gammaln(t-k+1/2)

def precise_co96(
    x: np.ndarray,
    delta: float,
) -> (np.array, np.array):
  """ PRECiSE-CO96
  
  Args:
    x: 1-D array of data points
    delta: error probability

  Returns:
    Tuple with two arrays (corresponding to LCB and UCB measurements) of length # data points.
  """

  prec = np.finfo(float).eps

  n = x.size
  lcblist = np.zeros(n)
  ucblist = np.zeros(n)

  m_lb_old=prec
  m_ub_old=1-prec

  for i in range(1, n+1):    
    loginvdelta = ln(1/delta) 

    mean_c = x[:i].mean()

    # upper confidence interval
    m_ub = m_ub_old
    m_lb = max(m_lb_old, mean_c)
    
    # calculate regret
    m_try = m_ub
    bmax = 1 / m_try
    bmin = -1 / (1-m_try)

    log_W_star, bet_star = find_mean_max_log_wealth_constrained(
      x[:i],
      bmin,
      bmax,
      m_try
    )
    b=(-bmin+bet_star)/(bmax-bmin);
    bound=max(regret_func(ceil(b*i-0.5)/i,ceil(b*i-0.5),i),regret_func(floor(mean_c*i+0.5)/i,floor(mean_c*i+0.5),i))
    #bound = ln(sqrt(pi)*gamma(i+1)/gamma(i+0.5)) # non-instance-dependent bound

    if log_W_star - bound >= loginvdelta:
      while (m_ub - m_lb) > 0.0001:
        m_try = (m_ub + m_lb)/2
        bmax = 1 / m_try
        bmin = -1 / (1-m_try)

        log_W_star, bet_star = find_mean_max_log_wealth_constrained(
          x[:i],
          bmin,
          bmax,
          m_try
        )
        if log_W_star - bound >= loginvdelta:
          m_ub = m_try
          b=(-bmin+bet_star)/(bmax-bmin)
          bound=max(regret_func(ceil(b*i-0.5)/i,ceil(b*i-0.5),i),regret_func(floor(mean_c*i+0.5)/i,floor(mean_c*i+0.5),i))
        else:
          m_lb = m_try

    ucblist[i-1] = m_ub
    m_ub_old = m_ub
    
    # lower confidence interval
    m_ub = min(m_ub, mean_c)
    m_lb = m_lb_old

    # calculate regret
    m_try = m_lb
    bmax = 1 / m_try
    bmin = -1 / (1-m_try)
    log_W_star, bet_star = find_mean_max_log_wealth_constrained(
        x[:i],
        bmin,
        bmax,
        m_try
    )
    b=min((-bmin+bet_star)/(bmax-bmin),1)
    bound=max(regret_func(floor(b*i+0.5)/i,floor(b*i+0.5),i),regret_func(ceil(mean_c*i-0.5)/i,ceil(mean_c*i-0.5),i))

    if log_W_star - bound >= loginvdelta:
      while (m_ub - m_lb) > 0.0001:
        m_try = (m_ub + m_lb)/2
        bmax = 1 / m_try
        bmin = -1 / (1-m_try)

        log_W_star, bet_star = find_mean_max_log_wealth_constrained(
          x[:i],
          bmin,
          bmax,
          m_try
        )
        if log_W_star-bound >= loginvdelta:
            m_lb = m_try
            b=min((-bmin+bet_star)/(bmax-bmin),1)
            bound=max(regret_func(floor(b*i+0.5)/i,floor(b*i+0.5),i),regret_func(ceil(mean_c*i-0.5)/i,ceil(mean_c*i-0.5),i))
        else:
            m_ub = m_try

    lcblist[i-1] = m_lb
    m_lb_old = m_lb

  return (lcblist, ucblist)

