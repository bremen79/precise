function [lcblist,ucblist] = precise_co96(x,delta)
%PRECISE_CO96    Portfolio REgret for Confidence SEquences using Cover and
%                Ordentlich [1996].
%   [L,U] = PRECISE_CO96(X,delta) produces two matrices, of the same 
%   dimension as X and with lower and upper confidence sequences with
%   probability of error delta. X is numer-of-samples by
%   number-of-repetitions, where each sample must be in [0,1].

%  This algorithm is described in  Orabona and Jun, "Tight Concentrations
%  and Confidence Sequences from the Regret of Universal Portfolio", ArXiv
%  2021.

func = @(b,k,t) k*log(b+eps)+(t-k)*log(1-b+eps)+gammaln(t+1)+2*gammaln(1/2)-gammaln(k+1/2)-gammaln(t-k+1/2);

[T,repetitions]=size(x);

lcblist = zeros(T,repetitions);
ucblist = zeros(T,repetitions);

loginvdelta=log(1/delta);

for j=1:repetitions

    tic
    fprintf('Repetition %d\n',j);
    
    c = x(:,j);
    
    m_lb_old=eps;
    m_ub_old=1-eps;
    
     
    for i=1:length(c)
        mean_c=mean(c(1:i));
        
        % upper confidence interval
        m_ub = m_ub_old;
        m_lb = max(m_lb_old,mean_c);
                
        % calculate regret
        m_try = m_ub ;
        bmax=1/m_try;
        bmin=-1/(1-m_try);
        [bet_star, log_W_star] = find_max_log_wealth_constrained(c(1:i)-m_try,bmin,bmax);
        b=(-bmin+bet_star)/(bmax-bmin);
        bound=max(func(ceil(b*i-0.5)/i,ceil(b*i-0.5),i),func(floor(mean_c*i+0.5)/i,floor(mean_c*i+0.5),i));
        if log_W_star - bound >= loginvdelta
            while (m_ub - m_lb)>0.0001
                m_try = (m_ub + m_lb)/2;
                bmax=1/m_try;
                bmin=-1/(1-m_try);
                [bet_star, log_W_star] = find_max_log_wealth_constrained(c(1:i)-m_try,bmin,bmax);
                if log_W_star-bound >= loginvdelta
                    m_ub = m_try;
                    b=(-bmin+bet_star)/(bmax-bmin);
                    bound=max(func(ceil(b*i-0.5)/i,ceil(b*i-0.5),i),func(floor(mean_c*i+0.5)/i,floor(mean_c*i+0.5),i));
                    %if bound2> bound+0.0001
                    %    disp 'err';
                    %end
                    %bound=min(bound2, bound);
                else
                    m_lb = m_try;
                end
            end
        end
        ucblist(i,j) = m_ub;
        m_ub_old=m_ub;
        

        % lower confidence interval
        m_ub = min(m_ub,mean_c);
        m_lb = m_lb_old;
        
        % calculate regret
        m_try = m_lb;
        bmax=1/m_try;
        bmin=-1/(1-m_try);
        [bet_star, log_W_star] = find_max_log_wealth_constrained(c(1:i)-m_try,bmin,bmax);
        b=min((-bmin+bet_star)/(bmax-bmin),1);
        bound=max(func(floor(b*i+0.5)/i,floor(b*i+0.5),i),func(ceil(mean_c*i-0.5)/i,ceil(mean_c*i-0.5),i));
        if log_W_star - bound >= loginvdelta
            while (m_ub - m_lb)>0.0001
                m_try = (m_ub + m_lb)/2;
                bmax=1/m_try;
                bmin=-1/(1-m_try);
                [bet_star, log_W_star] = find_max_log_wealth_constrained(c(1:i)-m_try,bmin,bmax);
                if log_W_star-bound >= loginvdelta
                    m_lb = m_try;
                    b=min((-bmin+bet_star)/(bmax-bmin),1);
                    bound=max(func(floor(b*i+0.5)/i,floor(b*i+0.5),i),func(ceil(mean_c*i-0.5)/i,ceil(mean_c*i-0.5),i));
                    %if bound2> bound+0.0001
                    %    disp 'err';
                    %end
                    %bound=min(bound2, bound);
                else
                    m_ub = m_try;
                end
            end
        end
        lcblist(i,j) = m_lb;
        m_lb_old=m_lb;
    end
    toc
end

function [betstar, fval] = find_max_log_wealth_constrained(g,bmin,bmax)
% finds maximum wealth of max_{beta in [bmin,bmax]} prod_{i=1}^t (1+beta * g(i))

myf = @(bet) log(prod(1 + g.*bet));
df = @(bet) sum(g./(1 + g.*bet));
df2 = @(bet) -sum((g./(1 + g.*bet)).^2);

[betstar, fval] = newton_1d_bnd(myf, df, df2, max(bmin,-1e10), min(bmax,1e10));
