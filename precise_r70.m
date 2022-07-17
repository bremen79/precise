function [lcblist,ucblist] = precise_r70(x,delta)
%PRECISE_R70    Portfolio REgret for Confidence SEquences using 
%               Robbins [1970].
%   [L,U] = PRECISE_R70(X,delta) produces two matrices, of the same 
%   dimension as X and with lower and upper confidence sequences with
%   probability of error delta. X is numer-of-samples by
%   number-of-repetitions.

%  This algorithm is described in  Orabona and Jun, "Tight Concentrations
%  and Confidence Sequences from the Regret of Universal Portfolio", ArXiv
%  2021.

[T,repetitions]=size(x);

lcblist = zeros(T,repetitions);
ucblist = zeros(T,repetitions);

loginvdelta=log(1/delta);
mn=inf;
mx=-inf;

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
        
        mn=min(c(i),mn);
        mx=max(c(i),mx); 
        
        % calculate regret
        m_try = m_ub ;
        [log_W_star] = find_max_log_wealth_constrained_lil(c(1:i)-m_try,mn-m_try,mx-m_try);
        if log_W_star >= loginvdelta
            while (m_ub - m_lb)>0.0001
                m_try = (m_ub + m_lb)/2;
                [log_W_star] = find_max_log_wealth_constrained_lil(c(1:i)-m_try,mn-m_try,mx-m_try);
                if log_W_star >= loginvdelta
                    m_ub = m_try;
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
        m_try = m_lb ;
        [log_W_star] = find_max_log_wealth_constrained_lil(c(1:i)-m_try,mn-m_try,mx-m_try);
        if log_W_star >= loginvdelta
            while (m_ub - m_lb)>0.0001
                m_try = (m_ub + m_lb)/2;
                [log_W_star] = find_max_log_wealth_constrained_lil(c(1:i)-m_try,mn-m_try,mx-m_try);
                if log_W_star >= loginvdelta
                    m_lb = m_try;
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

function [fval] = find_max_log_wealth_constrained_lil(g,mn,mx)
myf = @(bet) prod(1 + g.*bet);
df = @(bet) sum(g./(1 + g.*bet));
df2 = @(bet) -sum((g./(1 + g.*bet)).^2);

[betstar, fval] = newton_1d_bnd(myf, df, df2, -1, 1);

pdf = @(bet) log(log(6.6)+1)/(2*abs(bet)*(1+log(6.6/abs(bet)))*(log(1+log(6.6/abs(bet))))^2);

V=sum(g.^2);

if betstar>0
    s=mn;
else
    s=mx;
end

absbetstar=abs(betstar);
if absbetstar~=1
    delta=(1+min(s*betstar,0))/sqrt(V);
else
    delta=0;
end
delta=absbetstar-max(absbetstar-delta,0);
fval=log(max((fval-1)/(eps+log(fval))*abs(betstar),fval*exp(-1/2*delta^2/(1+min(s*betstar,0))^2*V)*delta)*pdf(absbetstar+eps));