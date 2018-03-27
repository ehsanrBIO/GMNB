function y=get_NB_likelihood(X, r, p)

R = repmat(r',size(X,1),1);
log_numerator = gammaln(R+X) + X .* log(p) + R .* log(1-p);
log_denominator = gammaln(X+1) + gammaln(R);
log_pr = log_numerator - log_denominator;
y = sum(sum(log_pr));

end
