function [Y]=mixednormpdf(xx,mu1,mu2,sig1,sig2,mix)
%this evaluates the mixture-of-normals pdf of xx (this could be a vector)
%the parameters are obvious...mix is the probability of normal #1
Y1 = normpdf(xx,mu1,sig1);
Y2 = normpdf(xx,mu2,sig2);
Y=mix*Y1+(1-mix)*Y2;