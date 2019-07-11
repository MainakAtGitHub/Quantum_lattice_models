function f=fermi_func(x)

large=(x>50);
small=(x<-50);
neg=(x<0);
other=(0*x+1-large-small-neg>=1);
f=0*x;
f(large)=x(large)*0;
f(small)=1+x(small)*0;
f(neg)=1./(1+exp(x(neg)));
f(other)= exp(-x(other))./(1+exp(-x(other)));